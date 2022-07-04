#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

import base64
from collections import defaultdict
from contextlib import ExitStack
import csv
import itertools
import os
import pickle
import re
import shutil
import signal
import sys
from tempfile import mkdtemp
from tempfile import mkstemp
from tempfile import TemporaryDirectory
import traceback

from Bio import SeqIO
from Bio.Emboss.Applications import StretcherCommandline
from Bio.Seq import Seq
from Bio.SeqUtils.CheckSum import seguid
from more_itertools import consecutive_groups
from tqdm import tqdm

from . import __version__
from .dbm import sonarDBManager

# COMPATIBILITY
SUPPORTED_DB_VERSION = 3


class sonarTimeout:
    """
    this class is a helper class raising a TimeoutError within a defined context

    Example
    --------

    >>> with sonarTimeout(1) as _:
    ... 	sleep(60)
    Traceback (most recent call last):
    ...
    TimeoutError: Timeout

    Parameters
    ----------

    seconds : int
            define time in seconds until TimeoutError is raised
            values below 1 deactivate the TimeoutError
    error_message: str
            define error message to be shown [ 'Timeout' ]

    Attributes
    ----------

    seconds : int
            time in seconds when a TimeoutError is raised
    error_message : str [ 'Timeout' ]
            error message to be shown
    """

    def __init__(self, seconds, error_message="Timeout"):
        self.seconds = seconds
        self.error_message = error_message

    def __enter__(self):
        if self.seconds > 0:
            signal.signal(signal.SIGALRM, self.handle_timeout)
            signal.alarm(self.seconds)

    def __exit__(self, type, value, traceback):
        if self.seconds > 0:
            signal.alarm(0)

    def handle_timeout(self, signum, frame):
        raise TimeoutError(self.error_message)


class sonarFiler:
    """
    this class is a helper class providing a (temporary) file handler for
    writing in a given context

    Notes
    -----
            Please consider, that an existing file will be overwritten.

    Examples
    --------

    >>> with sonarFiler() as handle:
    ... 	fname = handle.name
    ... 	os.path.isfile(fname)
    True
    >>> os.path.isfile(fname)
    False

    Parameters
    ----------
    fname : str [ None ]
            define designated file name to open (mode=w). If None, a temporary file is
            created and deleted after use.

    Attributes
    ----------
    name : str
            stores the given file name
    basename : str
            stores file basename
    path : str
            stores absolute file path
    tmp : bool
            stores True if it is a temporary file else False
    handle: file handler
            opened file handler
    """

    def __init__(self, fname=None):
        self.fname = fname
        self.tmp = True if fname is None else False

    def __enter__(self):
        if self.fname is None:
            self.handle, path = mkstemp()
        else:
            self.handle = open(self.fname, "w")
            path = self.fname
        self.name = os.path.abspath(path)
        self.basename = os.path.basename(path)
        self.path = os.path.dirname(path)
        return self

    def __exit__(self, type, value, traceback):
        if self.tmp:
            os.remove(self.name)


class sonarCDS(object):
    """
    this object stores information about a coding sequence (CDS)

    Notes
    -----
            Please note, that genomic coordinates are processed and returned 0-based
            by this object. While start or single coordinates are inclusive,
            end coordinates of ranges are exclusive, expressed in a mathematical
            notation: [start, end)

    Examples
    --------

    Initiating an sonarCDS object:

    >>> cds = sonarCDS("Loc1", "ORF1", [(155, 170)], ["ATGTTATGAATGGCC"], "+")

    Accessing amino acid sequence (genetic code 1):

    >>> cds.aa
    'ML*MA'

    Accessing CDS coordinates or genome range:

    >>> cds.coords
    (155, 170)
    >>> cds.range
    range(155, 170)

    Parameters
    ----------
    locus: str [ None ]
            define the gene locus accession
    symbol : str
            define the gene/protein symbol
            (e.g. ORF1b)
    coords : int
            define a sorted list of (start, end) tuples describing all exons of the
            gene (coordinates are 0-based, starts inclusive, ends exclusive,
            start always lower than end)
    seqs : list
            define a sorted list of exon nucleotide sequences (alswas forward strand
            sequence)
    strand : {'+', '-'}
            define the genomic strand the gene is encoded on (+ or -)
    translation_table : int [ 1 ]
            define the genetic code table used for in silico translation (see
            https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)

    Attributes
    ----------
    locus: str [ None ]
            sores the gene locus accession
    symbol : str
            stores the gene/protein symbol
    start : int
            stores the genomic start coordinate (0-based, inclusive). The start
            coordinate is always lower than the end coordinate.
    end : int
            stores the genomic end coordinate (0-based, exclusive). The end
            coordinate is always greater than the start coordinate.
    strand : str
            stores the coding strand as '+' or '-' of the respective CDS
    seqs : list
            stores ordered list of exon sequences of the respective CDS
    coords : tuple
            stores a tuple of CDS start and end coordinate
    coordlist : list
            stores a list of (start, end) tuples of of all exons (coordinates are
            0-based, starts inclusive, ends exclusive, start always lower than end)
    range : range
            stores a range from CDS start to end
    ranges : list
            stores a list of exon ranges
    nuc : str
            stores the coding nucleotide sequence (joined exons)
    aa : str
            stores the in silico translated amino acid sequence
    coding_positions: list
            stores all genomic positions part of the coding sequences of this gene
            as ordered list (positions might be redundant in case of ribosomal slippage)
    coding_positions_set: set
            stores all genomic positions part of the coding sequences of this gene
            as set
    translation_table : int [ 1 ]
            stores the genetic code table used for in silico translation (see
            https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
    length : int
            stores the aa length of the CDS product
    """

    def __init__(self, locus, symbol, coords, seqs, strand, translation_table=1):
        if not symbol:
            sys.exit("cds error: protein symbol missing for locus " + locus)
        if ":" in symbol:
            sys.exit("cds error: protein symbol " + symbol + " contains ':'.")
        self.symbol = symbol
        self.locus = locus
        self.start = coords[0][0]  # inclusive
        self.end = coords[-1][1]  # exclusive
        self.strand = strand
        self.seqs = seqs
        self.coordlist = coords
        self.ranges = [range(s, e) for s, e in coords]
        self.translation_table = translation_table
        self.__aa = None
        self.__nuc = None
        self.__coding_positions = None
        self.__coding_positions_set = None
        self.__length = None

    @property
    def nuc(self):
        if self.__nuc is None:
            self.__nuc = "".join(self.seqs)
        return self.__nuc

    @property
    def aa(self):
        if self.__aa is None:
            nuc = (
                self.nuc
                if self.strand == "+"
                else str(Seq(self.nuc).reverse_complement())
            )
            l = len(nuc)
            if l % 3 == 1:
                l = -1
            elif l % 3 == 2:
                l = -2
            self.__aa = str(Seq.translate(nuc[:l], table=self.translation_table))
        return self.__aa

    @property
    def coords(self):
        return self.start, self.end

    @property
    def range(self):
        return range(self.start, self.end)

    @property
    def coding_positions(self):
        if self.__coding_positions is None:
            self.__coding_positions = []
            for r in self.ranges:
                self.__coding_positions.extend(list(r))
        return self.__coding_positions

    @property
    def length(self):
        if self.__length is None:
            self.__length = len(self.__coding_positions) / 3
        return self.__length

    @property
    def coding_positions_set(self):
        if self.__coding_positions_set is None:
            self.__coding_positions_set = set(self.coding_positions)
        return self.__coding_positions_set

    def aa_to_nuc_pos(self, x):
        """
        function to return the respective lower-bound genome position for
        a given protein position

        Examples
        --------

        >>> cds=sonarCDS("loc1", "prot1", [(10, 24)], ['ATGTTATCCTGAAA'], "+")
        >>> cds.aa_to_nuc_pos(0)
        10
        >>> cds.aa_to_nuc_pos(2)
        16

        Parameters
        ----------
        x : int
                define the protein position to convert (0-based, inclusive)

        Returns
        -------
        int
                respective lower-bound genome position
        """
        return self.coding_positions[3 * x]

    def iter_coords(self, codon_wise=False):
        """
        function to iterate over genomic coordinate of all exons of an
        annotated coding sequence (CDS).

        Examples
        --------

        >>> cds=sonarCDS("loc1", "prot1", [(10, 15), (14, 18)], ['ATGTG', 'CTAATGA'], "+")
        >>> for i in cds.iter_coords():
        ... 	print(i)
        10
        11
        12
        13
        14
        14
        15
        16
        17
        >>> for i in cds.iter_coords(True):
        ... 	print(i)
        10
        13
        15

        Parameters
        ----------
        codon_wise : bool
                if true, iterator considers only the lower-bound coordinates of each codon

        Returns
        -------
        iterator
                iterator of genomic coordinates

        """
        j = (
            self.coding_positions
            if not codon_wise
            else zip(
                self.coding_positions[0::3],
                self.coding_positions[1::3],
                self.coding_positions[2::3],
            )
        )
        for i in j:
            yield i

    def is_exon(self, x, y=None):
        """
        function to check if a given genomic coordinate (range) overlaps with the
        coding part of this coding sequence (CDS).

        Examples
        --------

        >>> cds=sonarCDS("loc1", "prot1", [(10, 15), (25, 32)], ['ATGTG', 'CTAATGA'], "+")
        >>> cds.is_exon(10)
        True
        >>> cds.is_exon(15)
        False

        Parameters
        ----------
        x : int
                genomic (start) coordinate (0-based, inclusive)
        y : int [ None ]
                genomic end coordniate (0-based, exclusive)

        Returns
        -------
        bool
                True if coordinate(s) are within or overlapping with this exons of
                the CDS, False otherwise.

        Dev Note
        --------
        Working with intersection of ranges is for long sequences less performant
        why we don't use it here.

        """
        if y is None:
            y = x + 1
        for start, end in self.coordlist:
            if y >= start and end >= x:
                return True
        return False

    def is_cds(self, x, y=None):
        """
        function to check if a given genomic coordinate (range) overlaps with
        this coding sequence (CDS).

        Examples
        --------

        >>> gff=sonarCDS("loc1", "prot1", [(10, 15), (25, 32)], ['ATGTG', 'CTAATGA'], "+")
        >>> gff.is_cds(28)
        True
        >>> gff.is_cds(15)
        False

        Parameters
        ----------
        x : int
                genomic (start) coordinate (0-based, inclusive)
        y : int [ None ]
                genomic end coordinate (0-based, exclusive, greater than x)

        Returns
        -------
        bool
                True if coordinate(s) are within or overlapping with this CDS, False otherwise.

        Dev Note
        --------
        Working with intersection of ranges is for long sequences less performant

        """
        if y is None:
            y = x + 1
        return y >= self.start and self.end >= x

    def is_frameshift_del(self, x, y):
        """
        function to check if a deletion of a given genomic range (x to y)
        leads to an frameshift within this CDS.

        Examples
        --------

        >>> cds=sonarCDS("loc1", "prot1", [(10, 15), (25, 32)], ['ATGTG', 'CTAATGA'], "+")
        >>> cds.is_frameshift_del(11, 13)
        True
        >>> cds.is_frameshift_del(14,16)
        True
        >>> cds.is_frameshift_del(15,17)
        False
        >>> cds.is_frameshift_del(27,30)
        False

        >>> cds=sonarCDS("loc1", "prot1", [(10, 15), (15, 16), (15,20)], ['ATGTG', 'G', 'GATC'], "+")
        >>> cds.is_frameshift_del(15, 16)
        False
        >>> cds.is_frameshift_del(13, 16)
        True

        Parameters
        ----------
        x : int
                genomic (start) coordinate (0-based, inclusive)
        y : int
                genomic end coordinate (0-based, exclusive, greater than x)

        Returns
        -------
        bool
                True if deletion of the given genomic range causes a frameshift mutation
                within the CDS, False otherwise.
        """
        if (
            self.is_cds(x, y)
            and len([True for z in self.coding_positions if z < x or z >= y]) % 3 != 0
        ):
            return True
        return False

    def is_frameshift_in(self, x, l):
        """
        function to check if a insertion at given genomic position (x) with a given
        insertion length (l) leads to an frameshift within this CDS.

        Examples
        --------

        >>> cds=sonarCDS("loc1", "prot1", [(10, 16), (15, 21)], ['ATGTGC', 'GATNTC'], "+")
        >>> cds.is_frameshift_in(12, 3)
        False
        >>> cds.is_frameshift_in(12, 7)
        True
        >>> cds.is_frameshift_in(15, 4)
        True
        >>> cds.is_frameshift_in(15, 3)
        False

        Parameters
        ----------
        x : int
                genomic (start) coordinate (0-based, inclusive)
        l : int
                length of insertion (excluding anchor base)

        Returns
        -------
        bool
                True if deletion of the given genomic range causes a frameshift mutation
                within the CDS, False otherwise.
        """
        if l % 3 != 0 and x in self.coding_positions_set:
            return True
        return False


class sonarGFF(object):
    """
    this object stores CDS objects based on a GFF3 file.

    Notes
    -----
            Please note, that genomic coordinates are processed and returned 0-based
            by this object. While start or single coordinates are inclusive,
            end coordinates of ranges are exclusive, expressed in a mathematical
            notation: [start, end)

            Please note, that only single molecule genome annotations can be handled
            by this object.

    Examples
    --------

    Initiating an sonarGFF object. In this example the REF_GFF_FILE and REF_FASTA_FILE
    variable stores the path of an GFF3 and FASTA file containing the annotation
    and genomic sequence of the SARS-COV-2 NC_045512.2, respectively.

    >>> gff = sonarGFF(REF_GFF_FILE, REF_FASTA_FILE)

    Parameters
    ----------
    gff3 : str
            define a path to a valid GFF3 file storing genome annotation
    fna : str
            define a path to a valid FASTA file storing the nucleotide
            sequence of the annotated genome
    translation_table : int [ 1 ]
            define the genetic code table used for in silico translation of CDS (see
            https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)

    Attributes
    ----------
    translation_table : int
            stores the genetic code table used for in silico translation of CDS (see
            https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
    cds : list
            stores a list of sonarCDS objects (one per CDS of the given annotation)
    coords : dict
            stores a dictionary with protein symbol as keys and respective 0-based
            genomic coordinate tuples (start always lower than end coordinate,
            start coordinate inclusive, end coordinate exclusive)
    ranges : dict
            stores a dictionary with protein symbol as keys and respective list of
            coding ranges
    cds_positions : set
            stores a set of all genomic positions annotated within a coding gene
            (includes exons and introns).
    exon_positions : set
            stores a set of all genomic positions within annotated exons.
    symbols : list
            stores a list of protein symbols

    """

    def __init__(self, gff3, fna, translation_table=1):
        self.translation_table = translation_table
        self.cds = self.process_gff3(gff3, fna)
        self.coords = {x.symbol: (x.start, x.end) for x in self.cds}
        self.ranges = {x.symbol: x.ranges for x in self.cds}
        self.symbols = [x.symbol for x in self.cds]
        self.__cds_positions = None
        self.__exon_positions = None

    @property
    def cds_positions(self):
        if self.__cds_positions is None:
            positions = set()
            for ranges in self.ranges.values():
                for r in ranges:
                    positions.update(r)
            self.__cds_positions = positions
        return self.__cds_positions

    @property
    def exon_positions(self):
        if self.__exon_positions is None:
            positions = set()
            for ranges in self.ranges.values():
                for r in ranges:
                    positions.update(r)
            self.__exon_positions = positions
        return self.__exon_positions

    def in_any_exon(self, x, y=None):
        """
        function to check if a given genomic coordinate (range) overlaps with the
        coding region of any annotated coding sequence (CDS).

        Examples
        --------

        >>> gff=sonarGFF(REF_GFF_FILE, REF_FASTA_FILE)
        >>> gff.in_any_exon(21562)
        True
        >>> gff.in_any_exon(25384)
        False
        >>> gff.in_any_exon(25380, 25384)
        True

        Parameters
        ----------
        x : int
                genomic (start) coordinate (0-based, inclusive)
        y : int [ None ]
                genomic end coordniate (0-based, exclusive, greater than x)

        Returns
        -------
        bool
                True if coordinate(s) within CDS, False otherwise.

        """
        if y is None:
            y = x + 1
        for cds in self.cds:
            if cds.is_exon(x, y):
                return True
        return False

    def in_any_cds(self, x, y=None):
        """
        function to check if a given genomic coordinate (range) overlaps with an
        annotated coding sequence (CDS).

        Examples
        --------

        >>> gff=sonarGFF(REF_GFF_FILE, REF_FASTA_FILE)
        >>> gff.in_any_cds(21562)
        True
        >>> gff.in_any_cds(25384)
        False
        >>> gff.in_any_cds(25380, 25384)
        True

        Parameters
        ----------
        x : int
                genomic (start) coordinate (0-based, inclusive)
        y : int
                genomic end coordniate (0-based, exclusive) [ None ]

        Returns
        -------
        bool
                True if coordinate(s) within CDS, False otherwise.

        """
        if y is None:
            y = x + 1
        for cds in self.cds:
            if cds.is_cds(x, y):
                return True
        return False

    def any_frameshift_by_deletion(self, start, end):
        """
        function to check if a deletion causes a frameshift in any annotated
        CDS.

        Returns
        -------
        bool
                True if dna variant causes a frameshift, otherwise False.
        """
        for cds in self.cds:
            if cds.is_frameshift_del(start, end):
                return True
        return False

    def any_frameshift_by_insertion(self, pos, l):
        """
        function to check if a insertion causes a frameshift in any annotated
        CDS.

        Returns
        -------
        bool
                True if dna variant causes a frameshift, otherwise False.
        """
        if (l - 1) % 3 != 0 and self.in_any_exon(pos):
            return True
        else:
            return False

    def process_gff3(self, gff, fna):  # noqa: C901
        """
        function to parse CDS from a given GFF3 file

        Examples
        --------

        >>> os.chdir(os.path.dirname(os.path.realpath(__file__)))
        >>> gff = sonarGFF(REF_GFF_FILE, REF_FASTA_FILE)
        >>> gff.coords == {'ORF1a': (265, 13483), 'ORF1b': (265, 21555), 'S': (21562, 25384), 'ORF3a': (25392, 26220), 'E': (26244, 26472), 'M': (26522, 27191), 'ORF6': (27201, 27387), 'ORF7a': (27393, 27759), 'ORF7b': (27755, 27887), 'ORF8': (27893, 28259),'N': (28273, 29533), 'ORF10': (29557, 29674)}
        True

        Parameters
        ----------
        gff : str
                path to a valid GFF3 file storing the genome annotation
        fna : str
                path to a valid FASTA file storing the respective genome sequence

        Returns
        -------
        list
                list of sonarCDS objects for CDS annotated in the given GFF3 file
                sorted by CDS start (lower genomic coordinate).

        """

        symbol_regex = re.compile("gene=([^;]+)(?:;|$)")
        locus_regex = re.compile("locus_tag=([^;]+)(?:;|$)")
        id_regex = re.compile("ID=([^;]+)(?:;|$)")

        record = SeqIO.read(fna, "fasta")
        gseq = str(record.seq).upper()

        with open(gff, "r") as handle:
            cds = {}
            for line in handle:
                fields = line.rstrip("\r\n").split("\t")
                if line.startswith("#") or len(fields) < 7:
                    continue
                if fields[2] == "CDS":
                    id = id_regex.search(fields[-1]).groups(1)[0]
                    symbol = symbol_regex.search(fields[-1]).groups(1)[0]
                    locus = locus_regex.search(fields[-1]).groups(1)[0]
                    strand = fields[6]
                    s = int(fields[3]) - 1
                    e = int(fields[4])
                    if id not in cds:
                        cds[id] = {
                            "locus": locus,
                            "symbol": symbol,
                            "coords": [(s, e)],
                            "strand": strand,
                        }
                    elif id in cds:
                        if symbol != cds[id]["symbol"]:
                            sys.exit("gff error: multiple symbols for locus " + locus)
                        if strand != cds[id]["strand"]:
                            sys.exit("gff error: different strands for locus " + locus)
                        cds[id]["coords"].append((s, e))

        cdsobjs = []
        for locus, data in cds.items():
            seqs = []
            for s, e in data["coords"]:
                if data["strand"] == "+":
                    seqs.append(gseq[s:e])
                else:
                    seqs.append(str(Seq.reverse_complement(gseq[s:e])))
            cdsobjs.append(
                sonarCDS(
                    data["locus"],
                    data["symbol"],
                    data["coords"],
                    seqs,
                    data["strand"],
                    self.translation_table,
                )
            )

        return sorted(cdsobjs, key=lambda x: x.start)


class sonarALIGN(object):
    """
    this object performs a pairwise sequence alignment and provides/stores selected
    alignment functionalities/statistics.

    Notes
    -----
            Please note, that genomic coordinates are processed and returned 0-based
            by this object. While start or single coordinates are inclusive,
            end coordinates are exclusive, expressed as mathematical notation:
            [start, end)

            Please note, alignment is based on EMBOSS Stretcher.

    Example
    -------

    In this example the QRY_FASTA_FILE and REF_FASTA_FILE variables store
    the path of FASTA files containing the query and reference genome sequences,
    respectively.

    >>> algn = sonarALIGN(QRY_FASTA_FILE, REF_FASTA_FILE)

    Parameters
    ----------
    query_file : str
            define a path to a valid FASTA file storing the query genome sequence
    target_file : str
            define a path to a valid FASTA file storing the target genome sequence
            (= reference)
    out_file : str [ None ]
            define a path to an output file that will store the FASTA formatted
            alignment. Please consider, that an existing file will be overwritten!
            If None, a temporary file is used and deleted after processing.
    sonarGFFObj : object [ None ]
            define a sonarGFF object based on the reference genome annotation

    Attributes
    ----------
    aligned_query : str
            stores the aligned upper-case query sequence (U replaced by T)
    aligned_target : str
            stores the aligned upper-case target or reference sequence (U replaced by T)
    gff : object
            stores the sonarGFF object if provided, otherwise None
    dnadiff : list
            stores a list of tuples for each genomic variation (based on the alignment).
            Each tuple consists of:
             - reference base (or bases in case of deletions)
             - query base (or bases in case of insertions)
             - genomic coordinate (0-based, inclusive)
             - genomic end coordinate (in case of InDels 0-based and exlusive otherwise None)
             - None
             - None
            Accordingly to the VCF format, insertions are expressed considering the upstream
            base as anchor. As a special case, an insertion at the start of the sequence
            has no anchor and a genomic coordinate of -1. Deletions are expressed for each
            position they occur and not fused. The last two tuple elements
            are always None to keep the length according to tuples stored in aadiff.
    aadiff : list
            stores a list of tuples for each amino acid variation in an annotated protein.
            Each tuple consists of:
             - reference amino acid (or amino acids in case of deletions)
             - query amino acid (or amino acids in case of insertions)
             - protein position (0-based, inclusive)
             - protein end position (in case of InDels 0-based and exlusive otherwise None)
             - protein symbol
             - gene locus
            Accordingly to the VCF format, insertions are expressed considering the upstream
            base as anchor. As a special case, an insertion at the start of the sequence
            has no anchor and a genomic coordinate of -1. Deletions are expressed for each
            position they occur and not fused. The last two tuple elements
            are always None to keep the length according to tuples stored in aadiff.
    """

    def __init__(self, query_file, target_file, out_file=None, sonarGFFObj=None):
        self.aligned_query, self.aligned_target = self.align_dna(
            query_file, target_file, out_file
        )
        self.gff = sonarGFFObj
        self._insert_regex = re.compile("[^-]-+")
        self._gap_regex = re.compile("-+")
        self._codon_regex = re.compile("[^-]-*[^-]-*[^-]-*")
        self._leading_gap_regex = re.compile("^-+")
        self._tailing_gap_regex = re.compile("-+$")
        self._dnadiff = None
        self._aadiff = None
        self.__target_coords_matrix = None

    @property
    def _target_coords_matrix(self):
        if self.__target_coords_matrix is None:
            self.__target_coords_matrix = [
                len(x.group()) for x in re.finditer(".-*", self.aligned_target)
            ]
        return self.__target_coords_matrix

    def use_stretcher(
        self,
        query_file,
        target_file,
        out_file,
        gapopen=16,
        gapextend=4,
        right_align=True,
    ):
        """
        function to perform a pairwise aligment using EMBOSS Stretcher

        Parameters
        ----------
        query_file : str
                define a path to a valid FASTA file storing the query sequence
        target_file : str
                define a path to a valid FASTA file storing the target sequence
                (= reference)
        out_file : str
                define a path to a file that will store the alignment. Please consider,
                that an existing file will be overwritten.
        gapopen : int [ 16 ]
                define penalty for gap opening
        gapextend : int [ 4 ]
                define penalty for gap extension

        Returns
        -------
        list
          list of aligned query and target sequence, in that order
        """
        temp = True if not out_file else False
        if temp:
            handle, out_file = mkstemp()
        cline = StretcherCommandline(
            asequence=query_file,
            bsequence=target_file,
            gapopen=gapopen,
            gapextend=gapextend,
            outfile=out_file,
            aformat="fasta",
        )
        stdout, stderr = cline()
        alignment = [str(x.seq) for x in SeqIO.parse(out_file, "fasta")]
        if temp:
            os.remove(out_file)
        if right_align:
            alignment = self.left_align_gaps(*alignment)
        return alignment

    def left_align_gaps(self, query, target):
        """
        function to align gaps to the left in two aligned sequences

        Parameters
        ----------
        query : str
                define the query sequence in aligned form
        target : str
                define the target sequence (reference) in aligned form

        Returns
        -------
        list
          aligned query and target sequence strings with left-aligned gaps,
          in that order.
        """
        l = len(query) - 1
        for match in re.finditer("-+", query):
            s = match.start() - 1
            e = match.end() - 1
            g = "-" * (e - s)
            while s >= 0 and e < l and query[s] == target[e]:
                query = query[:s] + g + query[s] + query[e + 1 :]
                s -= 1
                e -= 1
        for match in re.finditer("-+", target):
            s = match.start() - 1
            e = match.end() - 1
            g = "-" * (e - s)
            while s >= 0 and e < l and target[s] == query[e]:
                target = target[:s] + g + target[s] + target[e + 1 :]
                s -= 1
                e -= 1
        return query, target

    def align_dna(
        self,
        query_file,
        target_file,
        out_file=None,
        gapopen=16,
        gapextend=4,
        right_align=True,
    ):
        """
        function to perform the default pairwise nucleotide aligment

        Parameters
        ----------
        query_file : str
                define a path to a valid FASTA file storing the query sequence
        target_file : str
                define a path to a valid FASTA file storing the target sequence
                (= reference)
        out_file : str
                define a path to a file that will store the alignment. Please consider,
                that an existing file will be overwritten.
        gapopen : int [ 16 ]
                define penalty for gap opening
        gapextend : int [ 4 ]
                define penalty for gap extension

        Returns
        -------
        list
          list of aligned query and target sequence
        """
        return self.use_stretcher(
            query_file, target_file, out_file, gapopen, gapextend, right_align
        )

    def real_pos(self, x):
        """
        function to convert an alignment position to the position in the
        unaligned target sequence (= reference).

        Example
        -------
        In this example the QRY_FASTA_FILE and REF_FASTA_FILE variables store
        the path of FASTA files containing the query and reference genome sequences,
        respectively.

        >>> algn = sonarALIGN(QRY_FASTA_FILE, REF_FASTA_FILE)
        >>> algn.real_pos(29282)
        29282

        Parameters
        ----------
        x : int
                define a position within the alignment (0-based)

        Returns
        -------
        int
          corresponding position (0-based) in the unaligned target/reference
          sequence
        """
        return x - self.aligned_target[: x + 1].count("-")

    def align_pos(self, x):
        """
        function to convert an target/reference position to the corresponding
        position in the alignment.

        Example
        -------

        >>> algn = sonarALIGN(QRY_FASTA_FILE, REF_FASTA_FILE)
        >>> algn.align_pos(29282)
        29282

        Parameters
        ----------
        x : int
                define a reference position (0-based)

        Returns
        -------
        int
          corresponding position of the sequence alignment
        """
        return sum(self._target_coords_matrix[:x])

    def iter_vars(self):  # noqa: C901
        """
        function to iterate variations on nucleotide level.

        Example
        -------

        In this example the QRY_FASTA_FILE and REF_FASTA_FILE variables store
        the path of FASTA files containing the query and reference genome sequences,
        respectively. The reference is NC_045512.2 while the query is a B.1.1.7
        prototype sequence.

        >>> algn = sonarALIGN(QRY_FASTA_FILE, REF_FASTA_FILE)
        >>> for x in algn.iter_dna_vars():
        ... 	print(x)
        ('C', 'T', 3266, None, None, None)
        ('C', 'A', 5387, None, None, None)
        ('T', 'C', 6953, None, None, None)
        ('T', '', 11287, None, None, None)
        ('C', '', 11288, None, None, None)
        ('T', '', 11289, None, None, None)
        ('G', '', 11290, None, None, None)
        ('G', '', 11291, None, None, None)
        ('T', '', 11292, None, None, None)
        ('T', '', 11293, None, None, None)
        ('T', '', 11294, None, None, None)
        ('T', '', 11295, None, None, None)
        ('T', '', 21764, None, None, None)
        ('A', '', 21765, None, None, None)
        ('C', '', 21766, None, None, None)
        ('A', '', 21767, None, None, None)
        ('T', '', 21768, None, None, None)
        ('G', '', 21769, None, None, None)
        ('T', '', 21990, None, None, None)
        ('T', '', 21991, None, None, None)
        ('A', '', 21992, None, None, None)
        ('A', 'T', 23062, None, None, None)
        ('C', 'A', 23270, None, None, None)
        ('C', 'A', 23603, None, None, None)
        ('C', 'T', 23708, None, None, None)
        ('T', 'G', 24505, None, None, None)
        ('G', 'C', 24913, None, None, None)
        ('C', 'T', 27971, None, None, None)
        ('G', 'T', 28047, None, None, None)
        ('A', 'G', 28110, None, None, None)
        ('G', 'C', 28279, None, None, None)
        ('A', 'T', 28280, None, None, None)
        ('T', 'A', 28281, None, None, None)
        ('C', 'T', 28976, None, None, None)

        Returns
        -------
        iterator of tuples
                each tuple represents a nucleotide level variation and consists of:
                         - target nucleotide
                         - query nucleotide(s)
                         - target or reference start position (0-based
                         - target or reference end position (0-based)
                         - None
                         - None
                Accordingly to the VCF format, insertions are expressed considering the upstream
                base as anchor. As a special case, an insertion at the start of the sequence
                has no anchor and a genomic coordinate of -1. Deletions are are expressed for
                each position they occur and not fused. The last two tuple elements
                are always None to keep the length according to tuples stored in aadiff.
        """
        target = self.aligned_target
        query = self.aligned_query

        # deletions
        for match in self._gap_regex.finditer(query):

            # nucleotide level
            s = self.real_pos(match.start())
            e = self.real_pos(match.end())
            ref = target[match.start() : match.end()]
            fs = 1 if self.gff and self.gff.any_frameshift_by_deletion(s, e) else 0
            yield "", ref, "", self.real_pos(s), self.real_pos(e), fs

            # protein level
            if self.gff:
                for cds in self.gff.cds:
                    if cds.is_cds(s, e):
                        dels = []
                        for j, codon_coords in enumerate(cds.iter_coords(True), 0):
                            if [x for x in codon_coords if x >= s and x <= e]:
                                dels.append(j)

                        for group in consecutive_groups(dels):
                            group = list(group)
                            yield cds.symbol, cds.aa[
                                group[0] : group[-1] + 1
                            ], "", group[0], group[-1], 0

        # insertions
        for match in self._gap_regex.finditer(target):

            # nucleotide level
            s = self.real_pos(match.start())
            e = s + 1
            ref = "" if s < 0 else target[match.start()]
            alt = query[match.start() : match.end()]
            fs = (
                1
                if self.gff and self.gff.any_frameshift_by_insertion(s, len(alt))
                else 0
            )
            yield "", ref, alt, s, e, fs

            # protein level
            if self.gff:
                for cds in self.gff.cds:
                    if cds.is_cds(s, e):
                        for j, codon_coords in enumerate(cds.iter_coords(True), 0):
                            if s in codon_coords:
                                ref = cds.aa[j]
                                alt = self.translate(alt, cds.translation_table)
                                if ref != alt:
                                    yield cds.symbol, cds.aa[j], alt, j, j + 1, 0

        # snps
        for i in [
            x
            for x in range(len(target))
            if target[x] != "-" and query[x] != "-" and target[x] != query[x]
        ]:

            # nucleotide level
            s = self.real_pos(i)
            yield "", target[s], query[s], s, s + 1, False

            # protein level
            if self.gff:
                for cds in self.gff.cds:
                    if cds.is_cds(s):
                        for j, codon_coords in enumerate(cds.iter_coords(True), 0):
                            if s in codon_coords:
                                ref = cds.aa[j]
                                alt = (
                                    query[codon_coords[0]]
                                    + query[codon_coords[1]]
                                    + query[codon_coords[2]]
                                )
                                alt = self.translate(alt, cds.translation_table)
                                if ref != alt:
                                    yield cds.symbol, cds.aa[j], alt, j, j + 1, 0

    @staticmethod
    def translate(seq, translation_table=1):
        """
        function to translate a nucleotide sequence.

        Notes
        -----
                If necessary, the given nucleotide sequence is shortened that its
                length is a multiple of 3.

        Example
        -------

        >>> algn = sonarALIGN(QRY_FASTA_FILE, REF_FASTA_FILE)
        >>> algn.translate("ATGTGAAA")
        'M*'

        Parameters
        ----------
        seq : str
                define the nucleotide sequence to translate
        translation_table : int
                define the genetic code table used for in silico translation (see
                https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) [ 1 ]

        Returns
        -------
        str
                translated amino acid sequence
        """
        l = len(seq)
        if l % 3 == 1:
            l = -1
        elif l % 3 == 2:
            l = -2
        return str(Seq.translate(seq[:l], table=translation_table))


class sonarDB(object):
    """
    this object provides sonarDB functionalities and intelligence

    Notes
    -----
            Please note, that genomic and protein coordinates are expected to be  and
            returned 0-based by this object, except for formatted profiles.
            While start or single coordinates are inclusive, end coordinates of
            ranges are exclusive, expressed in a mathematical notation: [start, end).
            Only in formatted profiles start and end coordinates are 1-based and both
            inclusive.

    Examples
    --------

    In this example the path to the database is stored in DOCTESTDB.

    >>> db = sonarDB(DOCTESTDB)

    Parameters
    ----------
    dbfile : str
            define a path to a non-existent or valid SONAR database file. If the
            file does not exist, a SONAR database is created.
    translation_table : int
            define the genetic code table used for in silico translation (see
            https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) [ 1 ]

    Attributes
    ----------
    db : str
            stores the absolute path to the used SONAR database file
    reffna : str
            stores the absolute path to the built-in FASTA file containing the reference
            genome sequence
    refgff : str
            stores the absolute path to the built-in GFF3 file containing the reference
            genome annotation
    translation_table : int
            stores the genetic code table used for in silico translation (see
            https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) [ 1 ]
    refseq : str
            stores the upper-case sequence of the built-in reference genome
    refdescr : str
            stores the FASTA header of the built-in reference genome
    refgffObj : object
            stores the sonarGFF object based on the built-in reference genome
            annotation
    iupac_nt_code : dict
            stores a dict with IUPAC one-letter nucleotide codes as keys and the
            respective set of matching explicit IUPAC one-letter nucleotide codes
            as values (e.g {"W": set('A', 'T')})
    iupac_explicit_nt_code : dict
            stores a set containing all non-ambiguous IUPAC one-letter nucleotide codes
    iupac_ambig_nt_code : set
            stores a set containing all ambiguous IUPAC one-letter nucleotide codes
    iupac_aa_code : dict
            stores a dict with IUPAC one-letter amino acid codes as keys and
            the respective set of matching IUPAC one-letter amino acids codes as values
    iupac_explicit_aa_code : dict
            stores a set containing all non-ambiguous IUPAC one-letter amino acid codes
    iupac_ambig_aa_code : dict
            stores a set containing all ambiguous IUPAC one-letter amino acid codes
    dna_var_regex : compiled re expression
            stores a compiled re expression that matches to nucleotide profiles but
            not to amino acid profiles
    aa_var_regex : compiled re expression
            stores a compiled re expression that matches to amino acid profiles but
            not to nucleotide profiles
    del_regex : compiled re expression
            stores a compiled re expression that matches to deletion profiles on
            nucleotide as well as on amino acid level.
    dnavar_grep_regex : compiled re expression
            stores a compiled re expression that matches to snp or dna insertion
            profiles with eference allele, genomic position and variant allele
            as groups.
    codedict : dict
            stores a dictionary with "dna" and "aa" containing the field name in the
            database that stores the profile data, the one letter code with and
            without ambiguities
    """

    def __init__(self, dbfile, translation_table=1):
        self.db = os.path.abspath(dbfile)
        self.__moduledir = os.path.dirname(os.path.realpath(__file__))
        self.reffna = os.path.join(self.__moduledir, "ref.fna")
        self.refgff = os.path.join(self.__moduledir, "ref.gff3")
        self.translation_table = translation_table
        self.__refseq = None
        self.__refdescr = None
        self.__refgffObj = None
        self.__iupac_nt_code = None
        self.__iupac_aa_code = None
        self.__iupac_explicit_nt_code = None
        self.__iupac_explicit_aa_code = None
        self.__iupac_ambig_nt_code = None
        self.__iupac_ambig_aa_code = None
        self.__terminal_letters_regex = re.compile("[A-Z]$")
        self.__dna_var_regex = None
        self.__aa_var_regex = None
        self.__del_regex = None
        self.__dnavar_grep_regex = None
        self.__codedict = None

    # PROPERTIES ON DEMAND

    @property
    def refseq(self):
        if not self.__refseq:
            record = SeqIO.read(self.reffna, "fasta")
            self.__refseq = self.harmonize(record.seq)
        return self.__refseq

    @property
    def refdescr(self):
        if not self.__refdescr:
            with open(self.reffna, "r") as handle:
                self.__refdescr = handle.readline().strip()[1:]
        return self.__refdescr

    @property
    def refgffObj(self):
        if not self.__refgffObj:
            self.__refgffObj = sonarGFF(
                self.refgff, self.reffna, self.translation_table
            )
        return self.__refgffObj

    @property
    def dna_var_regex(self):
        if self.__dna_var_regex is None:
            allowed_letters = "[" + "".join(self.iupac_nt_code.keys()) + "]"
            self.__dna_var_regex = re.compile(
                "^(?:(?:del:[0-9]+:[0-9]+)|(?:"
                + allowed_letters
                + "[0-9]+"
                + allowed_letters
                + "+))$"
            )
        return self.__dna_var_regex

    @property
    def dnavar_grep_regex(self):
        if self.__dnavar_grep_regex is None:
            self.__dnavar_grep_regex = re.compile("^([^0-9:]*)([0-9]+)([^0-9]*)$")
        return self.__dnavar_grep_regex

    @property
    def aa_var_regex(self):
        if self.__aa_var_regex is None:
            allowed_symbols = "(?:(?:" + ")|(?:".join(self.refgffObj.symbols) + "))"
            allowed_letters = (
                "[" + "".join(self.iupac_aa_code.keys()).replace("-", "") + "*~-" + "]"
            )
            self.__aa_var_regex = re.compile(
                "^"
                + allowed_symbols
                + ":(?:(?:del:[0-9]+:[0-9]+)|(?:"
                + allowed_letters
                + "[0-9]+"
                + allowed_letters
                + "+))$"
            )
        return self.__aa_var_regex

    @property
    def del_regex(self):
        if self.__del_regex is None:
            allowed_symbols = "(?:(?:" + ")|(?:".join(self.refgffObj.symbols) + "))"
            self.__del_regex = re.compile(
                "^(?:" + allowed_symbols + ":)?del:[0-9]+:[0-9]+$"
            )
        return self.__del_regex

    @property
    def iupac_nt_code(self):
        if self.__iupac_nt_code is None:
            self.__iupac_nt_code = {
                "A": set("A"),
                "C": set("C"),
                "G": set("G"),
                "T": set("T"),
                "R": set("AGR"),
                "Y": set("CTY"),
                "S": set("GCS"),
                "W": set("ATW"),
                "K": set("GTK"),
                "M": set("ACM"),
                "B": set("CGTB"),
                "D": set("AGTD"),
                "H": set("ACTH"),
                "V": set("ACGV"),
            }
            self.__iupac_nt_code["N"] = set(self.__iupac_nt_code.keys()) | set("N")
        return self.__iupac_nt_code

    @property
    def iupac_explicit_nt_code(self):
        if self.__iupac_explicit_nt_code is None:
            self.__iupac_explicit_nt_code = set(
                [x for x in self.iupac_nt_code if len(self.iupac_nt_code[x]) == 1]
            )
        return self.__iupac_explicit_nt_code

    @property
    def iupac_ambig_nt_code(self):
        if self.__iupac_ambig_nt_code is None:
            self.__iupac_ambig_nt_code = set(
                [x for x in self.iupac_nt_code if len(self.iupac_nt_code[x]) > 1]
            )
        return self.__iupac_ambig_nt_code

    @property
    def iupac_aa_code(self):
        if self.__iupac_aa_code is None:
            self.__iupac_aa_code = {
                "A": set("A"),
                "R": set("R"),
                "N": set("N"),
                "D": set("D"),
                "C": set("C"),
                "Q": set("Q"),
                "E": set("E"),
                "G": set("G"),
                "H": set("H"),
                "I": set("I"),
                "L": set("L"),
                "K": set("K"),
                "M": set("M"),
                "F": set("F"),
                "P": set("P"),
                "S": set("S"),
                "T": set("T"),
                "W": set("W"),
                "Y": set("Y"),
                "V": set("V"),
                "U": set("U"),
                "O": set("O"),
            }
            self.__iupac_aa_code.update(
                {
                    "B": set("DNB"),
                    "Z": set("EQZ"),
                    "J": set("ILJ"),
                    "Φ": set("VILFWYMΦ"),
                    "Ω": set("FWYHΩ"),
                    "Ψ": set("VILMΨ"),
                    "π": set("PGASπ"),
                    "ζ": set("STHNQEDKRζ"),
                    "+": set("KRH+"),
                    "-": set("DE-"),
                }
            )
            self.__iupac_aa_code["X"] = set(self.__iupac_aa_code.keys()) | set("X")
        return self.__iupac_aa_code

    @property
    def iupac_explicit_aa_code(self):
        if self.__iupac_explicit_aa_code is None:
            self.__iupac_explicit_aa_code = set(
                [x for x in self.iupac_aa_code if len(self.iupac_aa_code[x]) == 1]
            )
        return self.__iupac_explicit_aa_code

    @property
    def iupac_ambig_aa_code(self):
        if self.__iupac_ambig_aa_code is None:
            self.__iupac_ambig_aa_code = set(
                [x for x in self.iupac_aa_code if len(self.iupac_aa_code[x]) > 1]
            )
        return self.__iupac_ambig_aa_code

    @property
    def codedict(self):
        if self.__codedict is None:
            self.__codedict = {
                "dna": {
                    "field": "dna_profile",
                    "code": self.iupac_nt_code,
                    "explicit_code": self.iupac_explicit_nt_code,
                },
                "aa": {
                    "field": "aa_profile",
                    "code": self.iupac_aa_code,
                    "explicit_code": self.iupac_explicit_aa_code,
                },
            }

        return self.__codedict

    # DATA IMPORT

    @staticmethod
    def hash(seq):
        """
        static function to hash any sequence using SEGUID (SHA-1 hash of the
        upper-case sequence)

        Parameters
        ----------
        seq : str
                define a sequence to hash

        Returns
        -------
        str
                seguid

        """
        return seguid(seq)

    @staticmethod
    def harmonize(seq):
        """
        static function to return a sequence in upper case format and with T instead of U

        Parameters
        ----------
        seq : str
                define a sequence to harmonize

        Returns
        -------
        str
                sequence

        """
        return str(seq).strip().upper().replace("U", "T")

    def check_iupac_nt_code(self, seq):
        """
        returns a set of non-IUPAC characters present in a given sequence

        Parameters
        ----------
        seq : str
                define a sequence to check

        Returns
        -------
        str
                sequence

        """
        return set(seq).difference(self.iupac_nt_code.keys())

    def multi_process_fasta_wrapper(self, args):
        """
        wrapper function for sonarDB.process_fasta that accepts the needed
        parameters as list (which allows to be called by multiprocessing for
        parallelization) to add a genome sequences from a FASTA file. The FASTA
        file has to contain exactly one record.

        Parameters
        ----------
        args: list
                ordered list of the following arguments
        args[0] : str
                corresponds to fname in sonarDB.process_fasta
                define a valid FASTA file containing exactly one genome record to be
                added to the SONAR database
        args[1] : str
                corresponds to algnfile in sonarDB.process_fasta
                define a filename to permanently store the sequence alignment. Please
                consider, that an existing file will be overwritten. If None, a
                temporary file will be created and deleted after processing.
        args[2] : str
                corresponds to cache in sonarDB.process_fasta
                define a cache file (pickle format) that is used to permanently store
                processed data. Please consider, that an existing file will be
                overwritten. IfNone, a temporary file will be created and deleted after
                processing.
        args[3] : int
                timeout in seconds
                define a timeout in seconds for processing genomes
                integers below 1 deactivate the timeout.

        Returns
        -------
        tuple
                returns a tuple consisting of status and the hash of the processed
                genome sequence. Status False means TimeoutError (genome was not added
                to the database) while True means genome was successfully added.

        """
        fname, algnfile, picklefile, seqhash, timeout = args
        try:
            with sonarTimeout(seconds=timeout):
                self.process_fasta(fname, algnfile, picklefile)
        except TimeoutError:
            return False, seqhash
        else:
            return True, seqhash

    def process_fasta(self, fname, algnfile=None, pickle_file=None):
        """
        function to process a genome sequence from a single FASTA file, if
        the respective sequence is not in the database. The FASTA
        file has to contain exactly one record.

        Example
        -------

        In this example the path to the database is stored in DOCTESTDB.
        QRY_FASTA_FILE stores the path of a FASTA file conatining a
        B.1.1.7 prototype genome sequence.

        >>> a = os.remove(DOCTESTDB) if os.path.exists(DOCTESTDB) else None
        >>> db = sonarDB(DOCTESTDB)
        >>> data = db.process_fasta(QRY_FASTA_FILE)
        >>> data['acc']
        'b117'
        >>> data['descr']
        'b117 Ideal severe acute respiratory syndrome coronavirus 2 lineage B.1.1.7, complete genome'
        >>> data['dna_profile']
        'C3267T C5388A T6954C del:11288:9 del:21765:6 del:21991:3 A23063T C23271A C23604A C23709T T24506G G24914C C27972T G28048T A28111G G28280C A28281T T28282A C28977T'
        >>> data['prot_profile']
        'ORF1a:T1001I ORF1a:A1708D ORF1a:I2230T ORF1a:del:3675:3 ORF1b:T1001I ORF1b:A1708D ORF1b:I2230T ORF1b:del:3675:3 S:del:68:3 S:del:143:2 S:N501Y S:A570D S:P681H S:T716I S:S982A S:D1118H ORF8:Q27* ORF8:R52I ORF8:Y73C N:D3L N:S235F'

        Parameters
        ----------
        fname : str
                define a valid FASTA file containing exactly one genome record to be
                added to the SONAR database
        algnfile : str [ None ]
                define a filename to permanently store the sequence alignment. Please
                consider, that an existing file will be overwritten. If None, a
                temporary file will be created and deleted after processing.
        pickle_file : str [ None ]
                define a filname to store the dictionary in pickle format instead of
                returning it.  Please consider, that an existing file will be
                overwritten. If None, a temporary file will be created and deleted
                after processing.

        Returns
        -------
        dict
                if pickle_file is None a dictionary is returned, else there is no return
                value. The dictionary has following keys and values and can be directly
                used as input for the import_genome function of this class (**kwargs):
                        - acc: accession of processed genome
                        - descr: FASTA header of processed genome
                        - dnadiff: a list of nucleotide level variations (see sonarALIGN.dnadiff)
                        - aadiff: a list of amino acid level variations (see sonarALIGN.aadiff)
                        - dna_profile: the formatted nucleotide level profile (see sonarDB.build_profile)
                        - prot_profile: the formatted amino acid level profile (see sonarDB.build_profile)
                        - fs_profile: the dna_profile with frameshift mutations only
                        - seq: genome sequence
        """
        record = SeqIO.read(fname, "fasta")
        seq = self.harmonize(record.seq)
        seqhash = self.hash(seq)
        data = {"acc": record.id, "descr": record.description, "seqhash": seqhash}

        alignment = sonarALIGN(fname, self.reffna, algnfile, self.refgffObj)
        data["vars"] = [x for x in alignment.iter_vars()]

        if pickle_file:
            with open(pickle_file, "wb") as handle:
                pickle.dump(data, handle)
        else:
            data["seq"] = seq
            return data

    def import_genome_from_fasta_files(
        self, *fnames, dbm=None, msg=None, disable_progressbar=False
    ):
        """
        function to import genome sequence(s) from given FASTA file(s) to the
        SONAR database. Each FASTA file has to contain exactly one record.

        Example
        -------

        In this example the path to the database is stored in DOCTESTDB.
        QRY_FASTA_FILE stores the path of a FASTA file conatining a
        B.1.1.7 protoype genome sequence.

        >>> a = os.remove(DOCTESTDB) if os.path.exists(DOCTESTDB) else None
        >>> db = sonarDB(DOCTESTDB)
        >>> db.import_genome_from_fasta_files(QRY_FASTA_FILE, disable_progressbar=True)

        Parameters
        ----------
        *fnames : str
                define one or more valid FASTA files. Each file must contain
                exactly one genome record
        dbm : sonarDBManager object [ None ]
                define a sonarDBManager object to use for database transaction
        msg : str
                define a message used for the progress bar. If None, no progress
                bar is shown. [ None ]
        disable_progressbar : bool [ False ]
                define if the progress bar is shown (False) or not (True)
        """
        with ExitStack() as stack:
            if dbm is None:
                dbm = stack.enter_context(sonarDBManager(self.db))
            for i in tqdm(
                range(len(fnames)), desc=msg, disable_progressbar=disable_progressbar
            ):
                self.import_genome(**self.process_fasta(fnames[i]), dbm=dbm)

    def import_genome_from_cache(
        self, cachedir, acc_dict, dbm=None, msg=None, disable_progressbar=False
    ):
        """
        function to import data from a sonarCACHE directory to the SONAR database.

        Parameters
        ----------
        cachedir : str
                define a valid sonarCACHE directory
        acc_dict : dict
                define a dictionary (key: sequence hash, value: set of assigned accessions)
                to import to the database
        dbm : sonarDBManager object [ None ]
                define a sonarDBManager object to use for database transaction
        msg : str [ None ]
                define a message used for the progress bar. If None, no progress
                bar is shown
        disable_progressbar : bool [ False ]
                define if the progress bar is shown (False) or not (True)
        """
        seqhashes = list(acc_dict.keys())
        with ExitStack() as stack, sonarCache(cachedir) as cache:
            if dbm is None:
                dbm = stack.enter_context(sonarDBManager(self.db))
            for i in tqdm(range(len(seqhashes)), desc=msg, disable=disable_progressbar):
                seqhash = seqhashes[i]
                seq = cache.get_cached_seq(seqhash)
                preprocessed_data = cache.load_info(seqhash)
                for entry in acc_dict[seqhash]:
                    self.import_genome(
                        entry[0], entry[1], seqhash, preprocessed_data["vars"], seq, dbm
                    )

    def import_genome(self, acc, descr, seqhash, vars, seq=None, dbm=None):
        """
        function to import processed data to the SONAR database.

        Parameters
        ----------

        acc : str
                define the accession of the processed genome
        descr : str
                define the FASTA header of the processed genome
        seqhash : str
                define the hash (seguid) of the processed genome
        dnadiff : list
                define a sub list of nucleotide level variations (see sonarALIGN.dnadiff)
        aadiff : list
                define a sub list of amino acid level variations (see sonarALIGN.aadiff)
        dna_profile : str
                define the formatted nucleotide level profile (see sonarDB.build_profile)
        prot_profile : str
                define the formatted amino acid level profile (see sonarDB.build_profile)
        seq : str
                define the sequence of the processed genome (can be None, but then no paranoid test is done)
        dbm : sonarDBManager object [ None ]
                define a sonarDBManager object to use for database transaction
        """
        with ExitStack() as stack:
            if dbm is None:
                dbm = stack.enter_context(sonarDBManager(self.db))

            dbm.insert_genome(acc, descr, seqhash)
            dbm.insert_sequence(seqhash)

            for var in vars:
                dbm.insert_var(seqhash, *var)

            # if seq:
            # 	self.be_paranoid(acc, seq, auto_delete=True, dbm=dbm)

    # NOMENCLATURE

    def isdnavar(self, var):
        """
        function to validate nucleotide level profiles

        Examples
        --------

        >>> a = os.remove(DOCTESTDB) if os.path.exists(DOCTESTDB) else None
        >>> db = sonarDB(DOCTESTDB)
        >>> db.isdnavar("S:N501Y")
        False
        >>> db.isdnavar("A101T")
        True

        Parameters
        ----------

        var : str
                define the profile to validate

        Returns
        -------

        bool
                True if var is a valid nucleotide level profile otherwise False
        """
        return bool(self.dna_var_regex.match(var))

    def isaavar(self, var):
        """
        function to validate amino acid level profiles

        Examples
        --------

        >>> a = os.remove(DOCTESTDB) if os.path.exists(DOCTESTDB) else None
        >>> db = sonarDB(DOCTESTDB)
        >>> db.isaavar("S:N501Y")
        True
        >>> db.isaavar("A101T")
        False

        Parameters
        ----------

        var : str
                define the profile to validate

        Returns
        -------

        bool
                True if var is a valid amino acid level profile otherwise False
        """
        return bool(self.aa_var_regex.match(var))

    def isdel(self, var):
        """
        function to validate deletion profiles on both nucleotide and amino acid level

        Examples
        --------

        >>> a = os.remove(DOCTESTDB) if os.path.exists(DOCTESTDB) else None
        >>> db = sonarDB(DOCTESTDB)
        >>> db.isdel("del:100-118")
        False
        >>> db.isdel("del:100:18")
        True
        >>> db.isdel("ORF1b:del:5:2")
        True

        Parameters
        ----------

        var : str
                define the profile to validate

        Returns
        -------

        bool
                True if var is a deletion profile otherwise False
        """
        return bool(self.del_regex.match(var))

    # PROFILE BUILDING

    def build_profile(self, *vars):
        """
        function to build a valid variant profiles based on given variations

        Parameters
        ----------

        vars : list
                define for each variation to be considered by the profile a list
                with the following elements:
                 - reference nucleotide(s) or amino acid(s)
                 - alternative nucleotide(s) or amino acid(s)
                 - start position (0-based) related to the genome (nucleotide level profile) or
                   protein (amino acid level profile)
                 - end position (0-based) related to the genome (nucleotide level profile) or
                   protein (amino acid level profile) or None if single nucleotide/amino acid
                   polymorphism
                 - protein symbol (None in case of nucleotide level profiles)
                 - gene locus (None in case of nucleotide level profiles)

        Returns
        -------

        str
                valid variant profile
        """
        if len(vars) == 0:
            return ""
        profile = []
        if len(vars) == 1:
            this_ref, this_alt, this_start, this_end, this_protein, this_locus = vars[0]
            if this_alt == "" and this_end is None:
                this_end = this_start + len(this_ref)
        else:
            vars = sorted(vars, key=lambda x: (x[5], x[4], x[2]))
            for l in range(len(vars) - 1):
                (
                    this_ref,
                    this_alt,
                    this_start,
                    this_end,
                    this_protein,
                    this_locus,
                ) = vars[l]
                (
                    next_ref,
                    next_alt,
                    next_start,
                    next_end,
                    next_protein,
                    next_locus,
                ) = vars[l + 1]
                if this_alt != "":
                    var = self.format_var(
                        this_ref, this_alt, this_start, this_end, this_protein
                    )
                    profile.append(var)
                elif (
                    this_alt == ""
                    and next_alt == ""
                    and this_start + len(this_ref) == next_start
                    and this_protein == next_protein
                    and this_locus == next_locus
                ):
                    vars[l + 1] = (
                        this_ref + next_ref,
                        "",
                        this_start,
                        next_start + 1,
                        this_protein,
                        this_locus,
                    )
                else:
                    if this_alt == "" and this_end is None:
                        this_end = this_start + len(this_ref)
                    var = self.format_var(
                        this_ref,
                        this_alt,
                        this_start,
                        this_end,
                        this_protein,
                        this_locus,
                    )
                    profile.append(var)
            this_ref, this_alt, this_start, this_end, this_protein, this_locus = vars[
                l + 1
            ]
            if this_alt == "" and this_end is None:
                this_end = this_start + len(this_ref)
        var = self.format_var(
            this_ref, this_alt, this_start, this_end, this_protein, this_locus
        )
        if var not in profile:
            profile.append(var)

        return profile

    @staticmethod
    def format_var(ref, alt, start, end, protein=None, locus=None):
        """
        function to build a valid variant profile based on a single variation

        Parameters
        ----------

        ref : str
                define the reference nucleotide(s) or amino acid(s)
        alt : str
                define the alternative nucleotide(s) or amino acid(s)
        start : int
                define the start position (0-based) related to the genome (nucleotide
                level profile) or protein (amino acid level profile)
        end : int
                define the end position (0-based) related to the genome (nucleotide
                level profile) or protein (amino acid level profile) or None if
                single nucleotide/amino acid polymorphism
        protein : str
                define the protein symbol (None in case of nucleotide level profiles)
                [ None ]
        locus : str
                define the gene locus (None in case of nucleotide level profiles)
                [ None ]

        Returns
        -------

        str
                valid variant profile
        """
        if alt != "":
            coord = str(start + 1)
        else:
            ref = "del:"
            coord = str(start + 1) + ":" + str(end - start)
        protein = protein + ":" if protein else ""
        return protein + ref + coord + alt

    # FRAMESHIFT DETECTION

    def is_frameshift(self, dna_var):
        """
        function to check if a dna variant causes a frameshift in any annotated
        CDS.

        Returns
        -------
        bool
                True if dna variant causes a frameshift, otherwise False.
        """

        if dna_var.startswith("del:"):
            _, x, l = dna_var.split(":")
            x = int(x) - 1
            y = x + int(l)
            for cds in self.refgffObj.cds:
                if cds.is_frameshift_del(x, y):
                    return True
        else:
            match = self.dnavar_grep_regex.search(dna_var)
            x = int(match.group(2)) - 1
            l = len(match.group(3)) - 1
            if l % 3 != 0:
                for cds in self.refgffObj.cds:
                    if cds.is_frameshift_in(x, l):
                        return True
        return False

    def filter_frameshifts(self, dna_profile):
        """
        function to filter all frameshift mutations from a given dna_profile.

        Returns
        -------
        str
                dna_profile containing only frameshift mutations
        """
        if self.refgffObj and dna_profile.strip():
            return " ".join(
                [
                    x
                    for x in filter(None, dna_profile.split(" "))
                    if self.is_frameshift(x)
                ]
            )
        return ""

    # MATCHING

    def filter_ambig(self, profile, explicit_code, keep=None):
        """
        function to filter variations with ambiguities in the alternative allele
        from a valid nucleotide or amino acid level profile

        Parameters
        ----------

        profile : str
                valid nucleotide or amino acid level profile
        explicit_code : dict
                explicit IUPAC code dictionary to use (as provided by
                sonarDB.iupac_explicit_nt_code or sonarDB.iupac_explicit_aa_code)
        keep : list [ None ]
                list of single variation profiles to exclude from filtering

        Returns
        -------

        str
                valid variant profile
        """
        out = []
        keep = set(keep) if keep else set()
        for mutation in list(filter(None, profile.split(" "))):
            if mutation in keep or self.del_regex.search(mutation):
                out.append(mutation)
                continue
            match = self.__terminal_letters_regex.search(mutation)
            if (
                match
                and len(match.group(0)) == 1
                and match.group(0) not in explicit_code
            ):
                continue
            out.append(mutation)
        return " ".join(out)

    def pinpoint_mutation(self, mutation, code):
        """
        function to generate a set of all profiles consisting of
        non-ambiguous one-letter codes only that match to a given profile.
        If the given profile does not contain any ambiguities a list only
        containing the given profile is returned.

        Examples
        --------

        >>> a = os.remove(DOCTESTDB) if os.path.exists(DOCTESTDB) else None
        >>> db = sonarDB(DOCTESTDB)
        >>> sorted(db.pinpoint_mutation('A5001N', db.iupac_nt_code))
        ['A5001A', 'A5001B', 'A5001C', 'A5001D', 'A5001G', 'A5001H', 'A5001K', 'A5001M', 'A5001N', 'A5001R', 'A5001S', 'A5001T', 'A5001V', 'A5001W', 'A5001Y']
        >>> db.pinpoint_mutation('N501Y', db.iupac_aa_code)
        {'N501Y'}

        Parameters
        ----------

        mutation : str
                define a valid nucleotide or amino acid level profile that may contain
                ambiguities
        code : dict
                define the IUPAC code dictionary to use (as provided by
                sonarDB.iupac_nt_code or sonarDB.iupac_aa_code)

        Returns
        -------

        set
                set of profiles without ambiguities but matching to given profile
        """
        # extract ALT call from mutation profile
        match = self.__terminal_letters_regex.search(mutation)
        if not match:
            return {
                mutation,
            }
        match = match.group(0)

        # resolve ambiguities
        options = []
        for m in match:
            options.append(code[m])

        # generate the set of explicit mutations
        orig_stat = mutation[: -len(match)]
        return set(
            [mutation] + [orig_stat + "".join(x) for x in itertools.product(*options)]
        )

    def make_profile_explicit(self, profile):
        """
        function to replace ambiguous variants from a profile by the respective
        explicit variant descriptions and to sort profiles based on their level

        Parameters
        ----------

        profile : str
                define a valid nucleotide, amino acid or mixed  level profile that
                may contain ambiguities

        Returns
        -------

        dict
                dictionary with 'dna' or 'aa' as key and the respective list of
                explicit dna/protein level mutations
        """

        profile = set(profile)
        extended_profile = {"aa": [], "dna": []}
        for var in profile:
            key = "dna" if self.isdnavar(var) else "aa"
            extended_profile[key].extend(
                [v for v in self.pinpoint_mutation(var, self.codedict[key]["code"])]
            )
        return extended_profile

    def match(  # noqa: C901
        self,
        include_profiles=[],
        exclude_profiles=[],
        accessions=[],
        lineages=[],
        zips=[],
        dates=[],
        labs=[],
        sources=[],
        collections=[],
        technologies=[],
        platforms=[],
        chemistries=[],
        materials=[],
        software=None,
        software_version=None,
        min_ct=None,
        max_ct=None,
        ambig=False,
        count=False,
        frameshifts=0,
        debug=False,
        dbm=None,
    ):
        """
        function to search genomes in the SONAR database dependent on
        defined sequence and metadata profiles

        Parameters
        ----------

        include_profiles : list [ [] ]
                define a list of valid nucleotide, amino acid or mixed level profiles
                that may contain ambiguities to find genomes sharing respective
                profiles. Variations in each profile (sublist) are linked by AND operator
                while profiles from different sublists are linked by OR.
         exclude_profiles : list [ [] ]
                define a list of valid nucleotide, amino acid or mixed level profiles
                that may contain ambiguities to find genomes NOT sharing respective
                profiles. Variations in each profile (sublist) are linked by AND operator
                while profiles from different sublists are linked by OR.
        accessions : list [ [] ]
                list of accessions. Only genomes linked to accessions in this list
                will be matched. Accessions are negated when starting with ^. [ None ]
        lineages : list [ [] ]
                list of pangolin lineages. Only genomes assigend to a
                pangolin lineage in this list will be matched. Lineages are
                negated when starting with ^.
        zips : list [ [] ]
                list of zip codes. Only genomes linked to one of the given zip
                codes or whose linked zip code starts like one of the given
                zip codes are matched. zip codes are negated when starting with ^.
        dates : list [ [] ]
                define list of dates (YYYY-MM-DD) or date ranges (YYYY-MM-DD:YYYY-MM-DD).
                Only genomes linked to one of the given dates or date ranges are
                matched.
        sources : list [ [] ]
                list of data sources. Only genomes linked to a
                data source in this list will be matched. Data sources are
                negated when starting with ^.
        collections : list [ [] ]
                list of data collections. Only genomes linked to a
                data collection in this list will be matched. Data collections are
                negated when starting with ^.
        technologies : list [ [] ]
                list of sequencing technologies. Only genomes linked to a
                technology in this list will be matched. Technologies are
                negated when starting with ^.
        platforms : list [ [] ]
                list of sequencing platforms. Only genomes linked to a
                platform in this list will be matched. Platforms are
                negated when starting with ^.
        chemistries : list [ [] ]
                list of sequencing chemistries. Only genomes linked to a
                chemistry in this list will be matched. Chemistries are
                negated when starting with ^.
        software : str [ None ]
                software used for sequence reconstruction. Only genomes linked to the
                given software will be matched. Software is negated when starting with ^.
        software_version : str [ None ]
                software version used for sequence reconstruction. Only genomes linked
                to the given software version will be matched. Software version is
                negated when starting with ^. Needs software defined.
        materials : list [ [] ]
                list of sampling materials. Only genomes linked to a
                material in this list will be matched. Materials are
                negated when starting with ^.
        labs : list [ [] ]
                list of lab identifiers. Only genomes linked to a
                lab in this list will be matched. Labs are
                negated when starting with ^.
        min_ct : float [ None ]
                minimal ct value of genomes to match.
        max_ct : float [ None ]
                maximal ct value of genomes to match.
        ambig : bool [ False ]
                define if variant alleles including ambiguities should be shown (True)
                or not (False)
        count : bool [ False ]
                define if matched genomes should be counted (True) instead of collected
                (False).
        frameshifts : int [ 0 ]
                define if matched genomes have to conatin frameshift mutations (1)
                or have not to conatin frameshift mutations (-1) or frameshift mutations
                do not matter (0).
        debug : bool [ False ]
                activate debug mode for  sonarDBManager
        dbm : sonarDBManager object [ None ]
                define a sonarDBManager object to use for database transaction

        Returns
        -------

        list or int
                list of rows if count is False else number of rows as int. Each row
                represents a matching genome and is provided as dictionary with field
                names as keys.
        """

        # sanity check:
        check = []
        if include_profiles:
            check += [item for sublist in include_profiles for item in sublist]
        if exclude_profiles:
            check += [item for sublist in exclude_profiles for item in sublist]
        nonvalid = [x for x in check if not self.isdnavar(x) and not self.isaavar(x)]
        if nonvalid:
            sys.exit(
                "input error: Non-valid variant expression(s) entered: "
                + ", ".join(nonvalid)
            )

        if software_version and software is None:
            sys.exit(
                "input error: matching a given software version needs a software defined."
            )

        # adding conditions of profiles to include to where clause
        if include_profiles:
            include_profiles = [self.make_profile_explicit(x) for x in include_profiles]

        # adding conditions of profiles to exclude to where clause
        if exclude_profiles:
            exclude_profiles = [self.make_profile_explicit(x) for x in exclude_profiles]

        # adding accession, lineage, zips, and dates based conditions
        include_acc = [x for x in accessions if not x.startswith("^")]
        exclude_acc = [x[1:] for x in accessions if x.startswith("^")]

        include_lin = [x for x in lineages if not x.startswith("^")]
        exclude_lin = [x[1:] for x in lineages if x.startswith("^")]

        include_zip = [x for x in zips if not str(x).startswith("^")]
        exclude_zip = [x[1:] for x in zips if str(x).startswith("^")]

        include_dates = [x for x in dates if not str(x).startswith("^")]
        exclude_dates = [x[1:] for x in dates if str(x).startswith("^")]

        include_labs = [x for x in labs if not str(x).startswith("^")]
        exclude_labs = [x[1:] for x in labs if str(x).startswith("^")]

        include_source = [x for x in sources if not str(x).startswith("^")]
        exclude_source = [x[1:] for x in sources if str(x).startswith("^")]

        include_collection = [x for x in collections if not str(x).startswith("^")]
        exclude_collection = [x[1:] for x in collections if str(x).startswith("^")]

        include_technology = [x for x in technologies if not str(x).startswith("^")]
        exclude_technology = [x[1:] for x in technologies if str(x).startswith("^")]

        include_platform = [x for x in platforms if not str(x).startswith("^")]
        exclude_platform = [x[1:] for x in platforms if str(x).startswith("^")]

        include_chemistry = [x for x in chemistries if not str(x).startswith("^")]
        exclude_chemistry = [x[1:] for x in chemistries if str(x).startswith("^")]

        include_material = [x for x in materials if not str(x).startswith("^")]
        exclude_material = [x[1:] for x in materials if str(x).startswith("^")]

        if software:
            if not software.startswith("^"):
                include_software = software
                exclude_software = None
            else:
                include_software = None
                exclude_software = software[1:]
        else:
            include_software = None
            exclude_software = None

        if software_version:
            if not software_version.startswith("^"):
                include_software_version = software_version
                exclude_software_version = None
            else:
                include_software_version = None
                exclude_software_version = software_version[1:]
        else:
            include_software_version = None
            exclude_software_version = None

        # query
        with ExitStack() as stack:
            if dbm is None:
                dbm = stack.enter_context(sonarDBManager(self.db, readonly=True))
            dbm.debug = debug
            rows = dbm.match(
                include_profiles,
                exclude_profiles,
                include_acc,
                exclude_acc,
                include_lin,
                exclude_lin,
                include_zip,
                exclude_zip,
                include_dates,
                exclude_dates,
                include_labs,
                exclude_labs,
                include_source,
                exclude_source,
                include_collection,
                exclude_collection,
                include_technology,
                exclude_technology,
                include_platform,
                exclude_platform,
                include_chemistry,
                exclude_chemistry,
                include_material,
                exclude_material,
                include_software,
                exclude_software,
                include_software_version,
                exclude_software_version,
                min_ct,
                max_ct,
                count,
                frameshifts,
            )

        # remove ambiguities from database profiles if wished
        if not ambig and not count:
            keep = (
                [item for sublist in include_profiles for item in sublist]
                if include_profiles
                else None
            )
            for i in range(len(rows)):
                rows[i]["dna_profile"] = self.filter_ambig(
                    rows[i]["dna_profile"], self.iupac_explicit_nt_code, keep
                )
                rows[i]["aa_profile"] = self.filter_ambig(
                    rows[i]["aa_profile"], self.iupac_explicit_aa_code, keep
                )
        elif count:
            return rows[0]["count"]

        return rows

    # VALIDATION

    def restore_genome_using_dnavars(self, acc, dbm=None):
        """
        function to restore a genome sequence from the SONAR database using dna variation table

        Parameters
        ----------

        acc : str
                define the accesion of the genome that should be restored
        dbm : sonarDBManager object [ None ]
                define a sonarDBManager object to use for database transaction

        Raises
        ------

        Each variant site stored in the database is checked, if the linked reference
        nucleotide is correct. If not, program is terminated and an error shown.

        Returns
        -------

        tuple
                tuple of the FASTA header and sequence of the respective genome.
                None is returned if the given accession does not exist in the
                database.
        """

        with ExitStack() as stack:
            if dbm is None:
                dbm = stack.enter_context(sonarDBManager(self.db, readonly=True))
            rows = dbm.get_dna_vars(acc)
            if rows:
                prefix = ""
                qryseq = list(self.refseq)
                for row in rows:
                    if row["start"] is None:
                        continue
                    s = row["start"]
                    if s >= 0:
                        if row["ref"] != self.refseq[s]:
                            sys.exit(
                                "data error: data inconsistency found for '"
                                + acc
                                + "' ("
                                + row["ref"]
                                + " expected at position "
                                + str(s + 1)
                                + " of the reference sequence, got "
                                + self.refseq[s]
                                + ")."
                            )
                        qryseq[s] = row["alt"]
                    else:
                        prefix = row["alt"]
                return ">" + rows[0]["description"], prefix + "".join(qryseq)
            else:
                rows = dbm.get_genome(acc)
                if rows is None:
                    sys.exit("error: " + acc + " not found.")
                return ">" + rows["description"], self.refseq

    def restore_genome_using_dnaprofile(self, acc, dbm=None):
        """
        function to restore a genome sequence from the SONAR database using dna level profiles

        Parameters
        ----------

        acc : str
                define the accesion of the genome that should be restored
        dbm : sonarDBManager object [ None ]
                define a sonarDBManager object to use for database transaction

        Raises
        ------

        Each variant site stored in the database is checked, if the linked reference
        nucleotide is correct. If not, program is terminated and an error shown.

        Returns
        -------

        tuple
                tuple of the FASTA header and sequence of the respective genome.
                None is returned if the given accession does not exist in the
                database.
        """
        with ExitStack() as stack:
            if dbm is None:
                dbm = stack.enter_context(sonarDBManager(self.db, readonly=True))
            profile = dbm.get_dna_profile(acc)
            if profile:
                qryseq = list(self.refseq)
                prefix = ""
                for var in profile.strip().split(" "):
                    if var.startswith("del:"):
                        var = var.split(":")
                        s = int(var[1]) - 1
                        e = s + int(var[2])
                        for i in range(s, e):
                            qryseq[i] = ""
                    elif var:
                        match = self.dnavar_grep_regex.search(var)
                        pos = int(match.group(2)) - 1
                        ref = match.group(1)
                        alt = match.group(3)
                        if pos >= 0 and ref != self.refseq[pos]:
                            sys.exit(
                                "data error: data inconsistency found for '"
                                + acc
                                + "' ("
                                + ref
                                + " expected at position "
                                + str(pos + 1)
                                + " of the reference sequence, got "
                                + self.refseq[pos]
                                + ")."
                            )
                        if pos == -1:
                            prefix = alt
                        else:
                            qryseq[pos] = alt
                return prefix + "".join(qryseq)
            else:
                row = dbm.get_genome(acc)
                if row is None:
                    sys.exit("error: " + acc + " not found.")
                return ">" + row["description"], self.refseq

    def restore_alignment(self, acc, dbm=None):
        """
        function to restore a genome alignment from the SONAR database

        Parameters
        ----------

        acc : str
                define the accesion of the genome whose alignment versus the reference
                should be restored
        dbm : sonarDBManager object [ None ]
                define a sonarDBManager object to use for database transaction

        Raises
        ------

        Each variant site stored in the database is checked, if the linked reference
        nucleotide is correct. If not, program is terminated and an error shown.

        Returns
        -------

        tuple
                tuple of the FASTA header and aligned sequence of the respective genome
                followed by the FASTA header and aligned sequence of the reference genome.
                None is returned if the given accession does not exist in the
                database.
        """
        with ExitStack() as stack:
            if dbm is None:
                dbm = stack.enter_context(sonarDBManager(self.db, readonly=True))
            rows = dbm.get_dna_vars(acc)
        if rows:
            refseq = list(self.refseq)
            qryseq = refseq[:]
            for row in rows:
                if row["start"] is not None:
                    s = row["start"]
                    if s >= 0:
                        if row["ref"] != self.refseq[s]:
                            sys.exit(
                                "data error: data inconsistency found for '"
                                + acc
                                + "' ("
                                + row["ref"]
                                + " expected at position "
                                + str(s + 1)
                                + " of the reference sequence, got "
                                + refseq[s]
                                + ")."
                            )
                        qryseq[s] = "-" if not row["alt"] else row["alt"]
                        if len(row["alt"]) > 1:
                            refseq[s] += "-" * (len(row["alt"]) - 1)
                    else:
                        qryseq = [row["alt"]] + qryseq
                        refseq = ["-" * (len(row["alt"]))] + refseq
            return (
                ">" + rows[0]["description"],
                "".join(qryseq),
                ">" + self.dbobj.refdescr,
                "".join(refseq),
            )
        return None

    def be_paranoid(self, acc, orig_seq, auto_delete=False, dbm=None):  # noqa: C901
        """
        function to compare a given sequence with the respective sequence restored
        from the SONAR database

        Parameters
        ----------

        acc : str
                define the accesion of the genome that should be validated
        orig_seq : str
                define the sequence expected
        dbm : sonarDBManager object [ None ]
                define a sonarDBManager object to use for database transaction
        auto_delete : bool [ False ]
                define if the respective genome should be automatically deleted
                from the SONAR database if the test fails

        Returns
        -------

        bool
                True is returned if expected and restored sequences are not different
                otherwise False
        """
        orig_seq = self.harmonize(orig_seq)

        with ExitStack() as stack:
            if dbm is None:
                dbm = stack.enter_context(sonarDBManager(self.db))

            # dna table check
            s = self.restore_genome_using_dnavars(acc, dbm)[1]
            if orig_seq != s:
                if auto_delete:
                    dbm.delete_genome(acc)
                fd, path = mkstemp(suffix=".fna", prefix="paranoid_", dir=".")
                with open(path, "w") as handle:
                    handle.write(
                        ">original "
                        + acc
                        + "\n"
                        + orig_seq
                        + "\n"
                        + ">restored "
                        + acc
                        + "\n"
                        + orig_seq
                    )
                sys.exit(
                    "Good that you are paranoid: "
                    + acc
                    + " original and those restored from dna table do not match (sequences stored in "
                    + path
                    + ")."
                )

            # dna profile check
            s = self.restore_genome_using_dnaprofile(acc, dbm)
            if orig_seq != s:
                if auto_delete:
                    dbm.delete_genome(acc)
                fd, path = mkstemp(suffix=".fna", prefix="paranoid_", dir="./")
                with open(path, "w") as handle:
                    handle.write(
                        ">original "
                        + acc
                        + "\n"
                        + orig_seq
                        + "\n"
                        + ">restored "
                        + acc
                        + "\n"
                        + orig_seq
                    )
                sys.exit(
                    "Good that you are paranoid: "
                    + acc
                    + " original and those restored from its dna profile do not match (sequences stored in "
                    + path
                    + ")."
                )

            # frameshift checks
            row = self.match(accessions=[acc], ambig=True, dbm=dbm)[0]
            fs = set()
            for dna_var in row["dna_profile"].split(" "):
                if dna_var.strip() == "":
                    continue
                if self.is_frameshift(dna_var):
                    fs.add(dna_var)

            db_fs = set(filter(None, row["fs_profile"].split(" ")))
            missing_fs = [x for x in fs if x not in db_fs]
            wrong_fs = [x for x in db_fs if x not in fs]
            if wrong_fs:
                if auto_delete:
                    dbm.delete_genome(acc)
                fd, path = mkstemp(suffix=".csv", prefix="paranoid_", dir="./")
                with open(path, "w") as handle:
                    writer = csv.DictWriter(
                        handle, row.keys(), lineterminator=os.linesep
                    )
                    writer.writeheader()
                    writer.writerows([row])
                sys.exit(
                    "Good that you are paranoid: "
                    + ", ".join(wrong_fs)
                    + " not expected in frameshift profile of "
                    + acc
                    + " (profiles stored in "
                    + path
                    + ")."
                )

            if missing_fs:
                if auto_delete:
                    dbm.delete_genome(acc)
                fd, path = mkstemp(suffix=".csv", prefix="paranoid_", dir="./")
                with open(path, "w") as handle:
                    writer = csv.DictWriter(
                        handle, row.keys(), lineterminator=os.linesep
                    )
                    writer.writeheader()
                    writer.writerows([row])
                sys.exit(
                    "Good that you are paranoid: "
                    + ", ".join(missing_fs)
                    + " missing in frameshift profile of "
                    + acc
                    + " (profiles stored in "
                    + path
                    + ")."
                )

        return True

    @staticmethod
    def get_version():
        return __version__


class sonarCache:
    """
    this object manages permanent and temporary file caches

    Notes
    -----

    This class should be included via context manager to ensure that accession
    index is written and cleaning temporary objects is performed after abnormal
    program termination.

    In the SONAR cache for each unique sequence that has been cached a FASTA file
    containing the sequence. That files are named by the slugified hash of the
    sequence they contain while the used FASTA header represent the hash. Pre-processed
    data provided by the sonarDB.process_fasta is stored in info files als named by
    the slugified hash of the respective sequence they are related to (PICKLE format).
    The link between sequence hash and accession(s) is stored in the cache attribute and,
    when closing the cache, written to the index file (PICKLE format).

    Parameters
    ----------
    dir : str
            define a path to an non-existent, empty or valid SONAR cache directory.
            If None, a temporary cache directoryis created and deleted after use.
            [ None ]

    Attributes
    ----------
    dirname : str
            stores the absolute path to the cache directory
    temp : bool
            stores True if the cache is temporary and will be deleted after use
            otherwise False
    cache : dict
            stores a dictionary whose keys are hashes of cached genome sequences and
            and values tuples of linked accessions and FASTA headers

    """

    def __init__(self, dir=None):
        self.temp = not bool(dir)
        self.cache = defaultdict(set)
        self._fasta_ext = ".fasta"
        self._info_ext = ".info"
        self._algn_ext = ".algn"

        if self.temp:
            self.dirname = mkdtemp(prefix=".sonarCache_")
        else:
            self.dirname = os.path.abspath(dir)
            self.checkdir()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        if [exc_type, exc_value, exc_traceback].count(None) != 3:
            print("warning:", file=sys.stderr)
            print(traceback.format_exc(), file=sys.stderr)
        if os.path.isdir(self.dirname) and self.temp:
            shutil.rmtree(self.dirname)

    def checkdir(self):
        if not os.path.isdir(self.dirname):
            os.makedirs(self.dirname)

    @staticmethod
    def slugify(string):
        """
        function to provide a file-system- and collision-safe representation of
        a given string

        Parameters
        ----------

        string : str
                define the string to slugify

        Returns
        -------

        str
                a file-system- and collision-safe representation of
                the original string
        """
        return base64.urlsafe_b64encode(string.encode("UTF-8")).decode("UTF-8")

    @staticmethod
    def deslugify(string):
        return base64.urlsafe_b64decode(string).decode("utf-8")

    @staticmethod
    def get_seqhash_from_fasta_name(fname):
        return sonarCache.deslugify(os.path.basename(fname))

    def iter_fasta(self, fname):
        """
        function to iterate records of a given FASTA file

        Parameters
        ----------

        fname : str
                define the path to a valid FASTA file

        Returns
        -------

        tuple
                for each record a tuple is returned consisting of
                 - accession
                 - FASTA header
                 - upper-case sequence
        """
        for record in SeqIO.read(fname, "fasta"):
            yield record.id, record.description, str(record.seq).upper()

    def read_cached_fasta(self, seqhash):
        record = SeqIO.read(self.get_fasta_fname(seqhash), "fasta")
        return record.id, record.description[1:], str(record.seq).upper()

    def get_cached_filename(self, seqhash, ext=""):
        basename = self.slugify(seqhash)
        return os.path.join(self.dirname, basename[:2], basename + ext)

    def get_fasta_fname(self, seqhash):
        return self.get_cached_filename(seqhash, self._fasta_ext)

    def get_algn_fname(self, seqhash):
        return self.get_cached_filename(seqhash, self._algn_ext)

    def get_info_fname(self, seqhash):
        return self.get_cached_filename(seqhash, self._info_ext)

    def prep_cached_files(self, seqhash):
        fasta = self.get_fasta_fname(seqhash)
        algn = self.get_algn_fname(seqhash)
        info = self.get_info_fname(seqhash)
        os.makedirs(os.path.dirname(fasta), exist_ok=True)
        return fasta, algn, info

    def load_info(self, seqhash):
        with open(self.get_info_fname(seqhash), "rb") as handle:
            return pickle.load(handle, encoding="bytes")

    def write_info(self, seqhash, data={}):
        data["seqhash"] = seqhash
        with open(self.get_info_fname(seqhash), "wb") as handle:
            pickle.dump(data, handle)

    def add_seq(self, seqhash, seq):
        """
        function to add a sequence to the cache

        Parameters
        ----------

        seqhash : str
                define the seqhash of the sequence
        seq : str
                define the sequence

        """
        fasta, align, info = self.prep_cached_files(seqhash)

        # check for sequence hash collision
        if not os.path.isfile(fasta):
            with open(fasta, "w") as handle:
                handle.write(">" + seqhash + os.linesep + seq)
        elif seq != self.read_cached_fasta(seqhash)[2]:
            sys.exit("cache error: sequence hash collision for hash '" + seqhash + "'.")

    def get_cached_seqhashes(self):
        return set(self.cache.keys())

    def iter_cached_fasta_files(self):
        for x in self.cache:
            yield self.get_fasta_fname(x)

    def get_cached_seq(self, seqhash):
        return self.read_cached_fasta(seqhash)[-1]


if __name__ == "__main__":
    import doctest

    global DOCTESTDIR, DOCTESTDB, QRY_FASTA_FILE, REF_FASTA_FILE
    print("sonarDB", sonarDB.get_version())
    print("performing unit tests ...")
    with TemporaryDirectory() as tmpdirname:
        this_path = os.path.dirname(os.path.realpath(__file__))
        DOCTESTDIR = tmpdirname
        DOCTESTDB = os.path.join(DOCTESTDIR, "testdb")
        QRY_FASTA_FILE = os.path.join(this_path, "doctest_b117.fna")
        QRY_PICKLE_FILE = os.path.join(this_path, "doctest_b117.pickle")
        REF_FASTA_FILE = os.path.join(this_path, "ref.fna")
        REF_GFF_FILE = os.path.join(this_path, "ref.gff3")
        print(doctest.testmod(verbose=False))
