#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

import os
import re
import sys
import argparse
import sqlite3
from sqlite3 import Error
from Bio.SeqUtils.CheckSum import seguid
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Emboss.Applications import StretcherCommandline
from packaging import version
import shutil
import base64
from collections import OrderedDict, defaultdict
import pickle
from tqdm import tqdm
from urllib.parse import quote as urlquote
from math import floor, ceil
from tempfile import mkstemp, TemporaryDirectory, mkdtemp
import traceback
import itertools
import signal
import csv
from time import sleep

class sonarTimeout():
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
	error_message: str
		define error message to be shown [ 'Timeout' ]

    Attributes
    ----------

    seconds : int
        time in seconds when a TimeoutError is raised
	error_message : str
		error message to be shown [ 'Timeout' ]
	"""
	def __init__(self, seconds, error_message='Timeout'):
		self.seconds = seconds
		self.error_message = error_message

	def __enter__(self):
		signal.signal(signal.SIGALRM, self.handle_timeout)
		signal.alarm(self.seconds)

	def __exit__(self, type, value, traceback):
		signal.alarm(0)

	def handle_timeout(self, signum, frame):
		raise TimeoutError(self.error_message)


class sonarFiler():
	"""
	this class is a helper class providing a creating (temporary) file handler
	for writing in a given context

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
    fname : str
        define designated file name to open. If None, a temporary file is
		created and deleted after use. [ None ]

    Attributes
    ----------
    name : str
        stores file name
	basename : str
		stores file basename
	path : str
		stores absolute file path
	tmp : bool
		stores True if it is a temporary file else False
	"""
	def __init__(self, fname = None):
		self.fname = fname
		self.tmp = True if fname is None else False

	def __enter__(self):
		if 	self.fname is None:
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
	this object stores coding sequence information on one coding sequence(CDS)

	Notes
	-----
    	Please note, that genomic coordinates are processed and returned 0-based
		by this object. While start or single coordinates are inclusive,
		end coordinates of ranges are exclusive, expressed in a mathematical
		notation: [start, end)

	Examples
	--------

	Initiating an sonarCDS object:

	>>> cds = sonarCDS("ORF1", 155, 170, "+", "ATGTTATGAATGGCC")

	Accessing amino acid sequence (genetic code 1):

	>>> cds.aa
	'ML'

	Accessing CDS coordinates or genome range:

	>>> cds.coords
	(155, 170)
	>>> cds.range
	range(155, 170)

    Parameters
    ----------
    symbol : str
        define the symbol of the protein encoded
		(e.g. ORF1ab)
    start : int
        define the genomic start coordinate (0-based, inclusive)
    end : int
        define the genomic end coordinate (0-based, exclusive)
	strand : {'+', '-'}
		define the genomic strand encoded on
	seq : str
		define the nucleotide sequence
	locus: str
		define the gene locus accession number [ None ]
	translation_table : int
		define the genetic code table used for in silico translation (see
		https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) [ 1 ]

    Attributes
    ----------
    symbol : str
        stores the symbol of the protein encoded
	start : int
		stores the genomic start coordinate (0-based, inclusive). The start
		coordinate is always lower than the end coordinate.
	end : int
		stores the genomic end coordinate (0-based, exclusive). The end
		coordinate is always greater than the start coordinate.
	coords : tuple
		stores a tuple of start and end coordinate
	range :	range
		stores an range from start to end coordinate
	nuc : str
		stores the nucleotide sequence
	aa : str
		stores the in silico translated amino acid sequence from start to the
		first in-frame stop codon.
	translation_table : int
		stores the genetic code table used for in silico translation (see
		https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) [ 1 ]
	"""

	def __init__(self, symbol, start, end, strand, seq, locus=None, translation_table=1):
		self.symbol = symbol
		self.locus = locus
		self.start = min(start, end) # inclusive
		self.end = max(start, end) # exclusive
		self.strand = strand
		self.nuc = seq
		self.translation_table = translation_table
		self.__aa = None

	@property
	def aa(self):
		nuc = self.nuc if self.strand == "+" else str(Seq(self.nuc).reverse_complement())
		if self.__aa is None:
			l = len(nuc)
			if l%3 == 1:
				l = -1
			elif l%3 == 2:
				l = -2
			self.__aa = str(Seq.translate(nuc[:l], table=self.translation_table, to_stop=True))
		return self.__aa

	@property
	def coords(self):
		return self.start, self.end

	@property
	def range(self):
		return range(self.start, self.end)

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
    fna : int
        define a path to a valid FASTA file storing the nucleotide
		sequence of the annotated genome
	translation_table : int
		define the genetic code table used for in silico translation of CDS (see
		https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) [ 1 ]

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
	symbols : list
		stores a list of protein symbols

	"""

	def __init__(self, gff3, fna, translation_table=1):
		self.translation_table = translation_table
		self.cds = self.process_gff3(gff3, fna)
		self.coords = { x.symbol: (x.start, x.end) for x in self.cds }
		self.symbols = [x.symbol for x in self.cds ]

	def iscds(self, x, y=None):
		"""
		function to check if a given genomic coordinate (range) overlaps with an
		annotated coding sequence (CDS).

		Examples
		--------

		>>> gff=sonarGFF(REF_GFF_FILE, REF_FASTA_FILE)
		>>> gff.iscds(21562)
		True
		>>> gff.iscds(25384)
		False
		>>> gff.iscds(25380, 25384)
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
		r = range(x, y)
		for start, end in self.coords.values():
			if x < end and start < y:
				return True
		return False

	def convert_pos_to_dna(self, protein_symbol, x):
		"""
		function to convert a protein-based coordinate to a genomic coordinate

		Examples
		--------

		>>> os.chdir(os.path.dirname(os.path.realpath(__file__)))
		>>> gff=sonarGFF(REF_GFF_FILE, REF_FASTA_FILE)
		>>> gff.convert_pos_to_dna("S", 4)
		21574

		Parameters
		----------
		protein_symbol : str
			define the protein based on its symbol
		x : int
			define the position within the protein (0-based, inclusive)

		Returns
		-------
		int
			lower-bound genomic coordinate (0-based) where the respective protein
			position is encoded
		"""
		return x*3 + self.coords[protein_symbol][0]

	def convert_pos_to_aa(self, x, y=None):
		"""
		function to convert a genomic coordinate (range) to amino acid positions
		of the rexpective protein encoded

		Examples
		--------

		>>> os.chdir(os.path.dirname(os.path.realpath(__file__)))
		>>> gff=sonarGFF(REF_GFF_FILE, REF_FASTA_FILE)
		>>> gff.convert_pos_to_aa(265,274)
		('ORF1ab', 0, 3)

		Parameters
		----------
		x : int
			genomic (start) coordinate (0-based, inclusive)
		y : int
			genomic end coordniate (0-based, exclusive) [ None ]

		Returns
		-------
		tuple
		   tuple of respective protein symbol, 0-based start position within the
		   protein (inclusive), 0-based end position within the protein (exclusive)
		   or None if givben position is not in an annotated CDS.

		"""
		y = x + 1  if y is None else y
		for symbol, coords in self.coords.items():
			if x >= coords[0] and y <= coords[1]:
				x = floor((x - coords[0])/3)
				if y is None:
					return (symbol, x, None)
				y = floor((y - coords[0])/3)
				return (symbol, x, y)
		return None

	def process_gff3(self, gff, fna):
		"""
		function to parse CDS from a given GFF3 file

		Examples
		--------

		>>> os.chdir(os.path.dirname(os.path.realpath(__file__)))
		>>> gff = sonarGFF(REF_GFF_FILE, REF_FASTA_FILE)
		>>> gff.coords == {'ORF1ab': (265, 21555), 'S': (21562, 25384), 'ORF3a': (25392, 26220), 'E': (26244, 26472), 'M': (26522, 27191), 'ORF6': (27201, 27387), 'ORF7a': (27393, 27759), 'ORF7b': (27755, 27887), 'ORF8': (27893, 28259),'N': (28273, 29533), 'ORF10': (29557, 29674)}
		True

		Parameters
		----------
		gff : str
			path to a valid GFF3 file storing the genome annotation
		fna : str
			path to a valid FASTA file storing the respective genome sequence

		Returns
		-------
		dict
		  dictionary with protein symbols as keys and respective 0-based
		  genomic coordinate tuples (start always lower than end coordinate,
		  start coordinate inclusive, end coordinate exclusive)

		"""

		symbol_regex = re.compile("gene=([^;]+)(?:;|$)")
		locus_regex = re.compile("locus_tag=([^;]+)(?:;|$)")

		record = SeqIO.read(fna, "fasta")
		gseq = str(record.seq).upper()

		with open(gff, "r") as handle:
			cds = set()
			for line in handle:
				fields = line.rstrip("\r\n").split("\t")
				if line.startswith("#") or len(fields) < 7:
					continue
				if fields[2] == "CDS":
					symbol = symbol_regex.search(fields[-1])
					symbol = symbol.groups(1)[0] if symbol else None
					locus = locus_regex.search(fields[-1])
					locus = locus.groups(1)[0] if locus else None
					strand = fields[6]
					s = int(fields[3])-1
					e = int(fields[4])
					dna = gseq[s:e]
					if strand == "-":
						dna = str(Seq.reverse_complement(dna))
					cds.add((symbol, s, e, strand, dna, locus, self.translation_table))
		cds = sorted(cds, key=lambda x: x[1])
		return [sonarCDS(*x) for x in cds]

class sonarALIGN(object):
	"""
	this object performs pairwise sequence alignments and stores respective
	informations.

	Notes
	-----
    	Please note, that genomic coordinates are processed and returned 0-based
		by this object. While start or single coordinates are inclusive,
		end coordinates of ranges are exclusive, expressed in a mathematical
		notation: [start, end)

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
        define a path to a valid FASTA file storing the target or reference
		genome sequence
	out_file : str
		define a path to an output file that will store the FASTA formatted
		alignment. Please consider, that an existing file will be overwritten!
		If None, a temporary file is used and deleted after processing. [ None ]
	sonarGFFObj : object
		define a sonarGFF object storing the annotated CDS of the reference genome
		[ None ]

    Attributes
    ----------
	aligned_query : str
		stores the aligned upper-case query sequence
    aligned_target : str
        stores the aligned upper-case target or reference sequence
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
		Accordingly to the VCF format, InDels are expressed considering the upstream
		base as anchor. As a special case, an insertion at the start of the sequence
		has no anchor and a genomic coordinate of -1. The last two tuple elements
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
		Accordingly to the VCF format, InDels are expressed considering the upstream
		base as anchor. As a special case, an insertion at the start of the sequence
		has no anchor and a genomic coordinate of -1. The last two tuple elements
		are always None to keep the length according to tuples stored in aadiff.
	"""

	def __init__(self, query_file, target_file, out_file = None, sonarGFFObj = None):
		self.aligned_query, self.aligned_target = self.align_dna(query_file, target_file, out_file)
		self.gff = sonarGFFObj if sonarGFFObj else None
		self._indel_regex = re.compile("[^-]-+")
		self._codon_regex = re.compile(".-*.-*.-*")
		self._starting_gap_regex = re.compile("^-+")
		self.__dnadiff = None
		self.__aadiff = None
		self.__target_coords_matrix = None

	@property
	def dnadiff(self):
		if self.__dnadiff is None:
			self.__dnadiff = [ x for x in self.iter_dna_vars() ]
		return self.__dnadiff

	@property
	def aadiff(self):
		if self.__aadiff is None:
			self.__aadiff = [ x for x in self.iter_aa_vars() ]
		return self.__aadiff

	@property
	def _target_coords_matrix(self):
		if self.__target_coords_matrix is None:
			self.__target_coords_matrix = [len(x.group()) for x in re.finditer(".-*", self.aligned_target)]
		return self.__target_coords_matrix

	def use_stretcher(self, query_file, target_file, out_file, gapopen= 16, gapextend = 4, right_align = True):
		"""
		function to perform a pairwise aligment using EMBOSS Stretcher

		Parameters
		----------
		query_file : str
			define a path to a valid FASTA file storing the query sequence
		target_file : str
			define a path to a valid FASTA file storing the target or reference
			sequence
		out_file : str
			define a path to a file that will store the alignment. Please consider,
			that an existing file will be overwritten.
		gapopen : int
			define penalty for gap opening [ 16 ]
		gapextend : int
			define penalty for gap extension [ 4 ]

		Returns
		-------
		list
		  list of aligned query and target sequence, in that order
		"""
		temp = True if not out_file else False
		if temp:
			handle, out_file = mkstemp()
		cline = StretcherCommandline(asequence=query_file, bsequence=target_file, gapopen=gapopen, gapextend=gapextend, outfile=out_file, aformat="fasta")
		stdout, stderr = cline()
		alignment = [str(x.seq) for x in SeqIO.parse(out_file, "fasta")]
		if temp:
			os.remove(out_file)
		if right_align:
			alignment = self.right_align_gaps(*alignment)
		return alignment

	def right_align_gaps(self, query, target):
		"""
		function to align gaps to the right in two aligned sequences

		Parameters
		----------
		seq : str
			define an aligned sequence to align gaps

		Returns
		-------
		list
		  aligned query and target sequence with right-aligned gaps, in that order
		"""
		l = len(query)-1
		for match in re.finditer("-+", query):
			s = match.start()-1
			e = match.end()-1
			g = "-" * (e-s)
			while s >= 0 and e < l and query[s] == target[e]:
				query = query[:s] + g + query[s] + query[e+1:]
				s -= 1
				e -= 1
		for match in re.finditer("-+", target):
			s = match.start()-1
			e = match.end()-1
			g = "-" * (e-s)
			while s >= 0 and e < l and target[s] == query[e]:
				target = target[:s] + g + target[e:]
				s -= 1
				e -= 1
		return query, target

	def align_dna(self, query_file, target_file, out_file=None, gapopen = 16, gapextend = 4, right_align = True):
		"""
		function to perform the default pairwise nucleotide aligment

		Parameters
		----------
		query_file : str
			define a path to a valid FASTA file storing the query sequence
		target_file : str
			define a path to a valid FASTA file storing the target or reference
			sequence
		out_file : str
			define a path to a file that will store the alignment. Please consider,
			that an existing file will be overwritten.
		gapopen : int
			define penalty for gap opening [ 16 ]
		gapextend : int
			define penalty for gap extension [ 4 ]

		Returns
		-------
		list
		  list of aligned query and target sequence
		"""
		return self.use_stretcher(query_file, target_file, out_file, gapopen, gapextend, right_align)

	def real_pos(self, x):
		"""
		function to convert an alignment position to the real position in the
		target or reference sequence.

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
		  corresponding position (0-based) of the target/reference sequence
		"""
		return x - self.aligned_target[:x+1].count("-")

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

	def iter_dna_vars(self):
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
			has no anchor and a genomic coordinate of -1. The last two tuple elements
			are always None to keep the length according to tuples stored in aadiff.
		"""
		target = self.aligned_target
		query = self.aligned_query

		# query overhead in front
		match = self._starting_gap_regex.match(target)
		if match:
			yield "", query[:match.end()], -1, None, None, None

		# insertions
		isites = set()
		for match in self._indel_regex.finditer(target):
			isites.add(match.start())
			s = self.real_pos(match.start())
			yield match.group()[0], query[match.start():match.end()], s, None, None, None

		# deletions and snps
		for i, pair in enumerate(zip(target, query)):
			if pair[0] != "-" and pair[0] != pair[1] and i not in isites:
				s = self.real_pos(i)
				l = len(pair[1])
				e = None if l == 1 else s + l
				yield pair[0], pair[1].replace("-", ""), s, e, None, None

	def iter_aa_vars(self):
		"""
		function to iterate variations on amino acid level.

		Example
		-------

		In this example the QRY_FASTA_FILE, REF_FASTA_FILE, and REF_GFF_FILE
		variables store	the path of FASTA files containing the query and
		reference genome sequences as well as the reference genome annotation,
		in that order. The reference is NC_045512.2 while the query is a B.1.1.7
		prototype sequence without S:D613G variation.

		Please consider, that a sonarGFF is needed to consider annotation and
		deduce amino acid level profiles.

		>>> gff = sonarGFF(REF_GFF_FILE, REF_FASTA_FILE)
		>>> algn = sonarALIGN(QRY_FASTA_FILE, REF_FASTA_FILE, sonarGFFObj=gff)
		>>> for x in algn.iter_aa_vars():
		... 	print(x)
		('T', 'I', 1000, None, 'ORF1ab', 'GU280_gp01')
		('A', 'D', 1707, None, 'ORF1ab', 'GU280_gp01')
		('I', 'T', 2229, None, 'ORF1ab', 'GU280_gp01')
		('S', '', 3674, None, 'ORF1ab', 'GU280_gp01')
		('G', '', 3675, None, 'ORF1ab', 'GU280_gp01')
		('F', '', 3676, None, 'ORF1ab', 'GU280_gp01')
		('I', '', 67, None, 'S', 'GU280_gp02')
		('H', '', 68, None, 'S', 'GU280_gp02')
		('V', '', 69, None, 'S', 'GU280_gp02')
		('V', '', 142, None, 'S', 'GU280_gp02')
		('Y', '', 143, None, 'S', 'GU280_gp02')
		('N', 'Y', 500, None, 'S', 'GU280_gp02')
		('A', 'D', 569, None, 'S', 'GU280_gp02')
		('P', 'H', 680, None, 'S', 'GU280_gp02')
		('T', 'I', 715, None, 'S', 'GU280_gp02')
		('S', 'A', 981, None, 'S', 'GU280_gp02')
		('D', 'H', 1117, None, 'S', 'GU280_gp02')
		('Q', '*', 26, None, 'ORF8', 'GU280_gp09')
		('R', 'I', 51, None, 'ORF8', 'GU280_gp09')
		('Y', 'C', 72, None, 'ORF8', 'GU280_gp09')
		('D', 'L', 2, None, 'N', 'GU280_gp10')
		('S', 'F', 234, None, 'N', 'GU280_gp10')

		Returns
		-------
		iterator of tuples
			each tuple represents a amino acid level variation and consists of:
				 - target nucleotide
				 - query nucleotide(s)
				 - target or reference start position (0-based
				 - target or reference end position (0-based)
				 - protein symbol
				 - gene locus
			Accordingly to the VCF format, InDels are expressed considering the upstream
			base as anchor. As a special case, an insertion at the start of the sequence
			has no anchor and a genomic coordinate of -1. The last two tuple elements
			are always None to keep the length according to tuples stored in aadiff.
		"""
		if self.gff:
			for cds in self.gff.cds:
				s = self.align_pos(cds.start)
				e = self.align_pos(cds.end)
				if cds.strand == "+":
					query = self.aligned_query[s:e]
					target = self.aligned_target[s:e]
				else:
					query = str(Seq.reverse_complement(self.aligned_query[s:e]))
					target = str(Seq.reverse_complement(self.aligned_target[s:e]))

				for match in self._codon_regex.finditer(target):
					s = match.start()
					e = match.end()
					tcodon = match.group().replace("-", "")
					qcodon = query[s:e]
					taa = self.translate(tcodon, cds.translation_table)
					if "-" in qcodon:
						yield taa, "", int(s/3), None, cds.symbol, cds.locus
					else:
						qaa = self.translate(qcodon, cds.translation_table)
						if qaa != taa:
							e = None if len(qaa) == 1 else int(s/3) + len(qaa)
							yield (taa, qaa, int(s/3), e, cds.symbol, cds.locus)

	def translate(self, seq, translation_table=1):
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
		if l%3 == 1:
			l = -1
		elif l%3 == 2:
			l = -2
		return str(Seq.translate(seq[:l], table=translation_table))

class sonarDBManager():
	"""
	This object provides a sonarDB SQLite manager handler managing connections and
	providing context-safe transaction control.

	Notes
	-----
    	This object should be called using a context manager to ensure rollbacks
		on abnormal termination.

	Example
	-------

	In this example the DOCTESTDB variable store the path to a database file

	>>> with sonarDBManager(DOCTESTDB) as dbm:
	... 	pass

    Parameters
    ----------

    dbfile : str
        define a path to a non-existent or valid SONAR database file. If the
		file does not exist, a SONAR database is created.
	timeout : int
		define busy timeout [ -1 ],
	readonly : bool
		define if the connection should be read-only [ True ]

    Attributes
    ----------
    dbfile : str
        stores the path to the used SONAR database file.
	connection : object
		stores the SQLite3 connection
	cursor : method
		stores the SQLite3 cursor
	"""

	def __init__(self, dbfile, timeout=-1, readonly=False):
		self.dbfile = os.path.abspath(dbfile)
		self.connection = None
		self.cursor = None
		self.__timeout = timeout
		self.__mode = "ro" if readonly else "rwc"
		self.__uri = "file:" + urlquote(self.dbfile)

	def __enter__(self):
		if not os.path.isfile(self.dbfile) or os.stat(self.dbfile).st_size == 0:
			self.create_scheme()
		self.connection, self.cursor = self.connect()
		self.start_transaction()
		return self

	def __exit__(self, exc_type, exc_value, exc_traceback):
		if [exc_type, exc_value, exc_traceback].count(None) != 3:
			print("warning:", file=sys.stderr)
			print(traceback.format_exc(), file=sys.stderr)
			if self.__mode == "rwc":
				print("rollback", file=sys.stderr)
				self.rollback()
		elif self.__mode == "rwc":
			self.connection.commit()
		self.close()

	def __del__(self):
		if self.connection:
			self.close()

	def connect(self):
		con = sqlite3.connect(self.__uri + "?mode=" + self.__mode, self.__timeout, isolation_level = None, uri = True)
		con.row_factory = self.dict_factory
		cur = con.cursor()
		return con, cur

	def start_transaction(self):
		self.cursor.execute("BEGIN DEFERRED")

	def commit(self):
		self.connection.commit()

	def rollback(self):
		self.connection.rollback()

	def close(self):
		self.connection.close()

	def create_scheme(self):
		with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), "db.sqlite"), 'r') as handle:
			sql = handle.read()
		with sqlite3.connect(self.__uri + "?mode=rwc", uri = True) as con:
			con.executescript(sql)

	def select(self, table, fieldList=['*'], whereClause=None, whereVals=None, orderby=None, fetchone=False, print_sql=False):
		"""
		generic function to perform a SELECT query on the SONAR database.

		Parameters
		----------
		table : str
			define table to select data from
		fieldList : list
			define a list of field names to select. By default all fields are selected
			[ ['*'] ]
		whereClause : str
			define a where expression for conditional select. The expression has
			not to be introduced by 'where'. The actual values will be inserted
			syntax-safely and, thus, have to be replaced by ? within the clause
			(e.g. 'accession = ?'). [ None ]
		whereVals : list
			list of ordered values to be inserted into the whereClause. [ None ]
		orderby : str
			define an 'order by' expression to order results [ None ]
		fetchone : bool
			define wether only one (set True) or all (set False) rows should be
			returned. [ False ]
		print_sql : bool
			define if assembled query should be printed on screen
			(e.g. for debugging) [ False ]

		Returns
		-------
		list
			list of selected rows. Each row is represented as a dictionary with
			selected field names as keys.
		"""
		sql = "SELECT " + "".join(fieldList) + " FROM " + table
		if whereClause:
			sql += " WHERE " + whereClause
		else:
			whereVals = []
		if orderby:
			sql += " ORDER BY " + orderby
		if print_sql:
			print(sql)
		self.cursor.execute(sql, whereVals)
		if fetchone:
			return self.cursor.fetchone()
		return self.cursor.fetchall()

	def delete(self, table, whereClause=None, whereVals=None, print_sql=False):
		"""
		generic function to perform a DELETE request on the SONAR database.

		Parameters
		----------
		table : str
			define table to delete data from
		whereClause : str
			define a where expression for conditional deletions. The expression has
			not to be introduced by 'where'. The actual values will be inserted
			syntax-safely and, thus, have to be replaced by ? within the clause
			(e.g. 'accession = ?'). [ None ]
		whereVals : list
			list of ordered values to be inserted into the whereClause. [ None ]
		print_sql : bool
			define if assembled query should be printed on screen
			(e.g. for debugging) [ False ]
		"""
		sql = "DELETE FROM " + table
		if whereClause:
			sql += " WHERE " + whereClause
		else:
			whereVals = []
		if print_sql:
			print(sql)
		self.cursor.execute(sql, whereVals)

	def insert(self, table, fieldList, valList, ignore=False, print_sql=False):
		"""
		generic function to perform an INSERT request on the SONAR database.

		Parameters
		----------
		table : str
			define table to insert data to
		fieldList : list
			define a list of field names into which the data should be inserted
		valList : list
			define a list of values to be inserted. The values will be inserted
			in the same order as they appear and the field names have been defined
			in fieldList.
		ignore : bool
			define if INSERT can be ignored [ False ]
		print_sql : bool
			define if assembled query should be printed on screen
			(e.g. for debugging) [ False ]

		Returns
		-------
		int
			id of modified row
		"""
		ignore = "" if not ignore else "OR IGNORE "
		sql = "INSERT " + ignore + "INTO " + table + "(" + ", ".join(fieldList) + ") VALUES (" + ", ".join(['?']*len(valList)) + ")"
		if print_sql:
			print(sql)
		self.cursor.execute(sql, valList)
		if self.cursor.lastrowid:
			return self.cursor.lastrowid
		where = " AND ".join([x + "= ?" for x in fieldList])
		rowid = self.select(table, ['rowid'], where, valList)[0]['rowid']
		return rowid

	def update(self, table, fieldList, valList, whereClause, whereVals, print_sql=False):
		"""
		generic function to perform an UPDATE query on the SONAR database.

		Parameters
		----------
		table : str
			define table to select data from
		fieldList : list
			define a list of field names to be updated
		valList : list
			define a list of values to be inserted. The values will be inserted
			in the same order as they appear and the field names have been defined
			in fieldList.
		whereClause : str
			define a where expression. The expression has
			not to be introduced by 'where'. The actual values will be inserted
			syntax-safely and, thus, have to be replaced by ? within the clause
			(e.g. 'accession = ?').
		whereVals : list
			list of ordered values to be inserted into the whereClause.
		orderby : str
			define an 'order by' expression to order results [ None ]
		fetchone : bool
			define wether only one (set True) or all (set False) rows should be
			returned. [ False ]
		print_sql : bool
			define if assembled query should be printed on screen
			(e.g. for debugging) [ False ]
		"""
		sql = "UPDATE " + table + " SET " + ", ".join([ str(x) + "= ?" for x in fieldList ]) + " WHERE " + whereClause
		if print_sql:
			print(sql)
		self.cursor.execute(sql, valList + whereVals)

	@staticmethod
	def optimize(dbfile):
		with sqlite3.connect(dbfile) as con:
			con.executescript("VACUUM")

	@staticmethod
	def dict_factory(cursor, row):
		d = OrderedDict()
		for idx, col in enumerate(cursor.description):
			d[col[0]] = row[idx]
		return d

class sonarDB(object):
	"""
	this object provides sonarDB functionalities and intelligence

	Notes
	-----
    	Please note, that genomic and protein coordinates are processed and
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
		genome sequences
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
		stores a dict of IUPAC one-letter nucleotide codes (as keys) whereas
		the values are the respective set of matching IUPAC one-letter nucleotide
		codes (e.g "W" -> set('A', 'T'))
	iupac_explicit_nt_code : dict
		stores a set of non-ambiguous IUPAC one-letter nucleotide codes
	iupac_ambig_nt_code : set
		stores a set of ambiguous IUPAC one-letter nucleotide codes
	iupac_aa_code : dict
		stores a dict of IUPAC one-letter amnino acid codes (as keys) whereas
		the values are the respective set of matching IUPAC one-letter amino
		acids codes
	iupac_explicit_aa_code : dict
		stores a set of non-ambiguous IUPAC one-letter amino acid codes
	iupac_ambig_aa_code : dict
		stores a set of ambiguous IUPAC one-letter amino acid codes
	dna_var_regex : compiled re expression
		stores a compiled re expression that matches to nucleotide profiles but
		not to amino acid profiles
	aa_var_regex : compiled re expression
		stores a compiled re expression that matches to amino acid profiles but
		not to nucleotide profiles
	del_regex : compiled re expression
		stores a compiled re expression that matches to deletion profiles on
		nucleotide as well as on amino acid level.
	"""
	def __init__(self, dbfile, translation_table = 1):
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

	# PROPERTIES ON DEMAND

	@property
	def refseq(self):
		if not self.__refseq:
			with open(self.reffna, "r") as handle:
				self.__refseq = "".join([x.strip().upper() for x in handle.readlines()[1:]])
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
			self.__refgffObj = sonarGFF(self.refgff, self.reffna, self.translation_table)
		return self.__refgffObj

	@property
	def dna_var_regex(self):
		if self.__dna_var_regex is None:
			allowed_letters = "[" + "".join(self.iupac_nt_code.keys()) + "]"
			self.__dna_var_regex = re.compile("^(?:(?:del:[0-9]+:[0-9]+)|(?:" + allowed_letters + "[0-9]+" + allowed_letters + "+))$")
		return self.__dna_var_regex

	@property
	def aa_var_regex(self):
		if self.__aa_var_regex is None:
			allowed_symbols = "(?:(?:" + ")|(?:".join(self.refgffObj.symbols) + "))"
			allowed_letters = "[" + "".join(self.iupac_aa_code.keys()).replace("-", "") + "-" + "]"
			self.__aa_var_regex = re.compile("^" + allowed_symbols + ":(?:(?:del:[0-9]+:[0-9]+)|(?:" + allowed_letters + "[0-9]+" + allowed_letters + "+))$")
		return self.__aa_var_regex

	@property
	def del_regex(self):
		if self.__del_regex is None:
			allowed_symbols = "(?:(?:" + ")|(?:".join(self.refgffObj.symbols) + "))"
			self.__del_regex = re.compile("^(?:" + allowed_symbols + ":)?del:[0-9]+:[0-9]+$")
		return self.__del_regex

	@property
	def iupac_nt_code(self):
		if self.__iupac_nt_code is None:
			self.__iupac_nt_code = { "A": set("A"), "C": set("C"), "G": set("G"), "T": set("T"), "R": set("AGR"), "Y": set("CTY"), "S": set("GCS"), "W": set("ATW"), "K": set("GTK"), "M": set("ACM"), "B": set("CGTB"), "D": set("AGTD"), "H": set("ACTH"), "V": set("ACGV") }
			self.__iupac_nt_code['N'] = set(self.__iupac_nt_code.keys()) | set("N")
		return self.__iupac_nt_code

	@property
	def iupac_explicit_nt_code(self):
		if self.__iupac_explicit_nt_code is None:
 			self.__iupac_explicit_nt_code = set([ x for x in self.iupac_nt_code if len(self.iupac_nt_code[x]) == 1 ])
		return self.__iupac_explicit_nt_code

	@property
	def iupac_ambig_nt_code(self):
		if self.__iupac_ambig_nt_code is None:
 			self.__iupac_ambig_nt_code = set([ x for x in self.iupac_nt_code if len(self.iupac_nt_code[x]) > 1 ])
		return self.__iupac_ambig_nt_code

	@property
	def iupac_aa_code(self):
		if self.__iupac_aa_code is None:
			self.__iupac_aa_code = { "A": set("A"), "R": set("R"), "N": set("N"), "D": set("D"), "C": set("C"), "Q": set("Q"), "E": set("E"), "G": set("G"), "H": set("H"), "I": set("I"), "L": set("L"), "K": set("K"), "M": set("M"), "F": set("F"), "P": set("P"), "S": set("S"), "T": set("T"), "W": set("W"), "Y": set("Y"), "V": set("V"), "U": set("U"), "O": set("O") }
			self.__iupac_aa_code.update({"B": set("DNB"), "Z": set("EQZ"), "J": set("ILJ"), "Φ": set("VILFWYMΦ"), "Ω": set("FWYHΩ"), "Ψ": set("VILMΨ"), "π": set("PGASπ"), "ζ": set("STHNQEDKRζ"), "+": set("KRH+"), "-": set("DE-") })
			self.__iupac_aa_code['X'] = set(self.__iupac_aa_code.keys()) | set("X")
		return self.__iupac_aa_code

	@property
	def iupac_explicit_aa_code(self):
		if self.__iupac_explicit_aa_code is None:
 			self.__iupac_explicit_aa_code = set([ x for x in self.iupac_aa_code if len(self.iupac_aa_code[x]) == 1 ])
		return self.__iupac_explicit_aa_code

	@property
	def iupac_ambig_aa_code(self):
		if self.__iupac_ambig_aa_code is None:
 			self.__iupac_ambig_aa_code = set([ x for x in self.iupac_aa_code if len(self.iupac_aa_code[x]) > 1 ])
		return self.__iupac_ambig_aa_code

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

	def multi_process_fasta_wrapper(self, args):
		"""
		wrapper function for sonarDB.process_fasta that accepts the needed
		parameters as list (which allows to be called by multiprocessing for
		parallelization) to add a genome sequences from a FASTA file. The FASTA
		file has to contain exactly one record.

		Parameters
		----------
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

		Returns
		-------
		tuple
			returns a tuple consisting of status and the hash of the processed
			genome sequence. Status False means TimeoutError (genome was not added
			to the database) while True means genome was successfully added.

		"""
		fname, algnfile, cache, seqhash, timeout = args
		try:
			with sonarTimeout(seconds=timeout):
				self.process_fasta(fname, algnfile, cache)
		except TimeoutError:
			return False, seqhash
		else:
			return True, seqhash

	def process_fasta(self, fname, algnfile=None, cache=None):
		"""
		function to add a genome sequence from a single FASTA file. The FASTA
		file has to contain exactly one record.

		Example
		-------

		In this example the path to the database is stored in DOCTESTDB.
		QRY_FASTA_FILE stores the path of a FASTA file conatining a
		B.1.1.7 protoype genome sequence.

		>>> a = os.remove(DOCTESTDB) if os.path.exists(DOCTESTDB) else None
		>>> db = sonarDB(DOCTESTDB)
		>>> data = db.process_fasta(QRY_FASTA_FILE)
		>>> data[0]
		'b117'
		>>> data[1]
		'b117 Ideal severe acute respiratory syndrome coronavirus 2 lineage B.1.1.7, complete genome'
		>>> data[5]
		'C3267T C5388A T6954C del:11288:8 del:21765:5 del:21991:2 A23063T C23271A C23604A C23709T T24506G G24914C C27972T G28048T A28111G G28280C A28281T T28282A'
		>>> data[6]
		'N:D3L ORF8:Q27* ORF8:R52I S:del:68:2 ORF8:Y73C S:del:143:1 N:S235F S:N501Y S:A570D S:P681H S:T716I S:S982A ORF1ab:T1001I S:D1118H ORF1ab:A1708D ORF1ab:I2230T ORF1ab:del:3675:1'

		Parameters
		----------
		fname : str
			define a valid FASTA file containing exactly one genome record to be
			added to the SONAR database
		algnfile : str
			define a filename to permanently store the sequence alignment. Please
			consider, that an existing file will be overwritten. If None, a
			temporary file will be created and deleted after processing.
		cache : str
			define a cache file (pickle format) that is used to permanently store
			processed data. Please consider, that an existing file will be
			overwritten. IfNone, a temporary file will be created and deleted after
			processing.

		Returns
		-------
		list
			a list is returned only if data is not stored i a cache file. If cache
			is None a list is returned consisting of:
				- accession of processed genome
				- FASTA header of processed genome
				- a sub list of nucleotide level variations (see sonarALIGN.dnadiff)
				- a sub list of amino acid level variations (see sonarALIGN.aadiff)
				- the formatted nucleotide level profile (see sonarDB.build_profile)
				- the formatted amino acid level profile (see sonarDB.build_profile)
				- the complete sequence of the processed genome

		"""
		with sonarDBManager(self.db, readonly=True) as dbm:
			record = SeqIO.read(fname, "fasta")
			acc = record.id
			descr = record.description
			seq = str(record.seq.upper())
			seqhash = self.hash(seq)
			if not self.seq_exists(seqhash, dbm):
				alignment = sonarALIGN(fname, self.reffna, algnfile, self.refgffObj)
				dnadiff = alignment.dnadiff
				aadiff = alignment.aadiff
				dna_profile = self.build_profile(*dnadiff)
				prot_profile = self.build_profile(*aadiff)
			else:
				alignment = None
				dnadiff = None
				aadiff = None
				dna_profile = None
				prot_profile = None

		if cache is None:
			return [acc, descr, seqhash, dnadiff, aadiff, dna_profile, prot_profile, seq]

		with open(cache, "wb") as handle:
			pickle.dump([acc, descr, seqhash, dnadiff, aadiff, dna_profile, prot_profile], handle)

	def import_genome_from_fasta(self, *fnames, msg=None):
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
		>>> db.import_genome_from_fasta(QRY_FASTA_FILE)

		Parameters
		----------
		*fnames : str
			define one or more valid FASTA files. Each file must contain
			exactly one genome record
		msg : str
			define a message used for the progress bar. If None, no progress
			bar is shown [ None ]
		"""
		if not msg is None:
			rng = tqdm(range(len(fnames)), desc = msg)
		else:
			rng = range(len(fnames))

		with sonarDBManager(self.db) as dbm:
			for i in rng:
				data = self.process_fasta(fnames[i])
				self.import_genome(*data, dbm=dbm)

	def import_genome_from_pickle(self, *fnames, msg=None):
		"""
		function to import data from pre-processed pickle files (as created by
		sonarCACHE) to the SONAR database. Each PICKLE  file has to contain
		exactly one record.

		Example
		-------
		In this example the path to the database is stored in DOCTESTDB.
		QRY_PICKLE_FILE stores the pre-processed data obtained from a
		B.1.1.7 protoype genome sequence.

		>>> a = os.remove(DOCTESTDB) if os.path.exists(DOCTESTDB) else None
		>>> db = sonarDB(DOCTESTDB)
		>>> db.import_genome_from_pickle(QRY_PICKLE_FILE)

		Parameters
		----------
		*fnames : str
			define one or more valid sonarCHACHE PICKLE files. Each file must
			contain exactly one genome record
		msg : str
			define a message used for the progress bar. If None, no progress
			bar is shown [ None ]
		"""
		if not msg is None:
			rng = tqdm(range(len(fnames)), desc = msg)
		else:
			rng = range(len(fnames))

		with sonarDBManager(self.db) as dbm:
			for i in rng:
				with open(fnames[i], 'rb') as handle:
					data = pickle.load(handle, encoding="bytes")
				self.import_genome(*data, dbm=dbm)

	def import_genome_from_cache(self, cachedir, msg=None):
		"""
		function to import data from a sonarCACHE directory to the SONAR database.

		Parameters
		----------
		cachedir : str
			define a valid sonarCACHE directory
		msg : str
			define a message used for the progress bar. If None, no progress
			bar is shown [ None ]
		"""
		with sonarCache(cachedir) as cache:
			seqhashes = list(cache.cache.keys())
			if not msg is None:
				rng = tqdm(range(len(seqhashes)), desc = msg)
			else:
				rng = range(len(seqhashes))

			with sonarDBManager(self.db) as dbm:
				for i in rng:
					seqhash = seqhashes[i]
					seq = cache.get_cached_seq(seqhash)
					preprocessed_data = cache.load_info(seqhash)[2:]
					for entry in cache.cache[seqhash]:
						self.import_genome(*entry, *preprocessed_data, seq, dbm=dbm)

	def import_genome(self, acc, descr, seqhash, dnadiff, aadiff, dna_profile, prot_profile, seq, dbm):
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
			define the sequence of the processed genome
		dbm : str
			define a sonarDBManager object to use for database transaction
		"""
		if not self.genome_exists(acc, descr, seqhash, dbm):
			self.insert_genome(acc, descr, seqhash, dbm)
		if not self.seq_exists(seqhash, dbm):
			self.insert_sequence(seqhash, dbm)
			self.insert_profile(seqhash, dna_profile, prot_profile, dbm)
			for ref, alt, s, e, _, __ in dnadiff:
				varid = self.get_dna_varid(ref, alt, s, e, dbm)
				if varid is None:
					varid = self.insert_dna_var(ref, alt, s, e, dbm)
					self.insert_sequence2dna(seqhash, varid, dbm)

			for ref, alt, s, e, protein, locus in aadiff:
				varid = self.get_prot_varid(protein, locus, ref, alt, s, e, dbm)
				if varid is None:
					varid = self.insert_prot_var(protein, locus, ref, alt, s, e, dbm)
				self.insert_sequence2prot(seqhash, varid, dbm)

			self.be_paranoid(acc, seq, dbm, auto_delete=False)

	def insert_genome(self, acc, descr, seqhash, dbm):
		return dbm.insert('genome', ['accession', 'description', 'seqhash'], [acc, descr, seqhash])

	def insert_sequence(self, seqhash, dbm):
		return dbm.insert('sequence', ['seqhash'], [seqhash], ignore=True)

	def insert_profile(self, seqhash, dna_profile, aa_profile, dbm):
		dna_profile = " " + dna_profile.strip() + " "
		aa_profile = " " + aa_profile.strip() + " "
		return dbm.insert('profile', ['seqhash', 'dna_profile', 'aa_profile'], [seqhash, dna_profile, aa_profile], ignore=True)

	def insert_dna_var(self, ref, alt, start, end, dbm):
		return dbm.insert('dna', ['varid', 'ref', 'alt', 'start', 'end'], [None, ref, alt, start, end], ignore=True)

	def insert_sequence2dna(self, seqhash, varid, dbm):
		return dbm.insert('sequence2dna', ('seqhash', 'varid'), [seqhash, varid], ignore=True)

	def insert_prot_var(self, protein, locus, ref, alt, start, end, dbm):
		return dbm.insert('prot', ['varid', 'protein', 'locus', 'ref', 'alt', 'start', 'end'], [None, protein, locus, ref, alt, start, end], ignore=True)

	def insert_sequence2prot(self, seqhash, varid, dbm):
		return dbm.insert('sequence2prot', ['seqhash', 'varid'], [seqhash, varid], ignore=True)

	# DELETE DATA

	def delete_accession(self, acc, dbm):
		dbm.delete("genome", whereClause="accession = ?", whereVals=[acc], print_sql=False)

	# DATA SELECTS

	def genome_exists(self, acc, descr, seqhash, dbm):
		row = self.get_genome_data(acc, dbm=dbm)
		if not row or (descr is not None and row[0]['description'] != descr) or (seqhash is not None and row[0]['seqhash'] != seqhash):
			return False
		return True

	def seq_exists(self, seqhash, dbm):
		return dbm.select('sequence', fieldList=['COUNT(*)'], whereClause="seqhash = ?", whereVals=[seqhash], orderby=None, fetchone=True, print_sql=False)['COUNT(*)'] > 0

	def get_genome_data(self, *accs, dbm):
		where = " OR ".join(["accession = ?"] * len(accs))
		return dbm.select('genome', whereClause=where, whereVals=accs)

	def get_dna_varid(self, ref, alt, start, end, dbm):
		row = dbm.select('dna', whereClause="ref = ? AND alt = ? AND start = ? and end = ?", whereVals=[ref, alt, start, end], fetchone=True)
		if row:
			return row['varid']
		return None

	def get_prot_varid(self, protein, locus, ref, alt, start, end, dbm):
		row = dbm.select('prot', whereClause="protein = ? AND locus = ? AND ref = ? AND alt = ? AND start = ? AND end = ?", whereVals=[protein, locus, ref, alt, start, end], fetchone=True)
		if row:
			return row['varid']
		return None

	def get_dna_vars(self, acc, dbm):
		return dbm.select('dna_view', whereClause="accession = ?", whereVals=[acc], orderby="start DESC")

	def iter_table(self, table):
		sql = dbm.select('dna', whereClause="ref = ? AND alt = ? AND start = ? and end = ?", whereVals=[ref, alt, start, end], fetchone=True)
		for row in dbm.select(table):
			yield row

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
		>>> db.isdel("ORF1ab:del:5:2")
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
		else:
			vars = sorted(vars, key=lambda x: x[2])
			for l in range(len(vars)-1):
				this_ref, this_alt, this_start, this_end, this_protein, this_locus = vars[l]
				next_ref, next_alt, next_start, next_end, next_protein, next_locus = vars[l+1]
				if this_alt != "":
					var = self.format_var(this_ref, this_alt, this_start, this_end, this_protein)
					if var not in profile:
						profile.append(var)
				elif this_alt == "" and this_start + len(this_ref) == next_start and this_protein == next_protein and this_locus == next_locus:
					vars[l+1] = (this_ref + next_ref, vars[l+1][1], this_start, next_start, this_protein, this_locus)
				else:
					var = self.format_var(this_ref, this_alt, this_start, this_end, this_protein)
					if var not in profile:
						profile.append(var)
		var = self.format_var(this_ref, this_alt, this_start, this_end, this_protein)
		if var not in profile:
			profile.append(var)

		return " ".join(profile)

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
		if end is None:
			coord = str(start+1)
		else:
			ref = "del:"
			coord = str(start+1) + ":" + str(end-start)
		protein = protein + ":" if protein else ""
		return protein + ref + coord + alt

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
		keep : list
			list of single variation profiles to exclude from filtering [ None ]

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
			if match and len(match.group(0)) == 1 and  match.group(0) not in explicit_code:
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
			return mutation
		match = match.group(0)

		# resolve ambiguities
		options = []
		for m in match:
			options.append(code[m])

		# generate the set of explicit mutations
		orig_stat = mutation[:-len(match)]
		return set([mutation] + [ orig_stat + "".join(x) for x in itertools.product(*options) ])

	def match_builder(self, profile, exclusive=False, negate=False):
		"""
		function to build a where clause matching to nucleotide, amino
		acid or mixed level profiles.

		Parameters
		----------

		profile : str
			define a valid nucleotide or amino acid level profile that may contain
			ambiguities
		exclusive : bool
			define if additional variations are allowed for the matched genomes (False)
			or not (True) [ False ]
		negate : bool
			define if the query should be negated (True) or not (False) [ False ]

		Returns
		-------

		str
			where clause allowing matching of a given variant profile
		"""

		# configuring
		condition = " " if not negate else " NOT "
		config = {
				   "dna": {
							"field": "dna_profile",
					  		"code": self.iupac_nt_code,
					  		"explicit_code": self.iupac_explicit_nt_code
						  },
					"aa": {
							"field": "aa_profile",
							"code": self.iupac_aa_code,
							"explicit_code": self.iupac_explicit_aa_code
					}
				  }

		# generating clause
		clause = []
		vars = set(profile)
		dna_profile_length = 1
		aa_profile_length = 1
		for var in vars:
			if self.isdnavar(var):
				dna_profile_length += 1 + len(var)
				conf = config['dna']
			else:
				aa_profile_length += 1 + len(var)
				conf = config['aa']
			for v in self.pinpoint_mutation(var, conf['code']):
				clause.append(conf['field'] + condition + "LIKE '% " + v + " %'")

		clause = " AND ".join(clause)

		# add exclusiveness condition not allowing additional mutations
		if exclusive:
			if dna_profile_length > 1 and aa_profile_length > 1:
				sys.exit("input error: exclusive profiles must be defined on dna OR protein level only.")
			if dna_profile_length > 0:
				clause += " AND length(" + conf['field'] + ") = " + str(dna_profile_length)
			elif aa_profile_length > 0:
				clause += " AND length(" + conf['field'] + ") = " + str(aa_profile_length)
		return clause

	def match(self, include_profiles=None, exclude_profiles=None, accessions=None, lineages=None, zips=None, dates=None, exclusive=False, ambig=False, show_sql=False):
		"""
		function to match genomes in the SONAR database

		Parameters
		----------

		include_profiles : list
			define a list of valid nucleotide, amino acid or mixed level profiles
			that may contain ambiguities to find genomes sharing respective
			profiles. Variations in each profile are linked by AND operator
			while the different profiles are linked by OR. [ None ]
		 exclude_profiles : list
			define a list of valid nucleotide, amino acid or mixed level profiles
			that may contain ambiguities to find genomes NOT sharing respective
			profiles. Variations in each profile are linked by AND operator
			while the different profiles are linked by OR.
		accessions : list
			list of accessions. Only genomes of accessions in this list will be
			matched. Accessions are negated when starting with ^. [ None ]
		lineages : list
			list of pangolin lineages. Only genomes assigend to the respective
			pangolin lineage  in this list will be matched. Lineages are
			negated when starting with ^. [ None ]
		zips : list
			list of zip codes. Only genomes linked to one of the given zip
			codes or whose linked zip code starts like one of the given
			zip codes are matched. zip codes are negated when starting with ^.
			[ None ]
		dates : list
			define list of dates (YYYY-MM-DD) or date ranges (YYYY-MM-DD:YYYY-MM-DD).
			Only genomes linked to one of the given dates or date ranges are
			matched.
		exclusive : bool
			define if additional variations are allowed for the matched genomes (False)
			or not (True) [ False ]
		ambig : bool
			define if variations including ambiguities should be filtered (True)
			from teh profiles shown or not (False) [ False ]
		show_sql : bool
			define if the assembled SQLite query should be shown on screen (for
			debugging) [ False ]

		Returns
		-------

		list
			list of rows. Each row represents a matching genome and is provided as
			dictionary with field names as keys.
		"""

		clause = []
		vals =[]

		#sanity check:
		check = []
		if include_profiles:
			check += [item for sublist in include_profiles for item in sublist]
		if exclude_profiles:
			check += [item for sublist in exclude_profiles for item in sublist]
		nonvalid = [ x for x in check if not self.isdnavar(x) and not self.isaavar(x) ]
		if nonvalid:
			sys.exit("input error: Non-valid variant expression(s) entered: " + ", ".join(nonvalid))

		# adding conditions of profiles to include to where clause
		if include_profiles:
			includes = []
			for profile in include_profiles:
				includes.append(self.match_builder(profile, exclusive=exclusive, negate=False))
			if len(includes) > 1:
				includes = [ "(" + x + ")" if " AND " in x else x for x in includes ]
			clause.append(" OR ".join(includes))
			if len(includes) > 1:
				clause[-1] = "(" + clause[-1] + ")"

		# adding conditions of profiles to exclude to where clause
		if exclude_profiles:
			excludes = []
			for profile in exclude_profiles:
				excludes.append(self.match_builder(profile, exclusive=exclusive, negate=True))
			excludes = [ "(" + x + ")" if " AND " in x else x for x in excludes ]
			clause.append(" AND ".join(excludes))

		# adding accession, lineage, zips, and dates based conditions
		if accessions:
			include_acc = [x for x in accessions if not x.startswith("^")]
			exclude_acc = [x for x in accessions if x.startswith("^")]
			if include_acc:
				clause.append("accession IN (" + " , ".join(['?'] * len(include_acc)) + ")")
				vals.extend(include_acc)
			if exclude_acc:
				clause.append("accession NOT IN (" + " , ".join(['?'] * len(exclude_acc)) + ")")
				vals.extend(exclude_acc)
		if lineages:
			include_lin = [x for x in lineages if not x.startswith("^")]
			exclude_lin = [x for x in lineages if x.startswith("^")]
			if include_lin:
				clause.append("lineage IN (" + " ,".join(['?'] * len(include_lin))  + ")")
				vals.extend(include_lin)
			if exclude_acc:
				clause.append("lineage NOT IN (" + " ,".join(['?'] * len(exclude_acc))  + ")")
				vals.extend(exclude_acc)
		if zips:
			include_zip = [x for x in zips if not str(x).startswith("^")]
			exclude_zip = [x for x in zips if str(x).startswith("^")]
			if include_zip:
				z = []
				for zp in include_zip:
					z.append("zip LIKE '" + str(zp) + "%'")
				clause.append(" OR ".join(z))
				if len(z) > 1:
					clause[-1] = "(" + clause[-1] + ")"
			if exclude_zip:
				z = []
				for zp in exclude_zip:
					z.append("zip NOT LIKE '" + str(zp) + "%'")
				clause.append(" AND ".join(z))
		if dates:
			d = []
			for dt in dates:
				if ":" in dt:
					x, y = dt.split(":")
					d.append("(date BETWEEN '" + x + "' AND '" + y + "')")
				else:
					d.append("date = " + dt)
			clause.append(" OR ".join(d))
			if len(d) > 1:
				clause[-1] = "(" + clause[-1] + ")"

		# executing the query and storing results
		clause = " AND ".join(clause)
		with sonarDBManager(self.db, readonly=True) as dbm:
			rows = [x for x in dbm.select("essence", whereClause=clause, whereVals=vals)]

		# remove ambiguities from database profiles if wished
		if not ambig:
			keep = [item for sublist in include_profiles for item in sublist] if include_profiles else None
			for i in range(len(rows)):
				rows[i]['dna_profile'] = self.filter_ambig(rows[i]['dna_profile'], self.iupac_explicit_nt_code, keep)
				rows[i]['aa_profile'] = self.filter_ambig(rows[i]['aa_profile'], self.iupac_explicit_aa_code, keep)
		return rows

	# VALIDATION

	def restore_genome(self, acc, dbm):
		"""
		function to restore a genome sequence from the SONAR database

		Parameters
		----------

		acc : str
			define the accesion of the genome that should be restored
		dbm : object
			define a sonarDBManager handling the database connection

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
		rows = self.get_dna_vars(acc, dbm)
		if rows:
			qryseq = list(self.refseq)
			for row in rows:
				if row['start'] is not None:
					s = row['start']
					if s >= 0:
						if row['ref'] != self.refseq[s]:
							sys.exit("data error: data inconsistency found for '" + acc + "' (" + row['ref']+ " expected at position " + str(s+1) + " of the reference sequence, got " + refseq[s] + ").")
						qryseq[s] = "" if not row['alt'] else row['alt']
					else:
						qryseq = [row['alt']] + qryseq
			return (">" + rows[0]['description'], "".join(qryseq))
		return None

	def restore_alignment(self, acc, dbm):
		"""
		function to restore a genome alignment from the SONAR database

		Parameters
		----------

		acc : str
			define the accesion of the genome whose alignment versus the reference
			should be restored
		dbm : object
			define a sonarDBManager handling the database connection

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
		rows = self.get_dna_vars(acc, dbm)
		if rows:
			refseq = list(self.refseq)
			qryseq = refseq[:]
			for row in rows:
				if row['start'] is not None:
					s = row['start']
					if s >= 0:
						if row['ref'] != self.refseq[s]:
							sys.exit("data error: data inconsistency found for '" + acc + "' (" + row['ref']+ " expected at position " + str(s+1) + " of the reference sequence, got " + refseq[s] + ").")
						qryseq[s] = "-" if not row['alt'] else row['alt']
						if len(row['alt']) > 1:
							refseq[s] +=  "-" * (len(row['alt'])-1)
					else:
						qryseq = [row['alt']] + qryseq
						refseq = ["-" * (len(row['alt']))] + refseq
			return  (">" + rows[0]['description'], "".join(qryseq), ">" + self.dbobj.refdescr, "".join(refseq))
		return None

	def be_paranoid(self, acc, orig_seq, dbm, auto_delete=False):
		"""
		function to compare a given sequence with the respective sequence restored
		from the SONAR database

		Parameters
		----------

		acc : str
			define the accesion of the genome that should be validated
		orig_seq : str
			define the sequence expected
		dbm : object
			define a sonarDBManager handling the database connection
		auto_delete : bool
			define if the respective genome should be automatically deleted
			from the SONAR database if the test fails

		Returns
		-------

		bool
			True is returned if expected and restored sequences are not different
			otherwise False
		"""
		orig_seq = orig_seq.upper()
		restored_seq = self.restore_genome(acc, dbm=dbm)[1]
		if orig_seq != restored_seq:
			rows = self.get_dna_vars(acc, dbm)
			writer = csv.DictWriter(sys.stdout, rows[0].keys(), lineterminator=os.linesep)
			writer.writeheader()
			writer.writerows(rows)
			if auto_delete:
				self.delete_accession(acc, dbm)
			print("Good that you are paranoid: " + acc + " original and those restored from the database do not match.", file=sys.stderr)
			return False
		return True

	# OTHER

	@staticmethod
	def get_version():
		with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), ".version"), "r") as handle:
			return handle.read().strip()

class sonarCache():
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
		self._idx = None
		self.cache = defaultdict(set)

		if self.temp:
			self.dirname = mkdtemp(prefix=".sonarCache_")
		else:
			self.dirname = os.path.abspath(dir)
			self.checkdir()

	def __enter__(self):
		return self

	def __exit__(self, exc_type, exc_value, exc_traceback):
		self.write_idx(backup=True)
		if [exc_type, exc_value, exc_traceback].count(None) != 3:
			print("warning:", file=sys.stderr)
			print(traceback.format_exc(), file=sys.stderr)
		if os.path.isdir(self.dirname) and self.temp:
			shutil.rmtree(self.dirname)

	@property
	def idx(self):
		if self._idx is None:
			self._idx = os.path.join(self.dirname, "cache.idx")
		return self._idx

	def checkdir(self):
		if not os.path.isdir(self.dirname):
			os.makedirs(self.dirname)
		elif os.listdir(self.dirname):
			self.restore_cache()

	def restore_cache(self):
		if not os.path.isfile(self.idx):
			sys.exit("chache error: index file missing.")
		self.cache = self.load_idx()
		for seqhash in self.cache:
			files = self.get_cache_files(seqhash)
			if not os.path.isfile(files[0]) or not os.path.isfile(files[1]) :
				sys.exit("cache error: chache is corrupted.")

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
		return base64.urlsafe_b64encode(string.encode('UTF-8') ).decode('UTF-8')

	@staticmethod
	def deslugify(string):
		return base64.urlsafe_b64decode(string).decode("utf-8")

	@staticmethod
	def seqhash_from_fasta_name(fname):
		return sonarCache.deslugify(os.path.basename(fname))[:-5]

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

	def read_fasta(self, fname):
		record = SeqIO.read(fname, "fasta")
		return record.id, record.description[1:], str(record.seq).upper()

	def get_fasta_fname(self, seqhash):
		return os.path.join(self.dirname, self.slugify(seqhash) + ".fasta")

	def get_algn_fname(self, seqhash):
		return os.path.join(self.dirname, self.slugify(seqhash) + ".algn")

	def get_info_fname(self, seqhash):
		return os.path.join(self.dirname, self.slugify(seqhash) + ".pickle")

	def get_cache_files(self, seqhash):
		return self.get_fasta_fname(seqhash), self.get_algn_fname(seqhash), self.get_info_fname(seqhash)

	def load_info(self, seqhash):
		with open(self.get_info_fname(seqhash), 'rb') as handle:
			return pickle.load(handle, encoding="bytes")

	def write_info(self, seqhash, data):
		with open(self.get_info_fname(seqhash), 'wb') as handle:
			pickle.dump(data, handle)

	def load_idx(self):
		with open(self.idx, 'rb') as handle:
			return pickle.load(handle, encoding="bytes")

	def write_idx(self, backup):
		if backup and os.path.isfile(self.idx):
			shutil.copy(self.idx, self.idx + ".old")
		with open(self.idx, 'wb') as handle:
			pickle.dump(self.cache, handle)

	def add_fasta(self, fname):
		"""
		function to cache genomes from a valid FASTA file

		Parameters
		----------

		fname : str
			define the path to a valid FASTA file
		"""
		for record in SeqIO.parse(fname, "fasta"):
			acc = record.id
			descr = record.description
			seq = str(record.seq).upper()
			seqhash = sonarDB.hash(seq)
			fasta, align, info = self.get_cache_files(seqhash)

			#check or create cached fasta file
			if os.path.isfile(fasta):
				if seq != self.get_cached_seq(seqhash):
					sys.exit("cache error: sequence collision for hash '" + seqhash + "'.")
			else:
				with open(fasta, "w") as handle:
					handle.write(">" + seqhash + os.linesep + seq)

			# check for accession collision
			entry = (acc, descr)
			found_seqhashes = []
			for key in self.cache:
				if acc in [x[0] for x in self.cache[key]]:
					found_seqhashes.append(key)
			if found_seqhashes:
				if len(found_seqhashes) > 1 or found_seqhashes[0] != seqhash:
					sys.exit("cache error: accession collision for '" + acc + "'.")
			else:
				self.cache[seqhash].add(entry)

	def get_cached_seqhashes(self):
		return set(self.cache.keys())

	def iter_cached_fasta_files(self):
		for x in self.cache:
			yield self.get_fasta_fname(x)

	def get_cached_seq(self, seqhash):
		return self.read_fasta(self.get_fasta_fname(seqhash))[-1]

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
