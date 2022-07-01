#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

# DEPENDENCIES
import collections
import csv
import gzip
import logging
import lzma
import os
import re
import sys
from tempfile import mkstemp

from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid
import pandas as pd

from . import __version__
from .dbm import sonarDBManager


class Aliasor:
    def __init__(self, alias_file):
        aliases = pd.read_json(alias_file)
        self.alias_dict = {
            x: x if x.startswith("X") else aliases[x][0] for x in aliases.columns
        }
        self.alias_dict["A"] = "A"
        self.alias_dict["B"] = "B"
        self.realias_dict = {v: k for k, v in self.alias_dict.items()}

    def compress(self, name):
        name_split = name.split(".")
        if len(name_split) < 5:
            return name
        letter = self.realias_dict[".".join(name_split[0:4])]
        if len(name_split) == 5:
            return letter + "." + name_split[4]
        else:
            return letter + "." + ".".join(name_split[4:])

    def uncompress(self, name):
        name_split = name.split(".")
        letter = name_split[0]
        unaliased = self.alias_dict[letter]
        if len(name_split) == 1:
            return name
        if len(name_split) == 2:
            return unaliased + "." + name_split[1]
        else:
            return unaliased + "." + ".".join(name_split[1:])


# CLASS
class sonarBasics(object):
    """
    this object provides sonarBasics functionalities and intelligence

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

    >>> db = sonarBasics(DOCTESTDB)

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

    def __init__(self, dbfile, translation_table=1, debug=False):
        self.db = os.path.abspath(dbfile) if dbfile else mkstemp()[1]
        self.debug = debug
        self.__moduledir = self.get_module_base()
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
        self.__fasta_tag_regex = None
        self._covSonar_version = self.get_version()
        logging.basicConfig(format="%(asctime)s %(message)s")

    @staticmethod
    def get_version():
        return __version__

    # DB MAINTENANCE

    @staticmethod
    def get_module_base(*join_with):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), *join_with)

    @staticmethod
    def setup_db(  # noqa: C901
        fname, auto_create=False, default_setup=True, reference_gb=None, debug=False
    ):
        if os.path.isfile(fname):
            sys.exit("setup error: " + fname + " does already exist.")
        sonarDBManager.setup(fname, debug=debug)

        # loading default data
        if default_setup or auto_create:
            with sonarDBManager(fname, readonly=False, debug=debug) as dbm:
                # adding pre-defined sample properties
                dbm.add_property(
                    "imported",
                    "date",
                    "date",
                    "date sample has been imported to the database",
                )
                dbm.add_property(
                    "modified",
                    "date",
                    "date",
                    "date when sample data has been modified lastly",
                )
                dbm.add_property(
                    "nuc_profile",
                    "text",
                    "text",
                    "stores the nucleotide level profiles",
                )
                dbm.add_property(
                    "aa_profile", "text", "text", "stores the aa level profiles"
                )
                dbm.add_property(
                    "nuc_n_profile",
                    "text",
                    "text",
                    "stores the nucleotide(with N) level profiles",
                )
                dbm.add_property(
                    "aa_x_profile",
                    "text",
                    "text",
                    "stores the aa(with X) level profiles",
                )

                # if enable, create PREDEFINED properties
                if auto_create:
                    dbm.add_property("DATE_DRAW", "date", "date", "Sampling date")
                    dbm.add_property(
                        "PROCESSING_DATE", "date", "date", "Submission/Processing date"
                    )
                    dbm.add_property(
                        "COUNTRY", "text", "text", "Country where a sample belongs to"
                    )
                    dbm.add_property("HOST", "text", "text", "e.g., HUMAN")
                    dbm.add_property("ZIP", "text", "text", "zip code e.g., 33602")
                    dbm.add_property("LAB", "text", "text", "lab id e.g., 11069")
                    dbm.add_property("LINEAGE", "text", "text", "e.g., BA.2 or B.1.1.7")
                    dbm.add_property("TECHNOLOGY", "text", "text", "e.g., ILLUMINA")

                # adding reference
                if not reference_gb:
                    reference_gb = os.path.join(
                        os.path.dirname(os.path.abspath(__file__)), "data/ref.gb"
                    )
                records = [x for x in sonarBasics.iter_genbank(reference_gb)]
                ref_id = dbm.add_reference(
                    records[0]["accession"],
                    records[0]["description"],
                    records[0]["organism"],
                    1,
                    1,
                )

                # adding reference molecule and elements
                for i, record in enumerate(records):
                    gene_ids = {}
                    s = 1 if i == 0 else 0

                    mol_id = dbm.insert_molecule(
                        ref_id,
                        record["moltype"],
                        record["accession"],
                        record["symbol"],
                        record["description"],
                        i,
                        record["length"],
                        s,
                    )

                    # source handling
                    source_id = dbm.insert_element(
                        mol_id,
                        "source",
                        record["source"]["accession"],
                        record["source"]["symbol"],
                        record["source"]["description"],
                        record["source"]["start"],
                        record["source"]["end"],
                        record["source"]["strand"],
                        record["source"]["sequence"],
                        standard=1,
                        parts=record["source"]["parts"],
                    )
                    if record["source"]["sequence"] != dbm.get_sequence(source_id):
                        sys.exit(
                            "genbank error: could not recover sequence of '"
                            + record["source"]["accession"]
                            + "' (source)"
                        )

                    # gene handling
                    for elem in record["gene"]:
                        gene_ids[elem["accession"]] = dbm.insert_element(
                            mol_id,
                            "gene",
                            elem["accession"],
                            elem["symbol"],
                            elem["description"],
                            elem["start"],
                            elem["end"],
                            elem["strand"],
                            elem["sequence"],
                            standard=0,
                            parent_id=source_id,
                            parts=elem["parts"],
                        )
                        if elem["sequence"] != dbm.extract_sequence(
                            gene_ids[elem["accession"]]
                        ):
                            sys.exit(
                                "genbank error: could not recover sequence of '"
                                + elem["accession"]
                                + "' (gene)"
                            )

                    # cds handling
                    for elem in record["cds"]:
                        cid = dbm.insert_element(
                            mol_id,
                            "cds",
                            elem["accession"],
                            elem["symbol"],
                            elem["description"],
                            elem["start"],
                            elem["end"],
                            elem["strand"],
                            elem["sequence"],
                            0,
                            gene_ids[elem["gene"]],
                            elem["parts"],
                        )
                        if elem["sequence"] != dbm.extract_sequence(
                            cid, translation_table=1
                        ):
                            sys.exit(
                                "genbank error: could not recover sequence of '"
                                + elem["accession"]
                                + "' (cds)"
                            )
                print("Success: Database was successfully installed")

    # DATA IMPORT
    # genbank handling handling
    @staticmethod
    def process_segments(feat_location_parts, cds=False):
        base = 0
        div = 1 if not cds else 3
        segments = []
        for i, segment in enumerate(feat_location_parts, 1):
            segments.append(
                [int(segment.start), int(segment.end), segment.strand, base, i]
            )
            base += round((segment.end - segment.start - 1) / div, 1)
        return segments

    @staticmethod
    def iter_genbank(fname):
        gb_data = {}
        for gb_record in SeqIO.parse(fname, "genbank"):
            # adding general annotation
            gb_data["accession"] = (
                gb_record.name + "." + str(gb_record.annotations["sequence_version"])
            )
            gb_data["symbol"] = (
                gb_record.annotations["symbol"]
                if "symbol" in gb_record.annotations
                else ""
            )
            gb_data["organism"] = gb_record.annotations["organism"]
            gb_data["moltype"] = None
            gb_data["description"] = gb_record.description
            gb_data["length"] = None
            gb_data["segment"] = ""
            gb_data["gene"] = []
            gb_data["cds"] = []
            gb_data["source"] = None

            # adding source annotation
            source = [x for x in gb_record.features if x.type == "source"]
            if len(source) != 1:
                sys.exit(
                    "genbank error: expecting exactly one source feature (got: "
                    + str(len(source))
                    + ")"
                )
            feat = source[0]
            gb_data["moltype"] = (
                feat.qualifiers["mol_type"][0] if "mol_type" in feat.qualifiers else ""
            )
            gb_data["source"] = {
                "accession": gb_data["accession"],
                "symbol": gb_data["accession"],
                "start": int(feat.location.start),
                "end": int(feat.location.end),
                "strand": "",
                "sequence": sonarBasics.harmonize(feat.extract(gb_record.seq)),
                "description": "",
                "parts": sonarBasics.process_segments(feat.location.parts),
            }
            gb_data["length"] = len(gb_data["source"]["sequence"])
            if "segment" in feat.qualifiers:
                gb_data["segment"] = feat.qualifiers["segment"][0]

            for feat in gb_record.features:
                # adding gene annotation
                if feat.type == "gene":
                    gb_data["gene"].append(
                        {
                            "accession": feat.id
                            if feat.id != "<unknown id>"
                            else feat.qualifiers["gene"][0],
                            "symbol": feat.qualifiers["gene"][0],
                            "start": int(feat.location.start),
                            "end": int(feat.location.end),
                            "strand": feat.strand,
                            "sequence": sonarBasics.harmonize(
                                feat.extract(gb_data["source"]["sequence"])
                            ),
                            "description": "",
                            "parts": sonarBasics.process_segments(feat.location.parts),
                        }
                    )

                # adding cds annotation
                elif feat.type == "CDS":
                    gb_data["cds"].append(
                        {
                            "accession": feat.qualifiers["protein_id"][0],
                            "symbol": feat.qualifiers["gene"][0],
                            "start": int(feat.location.start),
                            "end": int(feat.location.end),
                            "strand": feat.strand,
                            "gene": feat.qualifiers["gene"][0],
                            "sequence": feat.qualifiers["translation"][0],
                            "description": feat.qualifiers["product"][0],
                            "parts": sonarBasics.process_segments(
                                feat.location.parts, True
                            ),
                        }
                    )
            yield gb_data

    # seq handling
    @staticmethod
    def hash(seq):
        """ """
        return seguid(seq)

    @staticmethod
    def harmonize(seq):
        """ """
        return str(seq).strip().upper().replace("U", "T")

    # matching
    @staticmethod
    def match(
        db,
        profiles=[],
        reserved_props_dict={},
        propdict={},
        reference=None,
        outfile=None,
        format="csv",
        debug="False",
    ):
        with sonarDBManager(db, debug=debug) as dbm:
            if format == "vcf" and reference is None:
                reference = dbm.get_default_reference_accession()
            cursor = dbm.match(
                *profiles,
                reserved_props=reserved_props_dict,
                properties=propdict,
                reference_accession=reference,
                format=format,
            )
            if format == "csv" or format == "tsv":
                tsv = True if format == "tsv" else False
                sonarBasics.exportCSV(
                    cursor, outfile=outfile, na="*** no match ***", tsv=tsv
                )
            elif format == "count":
                if outfile:
                    with open(outfile, "w") as handle:
                        handle.write(str(cursor.fetchone()["count"]))
                else:
                    print(cursor.fetchone()["count"])
            elif format == "vcf":
                sonarBasics.exportVCF(
                    cursor, reference=reference, outfile=outfile, na="*** no match ***"
                )
            else:
                sys.exit("error: '" + format + "' is not a valid output format")

    # restore
    @staticmethod
    def restore(
        db, *samples, reference_accession=None, aligned=False, outfile=None, debug=False
    ):
        with sonarDBManager(db, readonly=True, debug=debug) as dbm:
            handle = sys.stdout if outfile is None else open(outfile, "w")
            for sample in samples:
                prefixes = collections.defaultdict(str)
                molecules = {
                    x["element.id"]: {
                        "seq": list(x["element.sequence"]),
                        "mol": x["element.symbol"],
                    }
                    for x in dbm.get_alignment_data(sample, reference_accession=None)
                }
                gap = "-" if aligned else ""
                for vardata in dbm.iter_dna_variants(sample, *molecules.keys()):
                    if aligned and len(vardata["variant.alt"]) > 1:
                        vardata["variant.alt"] = (
                            vardata["variant.alt"][0]
                            + vardata["variant.alt"][1:].lower()
                        )
                    if vardata["variant.alt"] == " ":
                        for i in range(
                            vardata["variant.start"], vardata["variant.end"]
                        ):
                            molecules[vardata["element.id"]]["seq"][i] = gap
                    elif vardata["variant.start"] >= 0:
                        molecules[vardata["element.id"]]["seq"][
                            vardata["variant.start"]
                        ] = vardata["variant.alt"]
                    else:
                        prefixes[vardata["element.id"]] = vardata["variant.alt"]
                molecules_len = len(molecules)
                records = []
                for element_id in molecules:
                    if molecules_len == 1:
                        records.append(">" + sample)
                    else:
                        records.append(
                            ">"
                            + sample
                            + " [molecule="
                            + molecules[vardata["element.id"]]["mol"]
                            + "]"
                        )
                    records.append(
                        prefixes[element_id]
                        + "".join(molecules[vardata["element.id"]]["seq"])
                    )
                if len(records) > 0:
                    handle.write("\n".join(records) + "\n")

    # restore
    @staticmethod
    def delete(db, *samples, debug):
        with sonarDBManager(db, readonly=False, debug=debug) as dbm:
            before = dbm.count_samples()
            dbm.delete_samples(*samples)
            after = dbm.count_samples()
            print(
                before - after,
                "of",
                len(samples),
                "samples found and deleted.",
                after,
                "samples remain in the database.",
            )

    @staticmethod
    def show_db_info(db):
        with sonarDBManager(db, readonly=True) as dbm:
            print("covSonar Version:          ", sonarBasics.get_version())
            print("database path:             ", dbm.dbfile)
            print("database version:          ", dbm.get_db_version())
            print("database size:             ", dbm.get_db_size())
            print("unique sequences:          ", dbm.count_sequences())
            # print("Sample properties          ", dbm.get_earliest_import())
            # print("latest genome import:      ", dbm.get_latest_import())
            # print("earliest sampling date:    ", dbm.get_earliest_date())
            # print("latest sampling date:      ", dbm.get_latest_date())

    # output
    # profile generation
    @staticmethod
    def iter_formatted_match(cursor):
        nuc_profiles = collections.defaultdict(list)
        aa_profiles = collections.defaultdict(list)
        samples = set()
        logging.warning("getting profile data")
        logging.warning("processing profile data")
        for row in cursor:

            samples.add(row["samples"])
            if row["element.type"] == "cds":
                aa_profiles[row["samples"]].append(
                    (
                        row["element.id"],
                        row["variant.start"],
                        row["element.symbol"] + ":" + row["variant.label"],
                    )
                )
            else:
                nuc_profiles[row["samples"]].append(
                    (row["element.id"], row["variant.start"], row["variant.label"])
                )
        out = []
        print(len(samples))
        for sample in sorted(samples):
            print(sorted(nuc_profiles[sample], key=lambda x: (x[0], x[1])))

        logging.warning("assembling profile data")
        for sample in sorted(samples):
            print(sample)
            out.append(
                {
                    "sample.name": sample,
                    "nuc_profile": " ".join(
                        sorted(nuc_profiles[sample], key=lambda x: (x[0], x[1]))
                    ),
                    "aa_profile": " ".join(
                        sorted(aa_profiles[sample], key=lambda x: (x[0], x[1]))
                    ),
                }
            )
        print(out)
        return out

    # csv
    def exportCSV(cursor, outfile=None, na="*** no data ***", tsv=False):
        i = -1
        try:
            for i, row in enumerate(cursor):
                if i == 0:
                    outfile = sys.stdout if outfile is None else open(outfile, "w")
                    sep = "\t" if tsv else ","
                    writer = csv.DictWriter(
                        outfile, row.keys(), delimiter=sep, lineterminator=os.linesep
                    )  # extrasaction='ignore',
                    writer.writeheader()
                writer.writerow(row)
            if i == -1:
                print(na)
        except Exception:
            print("An exception occurred", row)
            raise

    # vcf
    def exportVCF(cursor, reference, outfile=None, na="*** no match ***"):  # noqa: C901
        records = collections.OrderedDict()
        all_samples = set()
        for row in cursor.fetchall():  # sonarBasics.iter_formatted_match(cursor):
            chrom, pos, ref, alt, samples = (
                row["molecule.accession"],
                row["variant.start"],
                row["variant.ref"],
                row["variant.alt"],
                row["samples"],
            )
            print(samples)
            if chrom not in records:
                records[chrom] = collections.OrderedDict()
            if pos not in records[chrom]:
                records[chrom][pos] = {}
            if ref not in records[chrom][pos]:
                records[chrom][pos][ref] = {}
            records[chrom][pos][ref][alt] = set(samples.split("\t"))
            all_samples.update(samples.split("\t"))
        print(all_samples)
        if len(records) != 0:
            all_samples = sorted(all_samples)
            if outfile is None:
                handle = sys.stdout
            else:
                # if not outfile.endswith(".gz"):
                # outfile += ".gz"
                handle = open(outfile, mode="w")  # bgzf.open(outfile, mode='wb')

            # vcf header
            handle.write("##fileformat=VCFv4.2\n")
            handle.write("##poweredby=covSonar2\n")
            handle.write("##reference=" + reference + "\n")
            handle.write(
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype"\n'
            )
            handle.write(
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                + "\t".join(all_samples)
                + "\n"
            )

            # records
            for chrom in records:
                for pos in records[chrom]:
                    for ref in records[chrom][pos]:
                        # snps and inserts (combined output)
                        alts = [x for x in records[chrom][pos][ref].keys() if x.strip()]
                        if alts:
                            alt_samples = set()
                            gts = []
                            for alt in alts:
                                samples = records[chrom][pos][ref][alt]
                                alt_samples.update(samples)
                                gts.append(
                                    ["1" if x in samples else "0" for x in all_samples]
                                )
                            gts = [
                                ["0" if x in alt_samples else "1" for x in all_samples]
                            ] + gts
                            record = [
                                chrom,
                                str(pos),
                                ref,
                                ",".join(alts),
                                ".",
                                ".",
                                ".",
                                "GT",
                            ] + ["/".join(x) for x in zip(*gts)]
                        # dels (individual output)
                        for alt in [
                            x for x in records[chrom][pos][ref].keys() if not x.strip()
                        ]:
                            samples = records[chrom][pos][ref][alt]
                            record = [
                                chrom,
                                str(pos),
                                ref,
                                ref[0],
                                ".",
                                ".",
                                ".",
                                "GT",
                            ] + ["0/1" if x in samples else "1/0" for x in all_samples]
                        handle.write("\t".join(record) + "\n")
            handle.close()
        else:
            print(na)

    # gerenal
    @staticmethod
    def open_file(fname, mode="r", compressed=False, encoding=None):
        if not os.path.isfile(fname):
            sys.exit("input error: " + fname + " does not exist.")
        if compressed == "auto":
            compressed = os.path.splitext(fname)[1][1:]
        try:
            if compressed == "gz":
                return gzip.open(fname, mode + "t", encoding=encoding)
            if compressed == "xz":
                return lzma.open(fname, mode + "t", encoding=encoding)
            else:
                return open(fname, mode, encoding=encoding)
        except Exception:
            sys.exit("input error: " + fname + " cannot be opened.")

    @staticmethod
    def set_key(dictionary, key, value):
        if key not in dictionary:
            dictionary[key] = value
        elif type(dictionary[key]) == list:
            dictionary[key].append(value)
        else:
            dictionary[key] = [dictionary[key], value]
        return dictionary
