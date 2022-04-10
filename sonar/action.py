#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

# DEPENDENCIES
import os
import re
import sys
from Bio.SeqUtils.CheckSum import seguid
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Emboss.Applications import StretcherCommandline
import pickle
from tqdm import tqdm
import tempfile
import itertools
from contextlib import ExitStack
from sonar import sonarDBManager, sonarAligner, sonartoVCF, sonartoVCF_v2
import csv
from tabulate import tabulate
from textwrap import fill
from tempfile import mkstemp, mkdtemp
import gzip
import lzma


# CLASS
class sonarActions(object):
	"""
	this object provides sonarActions functionalities and intelligence

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

	>>> db = sonarActions(DOCTESTDB)

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
	def __init__(self, dbfile, translation_table = 1, debug=False):
		self.db = os.path.abspath(dbfile)  if dbfile else mkstemp()[1]
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

	@staticmethod
	def get_version():
		with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), ".version"), "r") as handle:
			return handle.read().strip()

	# DB MAINTENANCE

	@staticmethod
	def get_module_base(*join_with):
		return os.path.join(os.path.dirname(os.path.realpath(__file__)), *join_with)

	@staticmethod
	def setup_db(fname, auto_create=False, default_setup=True, debug=False):
		if os.path.isfile(fname):
			sys.exit("setup error: " + fname + " does already exist.")
		sonarDBManager.setup(fname, debug=debug)
		## loading default data
		if default_setup:
			with sonarDBManager(fname, debug=debug) as dbm:
				### adding pre-defined sample properties
				dbm.add_property("imported", "date", "date", "date sample has been imported to the database")
				dbm.add_property("modified", "date", "date", "date when sample data has been modified lastly")
				##  if enable, create important properties
				if (auto_create):
					print('Enable automatic import property')
					dbm.add_property("DATE_DRAW", "date", "date", "Sampling date")
					dbm.add_property("PROCESSING_DATE ", "date", "date", "Submission date")
					dbm.add_property("country", "text", "text", "Country where a sample belongs to")
					dbm.add_property("host", "text", "text", "e.g., HUMAN")
					dbm.add_property("collection", "text", "text", "e.g., DESH_STICHPROBE")
					dbm.add_property("zip", "text", "text", "e.g., 33602")
					dbm.add_property("lab", "text", "text", "e.g., 11069")
					dbm.add_property("lineage", "text", "text", "e.g., BA.2 or B.1.1.7")
					dbm.add_property("technology", "text", "text", "e.g., ILLUMINA")

				### adding built-in reference (unsegmented genome)
				records = [x for x in sonarActions.iter_genbank(sonarActions.get_module_base("ref.gb"))]
				ref_id = dbm.add_reference(records[0]['accession'], records[0]['description'], records[0]['organism'], 1, 1)

				### adding reference molecule and elements
				for i, record in enumerate(records):
					gene_ids = {}
					s = 1 if i == 0 else 0

					mol_id = dbm.insert_molecule(ref_id, record['moltype'], record['accession'], record['symbol'], record['description'], i, record['length'], s)

					#### source handling
					source_id = dbm.insert_element(mol_id, "source", record['source']['accession'], record['source']['symbol'], record['source']['description'], record['source']['start'], record['source']['end'], record['source']['strand'], record['source']['sequence'], standard=1, parts=record['source']['parts'])
					if record['source']['sequence'] != dbm.get_sequence(source_id):
						sys.exit("genbank error: could not recover sequence of '" + record['source']['accession'] + "' (source)")

					#### gene handling
					for elem in record['gene']:
						gene_ids[elem['accession']] = dbm.insert_element(mol_id, "gene", elem['accession'], elem['symbol'], elem['description'], elem['start'], elem['end'], elem['strand'], elem['sequence'], standard=0, parent_id=source_id, parts=elem['parts'])
						if elem['sequence'] != dbm.extract_sequence(gene_ids[elem['accession']]):
							sys.exit("genbank error: could not recover sequence of '" + elem['accession'] + "' (gene)")

					#### cds handling
					for elem in record['cds']:
						cid = dbm.insert_element(mol_id, "cds", elem['accession'], elem['symbol'], elem['description'], elem['start'], elem['end'], elem['strand'], elem['sequence'], 0, gene_ids[elem['gene']], elem['parts'])
						if elem['sequence'] != dbm.extract_sequence(cid, translation_table=1):
							sys.exit("genbank error: could not recover sequence of '" + elem['accession'] + "' (cds)")
	# DATA IMPORT

	## genbank handling handling
	@staticmethod
	def process_segments(feat_location_parts, cds = False):
		base = 0
		div = 1 if not cds else 3
		segments = []
		for i, segment in enumerate(feat_location_parts, 1):
			segments.append([int(segment.start), int(segment.end), segment.strand, base,i])
			base += round((segment.end-segment.start-1)/div, 1)
		return segments

	@staticmethod
	def iter_genbank(fname):
		gb_data = {}
		for gb_record in SeqIO.parse(fname, "genbank"):

			## adding general annotation
			gb_data['accession'] = gb_record.name + "." + str(gb_record.annotations['sequence_version'])
			gb_data['symbol'] = gb_record.annotations['symbol'] if 'symbol' in gb_record.annotations else ""
			gb_data['organism'] = gb_record.annotations['organism']
			gb_data['moltype'] = None
			gb_data['description'] = gb_record.description
			gb_data['length'] = None
			gb_data['segment'] = ""
			gb_data['gene'] = []
			gb_data['cds'] = []
			gb_data['source'] = None

			## adding source annotation
			source = [x for x in gb_record.features if x.type == "source"]
			if len(source) !=  1:
				sys.exit("genbank error: expecting exactly one source feature (got: " + str(len(source)) + ")")
			feat = source[0]
			gb_data['moltype'] = feat.qualifiers["mol_type"][0] if "mol_type" in feat.qualifiers else ""
			gb_data['source'] = {
				"accession": gb_data['accession'],
				"symbol": gb_data['accession'],
				"start": int(feat.location.start),
				"end": int(feat.location.end),
				"strand": "",
				"sequence": sonarActions.harmonize(feat.extract(gb_record.seq)),
				"description": "",
				"parts": sonarActions.process_segments(feat.location.parts)
			}
			gb_data['length'] = len(gb_data['source']['sequence'])
			if "segment" in feat.qualifiers:
				gb_data['segment'] = feat.qualifiers["segment"][0]

			for feat in gb_record.features:
				## adding gene annotation
				if feat.type == "gene":
					gb_data['gene'].append({
						"accession": feat.id if feat.id != "<unknown id>" else feat.qualifiers['gene'][0],
						"symbol": feat.qualifiers['gene'][0],
						"start": int(feat.location.start),
						"end": int(feat.location.end),
						"strand": feat.strand,
						"sequence": sonarActions.harmonize(feat.extract(gb_data['source']['sequence'])),
						"description": "",
						"parts": sonarActions.process_segments(feat.location.parts)
					})


				## adding cds annotation
				elif feat.type == "CDS":
					gb_data['cds'].append({
						"accession": feat.qualifiers['protein_id'][0],
						"symbol": feat.qualifiers['gene'][0],
						"start": int(feat.location.start),
						"end": int(feat.location.end),
						"strand": feat.strand,
						"gene": feat.qualifiers['gene'][0],
						"sequence": feat.qualifiers['translation'][0],
						"description": feat.qualifiers['product'][0],
						"parts": sonarActions.process_segments(feat.location.parts, True)
					})
			yield gb_data


	## seq handling
	@staticmethod
	def hash(seq):
		"""
		"""
		return seguid(seq)

	@staticmethod
	def harmonize(seq):
		"""
		"""
		return str(seq).strip().upper().replace("U", "T")

	# matching
	def match(self, profiles=[], propdict={}, count=None):
		with sonarDBManager(self.db, debug=self.debug) as dbm:
			#  get ID
			rows = dbm.match(*profiles, properties=propdict)


			if count:
				print(len(rows))
			else:
				# print(rows)

				self.rows_to_csv(rows, na="*** no match ***", tsv=True)
	
	# show db infos
	def show_props(self):
		with sonarDBManager(self.db, debug=self.debug) as dbm:
			cols = ["name", "argument", "description", "data type", "query type", "default value"]
			rows = []
			if not dbm.properties:
				print("the database does not conatin any sample property.")
				exit()
			for prop in sorted(dbm.properties.keys()):
				rows.append([])
				dt = dbm.properties[prop]['datatype'] if dbm.properties[prop]['datatype'] != "float" else "floating point"
				rows[-1].append(prop)
				rows[-1].append("--" + prop)
				rows[-1].append(fill(dbm.properties[prop]['description'], width=25))
				rows[-1].append(dt)
				rows[-1].append(dbm.properties[prop]['querytype'])
				rows[-1].append(dbm.properties[prop]['standard'])
			print(tabulate(rows, headers=cols, tablefmt='fancy_grid'))
			print()
			print("DATE FORMAT")
			print("dates must comply with the following format: YYYY-MM-DD")
			print()
			print("OPERATORS")
			print("integer, floating point and date data types support the following operators prefixed directly to the respective value without spaces:")
			print("  > larger than (e.g. >1)")
			print("  < smaller than (e.g. <1)")
			print("  >= larger than or equal to (e.g. >=1)")
			print("  <= smaller than or equal to (e.g. <=1)")
			print("  != different than (e.g. !=1)")
			print()
			print("RANGES")
			print("integer, floating point and date data types support ranges defined by two values directly connected by a colon (:) with no space between them:")
			print("  e.g. 1:10 (between 1 and 10)")
			print("  e.g. 2021-01-01:2021-12-31 (between 1st Jan and 31st Dec of 2021)")
			print()
	
	def show_system_info(self):
		print("System info.")
		print("conSonar Version:      ", self._covSonar_version)
		#print("reference length:      ", str(len(self.db.refseq)) + "bp")
		#print("annotated proteins:    ", ", ".join(self.db.refgffObj.symbols))
		#print("used translation table:", self.db.translation_table)

	def show_db_info(self):
		with sonarDBManager(self.db, readonly=True) as dbm:
			print("database path:             ", dbm.dbfile)
			print("database version:          ", dbm.get_db_version())
			print("database size:             ", dbm.get_db_size())
			print("unique sequences:          ", dbm.count_sequences())
			print("earliest genome import:    ", dbm.get_earliest_import())
			print("latest genome import:      ", dbm.get_latest_import())
			print("earliest sampling date:    ", dbm.get_earliest_date())
			print("latest sampling date:      ", dbm.get_latest_date())
		
	
	# output
	def rows_to_csv(self, rows, file=None, na="*** no data ***", tsv=False):
		if not rows:
			print(na, file=sys.stderr)
		else:
			file = sys.stdout if file is None else open(file, "w")
			sep = "\t" if tsv else ","
			writer = csv.DictWriter(file, rows[0].keys(), delimiter=sep, lineterminator=os.linesep)
			writer.writeheader()
			writer.writerows(rows)

	# Call SonarToVCF

	def var2vcf(self, acc, props, output, cpu, betaV2):
		if betaV2:
			print('Warning: Under Construction, please remove --betaV2')
			return True 	# sonartoVCF_v2.export2VCF(self.db,output, cpu, props)
		else:
			return	sonartoVCF.export2VCF(self.db,output, cpu, props, acc)


	def open_file(self, fname, mode="r", compressed=False, encoding=None):
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
		except:
			sys.exit("input error: " + fname + " cannot be opened.")
