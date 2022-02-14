#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

VERSION = "2.0.0"
import os
import sys
import csv
import argparse
import gzip
import lzma
import sonar
from sonar.action import sonarActions
from sonar.dbm import sonarDBManager
from Bio import SeqIO
from tempfile import mkstemp
from collections import defaultdict
import re
from tqdm import tqdm
from multiprocessing import Pool

class arg_namespace(object):
    pass

def parse_args():
	user_namespace = arg_namespace()
	parser = argparse.ArgumentParser(prog="sonar.py", description="")
	subparsers = parser.add_subparsers(help='detect, store, and screen for mutations in SARS-CoV-2 genomic sequences')
	subparsers.dest = 'tool'
	subparsers.required = True

	#parent parser: db input
	general_parser = argparse.ArgumentParser(add_help=False)
	general_parser.add_argument('--db', metavar="DB_DIR", help="sonar database directory", type=str, required=True)
	general_parser.add_argument('--cpus', metavar="int", help="number of cpus to use (default: 1)", type=int, default=1)
	general_parser.add_argument('--debug', help="activate debugging mode showing all sqllite queries on screen", action="store_true")

	# create the parser for the "setup" command
	parser_add = subparsers.add_parser('setup', parents=[general_parser], help='setup a new database.')

	# create the parser for the "add" command
	parser_add = subparsers.add_parser('add', parents=[general_parser], help='add genome sequences to the database.')
	parser_add_input = parser_add.add_mutually_exclusive_group()
	parser_add_input.add_argument('-f', '--file', metavar="FILE", help="fasta file(s) containing DNA sequences to add", type=str, nargs="+", default=[])
	parser_add_input.add_argument('-d', '--dir', metavar="DIR", help="add all fasta files (ending with \".fasta\" or \".fna\") from a given directory or directories", type=str, nargs="+", default=None)
	parser_add.add_argument('-c', '--cache', metavar="DIR", help="use (and restore data from) a given cache (if not set, a temporary cache is used and deleted after import)", type=str, default=None)
	parser_add.add_argument('-t', '--timeout', metavar="INT", help="timout for aligning sequences in seconds (default: 600)", type=int, default=600)
	parser_add.add_argument('--compressed', help="compression of input file format ('none', 'gz', 'xz', default: 'auto')", choices=['none', 'gz', 'xz', 'auto'], default='auto')
	parser_add.add_argument('--force', help="force updating of accessions if description or sequence has changed", action="store_true")
	parser_add.add_argument('--noprogress', '-p', help="do not show any progress bar", action="store_true")
	parser_add.add_argument('--source', help="define a common data source for all genomes", type=str, default=None)
	parser_add.add_argument('--collection', help="define a common data collection for all genomes", type=str, default=None)
	parser_add.add_argument('--lab', help="define a common lab for all genomes", type=str, default=None)
	parser_add.add_argument('--quiet', '-q', help="do not show any output", action="store_true")

	# create the parser for the "addprop" command
	parser_newprop = subparsers.add_parser('addprop', parents=[general_parser], help='add a new sample property to the database.')
	parser_newprop.add_argument('--name', metavar="STR", help="unique name of the new property (only alphanumeric characters and _ as only special character)", type=str, required=True)
	parser_newprop.add_argument('--descr', metavar="STR", help="description of the new property", type=str, required=True)
	parser_newprop.add_argument('--dtype', metavar="STR", help="data type of the new property", type=str, choices=['integer', 'float', 'text', 'date'], required=True)
	parser_newprop.add_argument('--qtype', metavar="STR", help="query type of the new property", type=str, choices=['numeric', 'float', 'date', 'text', 'zip'], required=True)
	parser_newprop.add_argument('--default', metavar="VAR", help="default value of the new property (none by default)", type=str, default=None)

	# create the parser for the "delprop" command
	parser_delprop = subparsers.add_parser('delprop', parents=[general_parser], help='delete property within the database.')
	parser_delprop.add_argument('--name', '-n', metavar="STR", help="name of the property to delete", type=str, required=True)

	# create the parser for the "match" command
	parser_match = subparsers.add_parser('match', parents=[general_parser], help='get mutations profiles for given accessions.')
	parser_match.add_argument('--profile', '-p', metavar="STR", help="match genomes sharing the given mutation profile", type=str, action='append', nargs="+", default=[])
	#parser_match.add_argument('--exclude', '-e', metavar="STR", help="match genomes not containing the mutation profile", type=str, action='append', nargs="+", default=[])
	##parser_match.add_argument('--lineage', metavar="STR", help="match genomes of the given pangolin lineage(s) only", type=str, nargs="+", default=[])
	#parser_match.add_argument('--acc', metavar="STR", help="match specific genomes defined by acession(s) only", type=str, nargs="+", default=[])
	#parser_match.add_argument('--zip', metavar="INT", help="only match genomes of a given region(s) defined by zip code(s)", type=str,  nargs="+", default=[])
	##parser_match.add_argument('--date', help="only match genomes sampled at a certain sampling date or time frame. Accepts single dates (YYYY-MM-DD) or time spans (YYYY-MM-DD:YYYY-MM-DD).", nargs="+", type=str, default=[])
	#parser_match.add_argument('--lab', metavar="STR", help="match genomes of the given lab only", type=str, nargs="+", default=[])
	#parser_match.add_argument('--source', metavar="STR", help="match genomes of the given data source only", type=str, nargs="+", default=[])
	#parser_match.add_argument('--collection', metavar="STR", help="match genomes of the given data collection only", type=str, nargs="+", default=[])
	#parser_match.add_argument('--technology', metavar="STR", help="match genomes of the given sequencing technology only", type=str, nargs="+", default=[])
	#parser_match.add_argument('--platform', metavar="STR", help="match genomes of the given sequencing platform only", type=str, nargs="+", default=[])
	#parser_match.add_argument('--chemistry', metavar="STR", help="match genomes of the given sequencing chemistry only", type=str, nargs="+", default=[])
	#parser_match.add_argument('--software', metavar="STR", help="software used for genome reconstruction", type=str, default=None)
	#parser_match.add_argument('--version', metavar="STR", help="software version used for genome reconstruction", type=str, default=None)
	#parser_match.add_argument('--material', metavar="STR", help="match genomes of the given sequencing chemistry only", type=str, nargs="+", default=[])
	#parser_match.add_argument('--min_ct', metavar="STR", help="minimal ct value of samples resulting genomes are matched to", type=float, default=None)
	#parser_match.add_argument('--max_ct', metavar="STR", help="maximal ct value of samples resulting genomes are matched to", type=float, default=None)
	#parser_match_g1 = parser_match.add_mutually_exclusive_group()
	parser_match.add_argument('--count', help="count instead of listing matching genomes", action="store_true")
	#parser_match_g1.add_argument('--ambig', help="include ambiguos sites when reporting profiles (no effect when --count is used)", action="store_true")
	#parser_match_g2 = parser_match.add_mutually_exclusive_group()
	#parser_match_g2.add_argument('--only_frameshifts', help="show only genomes containing one or more frameshift mutations", action="store_true")
	#parser_match_g2.add_argument('--no_frameshifts', help="show only genomes containing no frameshift mutation", action="store_true")
	#parser_match_g2.add_argument('--tsv', help="use tsv instead of csv output", action="store_true")

	#create the parser for the "restore" command
	parser_restore = subparsers.add_parser('restore', parents=[general_parser], help='restore sequence(s) from the database.')
	parser_restore.add_argument('--acc', metavar="STR", help="acession(s) whose sequences are to be restored", type=str, default=[], nargs = "+")
	parser_restore.add_argument('--file', '-f', metavar="STR", help="file containing acession(s) whose sequences are to be restored (one accession per line)", type=str, default=None)

	# create the parser for the "update" command
	parser_update = subparsers.add_parser('update', parents=[general_parser], help='add or update meta information.')
	parser_update_input = parser_update.add_mutually_exclusive_group()
	parser_update_input.add_argument('--pangolin', metavar="FILE", help="import linegae information from csv file created by pangolin", type=str, default=None)
	parser_update_input.add_argument('--csv', metavar="FILE", help="import metadata from a csv file", type=str, default=None)
	parser_update_input.add_argument('--tsv', metavar="FILE", help="import metadata from a tsv file", type=str, default=None)
	parser_update.add_argument('--fields', metavar="STR", help="if --csv or --tsv is used, define relevant columns like \"pango={colname_in_cs} zip={colname_in_cs} date={colname_in_csv}\"", type=str, nargs="+", default=None)
	parser_update.add_argument('--compressed', help="compression of input file format ('none', 'gz', 'xz', default: 'auto')", choices=['none', 'gz', 'xz', 'auto'], default='auto')

	# create the parser for the "info" command
	parser_info= subparsers.add_parser('info', help='show info')
	parser_info.add_argument('--db', metavar="DB_DIR", help="sonar database directory (optional)", type=str, default=None)

	# create the parser for the "props" command
	parser_props= subparsers.add_parser('props', parents=[general_parser], help='show database specific sample properties')

	#create the parser for the "optimize" command
	parser_opt = subparsers.add_parser('optimize', parents=[general_parser], help='optimizes the database.')

	#create the parser for the "dev" command
	parser_dev = subparsers.add_parser('dev', parents=[general_parser])

	# version
	parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
	args = parser.parse_known_args(namespace=user_namespace)
	if user_namespace.tool == "match" and hasattr(args[0], 'db'):
		with sonarDBManager(args[0].db, readonly = True, debug=args[0].debug) as dbm:
			for prop in dbm.properties.values():
				if prop['datatype'] == "integer":
					t = int
				elif prop['datatype'] == "float":
					t = float
				else:
					t = str
				parser_match.add_argument('--' + prop['name'], type=t, nargs = '+', default=argparse.SUPPRESS)

	#return parser.parse_args()
	return parser.parse_args(namespace=user_namespace)

class sonar():
	def __init__(self, db, gff=None, debug=False):
		self.dbfile = db if db else mkstemp()[1]
		self.db = sonardb.sonarDB(self.dbfile)
		self.gff = gff
		self.debug = debug

	def update_metadata(self, fname, accCol=None, lineageCol=None, zipCol=None, dateCol=None, gisaidCol=None, enaCol=None, labCol=None, sourceCol=None, collectionCol=None, technologyCol=None, platformCol=None, chemistryCol=None, softwareCol = None, versionCol = None, materialCol=None, ctCol=None, sep=",", pangolin=False, compressed=False):
		'''
		This function takes a file name and a list of column names and returns a dictionary of dictionaries.
		The outer dictionary is keyed by the accession number and the inner dictionary is keyed by the
		column name
		
		:param fname: the name of the file containing the metadata
		:param accCol: the column in the metadata file that contains the accession number
		:param lineageCol: the column in the metadata file that contains the lineage
		:param zipCol: the column in the metadata file that contains the zip code
		:param dateCol: the column in the metadata file that contains the date of sampling
		:param gisaidCol: the column in the metadata file that contains the GISAID accession number
		:param enaCol: the column in the metadata file that contains the ENA accession
		:param labCol: the column in the metadata file that contains the lab of origin
		:param sourceCol: the column in the metadata file that contains the source of the sequence (e.g.
		"human", "environmental", "animal")
		:param collectionCol: the column in the metadata file that contains the collection name
		:param technologyCol: the name of the column in the metadata file that contains the sequencing
		technology used to generate the reads
		:param platformCol: the column in the metadata file that contains the platform
		:param chemistryCol: the column in the metadata file that contains the chemistry used
		:param softwareCol: the name of the column containing the software used to generate the fasta file
		:param versionCol: the column in the metadata file that contains the version number
		:param materialCol: the column in the metadata file that contains the material type (e.g. cell,
		tissue, etc.)
		:param ctCol: the column in the metadata file that contains the ct value
		:param sep: the delimiter used in the metadata file, defaults to , (optional)
		:param pangolin: if True, the file is a pangolin metadata file, defaults to False (optional)
		:param compressed: If True, the file is compressed with gzip, defaults to False (optional)
		'''
		updates = defaultdict(dict)
		if pangolin:
			with self.open_file(fname, compressed=compressed, encoding='utf-8-sig') as handle:
				lines = csv.DictReader(handle, delimiter = ',', quoting=csv.QUOTE_MINIMAL)
				for line in lines:
					acc = line['Sequence name'].split(" ")[0]
					updates[acc]['lineage'] = line['Lineage']
		elif accCol:
			with self.open_file(fname, compressed=compressed) as handle:
				lines = csv.DictReader(handle, delimiter = sep)
				for line in lines:
					acc = line[accCol]
					if lineageCol and (acc not in updates or 'lineage' not in updates[acc]):
						updates[acc]['lineage'] = line[lineageCol].upper()
					if zipCol and line[zipCol]:
						updates[acc]['zip'] = line[zipCol]
					if dateCol and line[dateCol]:
						updates[acc]['date'] = line[dateCol]
					if gisaidCol and line[gisaidCol]:
						updates[acc]['gisaid'] = line[gisaidCol]
					if enaCol and line[enaCol]:
						updates[acc]['ena'] = line[enaCol]
					if collectionCol and line[collectionCol]:
						updates[acc]['collection'] = line[collectionCol].upper()
					if sourceCol and line[sourceCol]:
						updates[acc]['source'] = line[sourceCol].upper()
					if labCol and line[labCol]:
						updates[acc]['lab'] = line[labCol].upper()
					if technologyCol and line[technologyCol]:
						updates[acc]['technology'] = line[technologyCol].upper()
					if chemistryCol and line[chemistryCol]:
						updates[acc]['chemistry'] = line[chemistryCol].upper()
					if platformCol and line[platformCol]:
						updates[acc]['platform'] = line[platformCol].upper()
					if softwareCol and line[softwareCol]:
						updates[acc]['software'] = line[softwareCol].upper()
					if versionCol and line[versionCol]:
						updates[acc]['version'] = line[versionCol].upper()
					if materialCol:
						updates[acc]['material'] = line[materialCol].upper()
					if ctCol and line[ctCol]:
						try:
							updates[acc]['ct'] = float(line[ctCol])
						except:
							sys.exit("metadata error: " + line[ctCol] + " is not a valid ct value (accession: " + acc + ")")
		with sonardb.sonarDBManager(self.dbfile) as dbm:
			for acc, update in updates.items():
				dbm.update_genome(acc, **update)

	def restore(self, acc):
		'''
		It takes a genome accession number and returns a genome object.
		
		:param acc: the accession number of the genome to be restored
		:return: A list of dictionaries. Each dictionary is a genome.
		'''
		return self.db.restore_genome_using_dnavars(acc)

	def view(self, acc):
		'''
		Self.rows_to_csv(self.db.get_dna_vars(acc, dbm=dbm))
		
		This function is the heart of the program. It is called by the view() function.
		It takes a single argument, acc, which is the accession number of the sequence
		of interest. It returns a list of rows, which is a list of lists. Each list
		within the list is a row of data
		
		:param acc: the accession number of the sequence to be viewed
		'''
		with sonardb.sonarDBManager(self.dbfile, readonly=True) as dbm:
			self.rows_to_csv(self.db.get_dna_vars(acc, dbm=dbm))

	def show_system_info(self):
		'''
		Prints out the version of the database, the name of the reference genome, the length of the
		reference genome, the names of the annotated proteins, and the translation table used
		'''
		print("sonarDB version:       ", self.db.get_version())
		print("reference genome:      ", self.db.refdescr)
		print("reference length:      ", str(len(self.db.refseq)) + "bp")
		print("annotated proteins:    ", ", ".join(self.db.refgffObj.symbols))
		print("used translation table:", self.db.translation_table)

	def show_db_info(self):
		with sonardb.sonarDBManager(self.dbfile, readonly=True) as dbm:
			print("database path:             ", dbm.dbfile)
			print("database version:          ", dbm.get_db_version())
			print("database size:             ", self.get_db_size())
			g = dbm.count_genomes()
			print("genomes:                   ", g)
			print("unique sequences:          ", dbm.count_sequences())
			print("labs:                      ", dbm.count_labs())
			print("earliest genome import:    ", dbm.get_earliest_import())
			print("latest genome import:      ", dbm.get_latest_import())
			print("earliest sampling date:    ", dbm.get_earliest_date())
			print("latest sampling date:      ", dbm.get_latest_date())
			print("metadata:          ")
			fields = sorted(['lab', 'source', 'collection', 'technology', 'platform', 'chemistry', 'software', 'software_version', 'material', 'ct', 'gisaid', 'ena', 'lineage', 'zip', 'date'])
			maxlen = max([len(x) for x in fields])
			for field in fields:
				if g == 0:
					c = 0
					p = 0
				else:
					c = dbm.count_metadata(field)
					p = c/g*100
				spacer = " " * (maxlen-len(field))
				print("   " + field + " information:" + spacer, f"{c} ({p:.{2}f}%)")

	def rows_to_csv(self, rows, file=None, na="*** no data ***", tsv=False):
		'''
		The rows_to_csv function takes a list of dictionaries and writes them to a file in CSV format
		
		:param rows: the list of dictionaries to write to the file
		:param file: The file to write to. If None, the output is written to sys.stdout
		:param na: the value to print if there are no rows, defaults to *** no data *** (optional)
		:param tsv: If True, the output will be written in tab-separated format, defaults to False
		(optional)
		'''
		if len(rows) == 0:
			print(na, file=sys.stderr)
		else:
			file = sys.stdout if file is None else open(file, "w")
			sep = "\t" if tsv else ","
			writer = csv.DictWriter(file, rows[0].keys(), delimiter=sep, lineterminator=os.linesep)
			writer.writeheader()
			writer.writerows(rows)

	def get_db_size(self, decimal_places=3):
		size = os.path.getsize(self.dbfile)
		for unit in ['B','KiB','MiB','GiB','TiB']:
			if size < 1024.0:
				break
			size /= 1024.0
		return f"{size:.{decimal_places}f}{unit}"

def process_update_expressions(expr):
	allowed = {"accession": "accCol", "lineage": "lineageCol", "date": "dateCol", "zip": "zipCol", "gisaid": "gisaidCol", "ena": "enaCol", "collection": "collectionCol", "technology": "technologyCol", "platform": "platformCol", "chemistry": "chemistryCol", "software": "softwareCol", "version": "versionCol", "material": "materialCol", "ct": "ctCol", "source": "sourceCol", "lab": "labCol"}
	fields = {}
	for val in expr:
		val = val.split("=")
		if val[0] not in allowed or len(val) == 1:
			sys.exit("input error: " + val[0] + " is not a valid expression")
		key = allowed[val[0]]
		if key in fields:
			sys.exit("input error: multiple assignments for " + val[0])
		fields[key] = "=".join(val[1:])
		if 'accCol' not in fields:
			sys.exit("input error: an accession column has to be defined.")
	return fields

if __name__ == "__main__":
	args = parse_args()

	# setup
	if args.tool == "setup":
		sonarActions.setup_db(args.db, debug=args.debug)
		exit(0)

	snr = sonarActions(args.db, debug=args.debug)
	with sonarDBManager(args.db, readonly=True) as dbm:
		dbm.check_db_compatibility()

	# match
	if args.tool == "match":
		with sonarDBManager(args.db, readonly=True) as dbm:
			props = {}
			for pname in dbm.properties:
				if hasattr(args, pname):
					props[pname] = getattr(args, pname)
		snr.match(args.profile, props, count=args.count)

	# newprop
	if args.tool == "addprop":
		with sonarDBManager(args.db, debug=args.debug) as dbm:
			dbm.add_property(args.name, args.dtype, args.qtype, args.descr, args.default)

	# delprop
	if args.tool == "delprop":
		with sonarDBManager(args.db, debug=args.debug) as dbm:
			if args.name not in dbm.properties:
				sys.exit("input error: unknown property.")
			a = dbm.count_property(args.name)
			b = dbm.count_property(args.name, ignore_standard = True)
			print("WARNING: There are", a, "samples with content for this property. Amongst those,", b, "samples do not share the default value of this property.")
			decision = ""
			while decision not in ("YES", "no"):
				decision = input("Do you really want to delete this property? [YES/no]: ")
			if decision == "YES":
				dbm.delete_property(args.name)
				print("property deleted.")
			else:
				print("property not deleted.")

	# delprop
	if args.tool == "props":
		snr.show_props()
	# update
	if args.tool == "update":
		fields={}
		if args.csv:
			cols = process_update_expressions(args.fields)
			snr.update_metadata(args.csv, **cols, sep=",", pangolin=False, compressed=args.compressed)
		elif args.tsv:
			cols = process_update_expressions(args.fields)
			snr.update_metadata(args.tsv, **cols, sep="\t", pangolin=False, compressed=args.compressed)
		elif args.pangolin:
			snr.update_metadata(args.pangolin, pangolin=True, compressed=args.compressed)
		else:
			print("nothing to update.")

	# restore
	if args.tool == "restore":
		args.acc = set([x.strip() for x in args.acc])
		if args.file:
			if not os.path.isfile(args.file):
				sys.exit("input error: file " + args.file + " does not exist.")
			with snr.open_file(fname, compressed=args.file,) as handle:
				for line in handle:
					args.acc.add(line.strip())
		if len(args.acc) == 0:
			sys.exit("input error: nothing to restore.")
		for acc in filter(None, args.acc):
			print("\n".join(snr.restore(acc)))

	# view
	if args.tool == "view":
		snr.view(args.acc)

	# info
	if args.tool == "info":
		snr.show_system_info()
		if args.db:
			print()
			snr.show_db_info()

	# optimize
	if args.tool == "optimize":
		sonardb.sonarDBManager.optimize(args.db)

	# optimize
	if args.tool == "dev":
		print("***dev mode***")
		with sonarDBManager(args.db, debug=args.debug) as dbm:
			for feature in dbm.get_annotation():
				print(())
