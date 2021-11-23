#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

VERSION = "1.1.2"
import os
import sys
import csv
import argparse
import gzip
import lzma
from lib import sonardb, sonartogo
from Bio import SeqIO
from tempfile import mkstemp
from collections import defaultdict
import re
from tqdm import tqdm
from multiprocessing import Pool

def parse_args():
	parser = argparse.ArgumentParser(prog="sonar.py", description="")
	subparsers = parser.add_subparsers(help='detect, store, and screen for mutations in SARS-CoV-2 genomic sequences')
	subparsers.dest = 'tool'
	subparsers.required = True

	#parent parser: db input
	general_parser = argparse.ArgumentParser(add_help=False)
	general_parser.add_argument('--db', metavar="DB_DIR", help="sonar database directory", type=str, required=True)
	general_parser.add_argument('--cpus', metavar="int", help="number of cpus to use (default: 1)", type=int, default=1)

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

	# create the parser for the "match" command
	parser_match = subparsers.add_parser('match', parents=[general_parser], help='get mutations profiles for given accessions.')
	parser_match.add_argument('--include', '-i', metavar="STR", help="match genomes sharing the given mutation profile", type=str, action='append', nargs="+", default=[])
	parser_match.add_argument('--exclude', '-e', metavar="STR", help="match genomes not containing the mutation profile", type=str, action='append', nargs="+", default=[])
	parser_match.add_argument('--lineage', metavar="STR", help="match genomes of the given pangolin lineage(s) only", type=str, nargs="+", default=[])
	parser_match.add_argument('--acc', metavar="STR", help="match specific genomes defined by acession(s) only", type=str, nargs="+", default=[])
	parser_match.add_argument('--zip', metavar="INT", help="only match genomes of a given region(s) defined by zip code(s)", type=str,  nargs="+", default=[])
	parser_match.add_argument('--date', help="only match genomes sampled at a certain sampling date or time frame. Accepts single dates (YYYY-MM-DD) or time spans (YYYY-MM-DD:YYYY-MM-DD).", nargs="+", type=str, default=[])
	parser_match.add_argument('--lab', metavar="STR", help="match genomes of the given lab only", type=str, nargs="+", default=[])
	parser_match.add_argument('--source', metavar="STR", help="match genomes of the given data source only", type=str, nargs="+", default=[])
	parser_match.add_argument('--collection', metavar="STR", help="match genomes of the given data collection only", type=str, nargs="+", default=[])
	parser_match.add_argument('--technology', metavar="STR", help="match genomes of the given sequencing technology only", type=str, nargs="+", default=[])
	parser_match.add_argument('--platform', metavar="STR", help="match genomes of the given sequencing platform only", type=str, nargs="+", default=[])
	parser_match.add_argument('--chemistry', metavar="STR", help="match genomes of the given sequencing chemistry only", type=str, nargs="+", default=[])
	parser_match.add_argument('--software', metavar="STR", help="software used for genome reconstruction", type=str, default=None)
	parser_match.add_argument('--version', metavar="STR", help="software version used for genome reconstruction", type=str, default=None)
	parser_match.add_argument('--material', metavar="STR", help="match genomes of the given sequencing chemistry only", type=str, nargs="+", default=[])
	parser_match.add_argument('--min_ct', metavar="STR", help="minimal ct value of samples resulting genomes are matched to", type=float, default=None)
	parser_match.add_argument('--max_ct', metavar="STR", help="maximal ct value of samples resulting genomes are matched to", type=float, default=None)
	parser_match_g1 = parser_match.add_mutually_exclusive_group()
	parser_match_g1.add_argument('--count', help="count instead of listing matching genomes", action="store_true")
	parser_match_g1.add_argument('--ambig', help="include ambiguos sites when reporting profiles (no effect when --count is used)", action="store_true")
	parser_match_g2 = parser_match.add_mutually_exclusive_group()
	parser_match_g2.add_argument('--only_frameshifts', help="show only genomes containing one or more frameshift mutations", action="store_true")
	parser_match_g2.add_argument('--no_frameshifts', help="show only genomes containing no frameshift mutation", action="store_true")
	parser_match_g2.add_argument('--tsv', help="use tsv instead of csv output", action="store_true")
	parser_match.add_argument('--debug', help="show database query for debugging", action="store_true")

	#create the parser for the "restore" command
	parser_restore = subparsers.add_parser('restore', parents=[general_parser], help='restore sequence(s) from the database.')
	parser_restore.add_argument('--acc', metavar="STR", help="acession(s) whose sequences are to be restored", type=str, default=[], nargs = "+")
	parser_restore.add_argument('--file', '-f', metavar="STR", help="file containing acession(s) whose sequences are to be restored (one accession per line)", type=str, default=None)

	#create the parser for the "Var2Vcf" command
	parser_var2vcf = subparsers.add_parser('var2vcf', parents=[general_parser], help='export variants from the database to vcf format.')
	parser_var2vcf.add_argument('--acc', metavar="STR", help="acession(s) whose sequences are to be exported", type=str, default=[], nargs = "+")
	parser_var2vcf.add_argument('--file', '-f', metavar="STR", help="file containing acession(s) whose sequences are to be exported (one accession per line)", type=str, default=None)
	parser_var2vcf.add_argument('--date', help="only match genomes sampled at a certain sampling date or time frame. Accepts single dates (YYYY-MM-DD) or time spans (YYYY-MM-DD:YYYY-MM-DD).", nargs="+", type=str, default=[])
	parser_var2vcf.add_argument('--output', '-o', metavar="STR", help="output file (merged vcf)", type=str, default=None, required=True)
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

	#create the parser for the "optimize" command
	parser_opt = subparsers.add_parser('optimize', parents=[general_parser], help='optimizes the database.')

	# version
	parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)

	return parser.parse_args()

class sonar():
	def __init__(self, db, gff=None, debug=False):
		self.dbfile = db if db else mkstemp()[1]
		self.db = sonardb.sonarDB(self.dbfile)
		self.gff = gff
		self.debug = debug

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


	def add(self, fnames, cachedir=None, cpus=1, timeout=600, force=False, paranoid=True, quiet=False, noprogress=False, compressed=False, source=None, collection=None, lab=None):
		'''
		Adds genome sequence(s) from given FASTA file(s) to the database.
		If cachedir is not defined, a temporary directory will be used as cache.
		'''

		# set display options
		disable_progressbar = False if not quiet and not noprogress else True
		print_steps = True if not quiet and noprogress else False

		# set global update vals
		updates = {}
		if source:
			updates['source'] = source
		if lab:
			updates['lab'] = lab
		if collection:
			updates['collection'] = collection

		# create db if necessary
		step = 0
		if cachedir and os.path.isdir(cachedir):
			step += 1
			if print_steps:
				print("[step", str(step) + "] restoring ... ")

		with sonardb.sonarCache(cachedir) as cache, sonardb.sonarDBManager(self.dbfile) as dbm:
			# db status
			if not quiet:
				dbstatus ={
						'genomes': dbm.count_genomes(),
						'seqs': dbm.count_sequences(),
						'labs': dbm.count_labs(),
				}

			# add fasta files to cache
			step += 1
			msg = "[step " + str(step) + "] caching ...   "
			if print_steps:
				print(msg)

			to_process = []
			to_import = defaultdict(set)
			to_update = set()

			for i in tqdm(range(len(fnames)), desc = msg, disable = disable_progressbar):
				with self.open_file(fnames[i], compressed=compressed) as handle:
					for record in SeqIO.parse(handle, "fasta"):
						acc = record.id
						descr = record.description
						seq = self.db.harmonize(record.seq)

						if len(seq) == 0:
							continue
						seqhash = self.db.hash(seq)
						genome_data = dbm.get_genomes(acc)

						if updates:
							to_update.add(acc)

						if genome_data:
							if genome_data['seqhash'] != seqhash:
								if not force:
									sys.exit("database error: " + acc + " exists in the database with a different sequence (use --force to allow updating)")
								dbm.delete_genome(acc)
							elif genome_data['description'] != descr:
								if not force:
									sys.exit("database error: " + acc + " exists in the database with a different description (use --force to allow updating)")
								dbm.update_genome(acc, description = descr)
								continue
							else:
								continue

						if dbm.seq_exists(seqhash):
							cache.prep_cached_files(seqhash)
							cache.write_info(seqhash)
							cache.add_seq(seqhash, seq)
						elif seqhash not in to_import:
							algn = cache.get_algn_fname(seqhash)
							fasta = cache.get_fasta_fname(seqhash)
							info = cache.get_info_fname(seqhash)

							if not os.path.isfile(fasta):
								unvalid_letters =  sorted(self.db.check_iupac_nt_code(seq))
								if unvalid_letters:
									sys.exit("input error: " + acc + " contains non-IUPAC characters (found: " + ", ".join(unvalid_letters) + ")")
								cache.add_seq(seqhash, seq)
								to_process.append([fasta, algn, info, seqhash, timeout])
							elif SeqIO.read(fasta, "fasta").seq != seq:
									sys.exit("cache error: sequence hash " + seqhash + " exists in cache but refers to a different sequence")
							elif not os.path.isfile(info):
								to_process.append([fasta, algn, info, seqhash, timeout])

						to_import[seqhash].add((acc, descr))

			step += 1
			msg = "[step " + str(step) + "] processing ..."
			if print_steps:
				print(msg)
			pool = Pool(processes=cpus)
			failed = set()
			for status, seqhash in tqdm(pool.imap_unordered(self.db.multi_process_fasta_wrapper, to_process), total=len(to_process), desc = msg, disable = disable_progressbar):
				if not status:
					failed.update([x[1] for x in cache.cache[seqhash]])
				if failed:
					print("timeout warning: following genomes were not added to the database since the respective sequence produced an timeout while aligning:", file=sys.stderr)
					for f in failed:
						print(f, file=sys.stderr)

			step += 1
			msg = "[step " + str(step) + "] importing ... "
			if print_steps:
				print(msg)
			self.db.import_genome_from_cache(cache.dirname, to_import, msg=msg, dbm=dbm, disable_progressbar=disable_progressbar)

			if updates:
				step += 1
				msg = "[step " + str(step) + "] updating ...  "
				if print_steps:
					print(msg)
				for i in tqdm(range(len(to_update)), desc = msg, disable = disable_progressbar):
					dbm.update_genome(to_update.pop(), **updates)

			# db status
			if not quiet:
				new_dbstatus ={
						'genomes': dbm.count_genomes(),
						'seqs': dbm.count_sequences(),
						'labs': dbm.count_labs(),
				}

				print("number of genomes:")
				print("\twas:   " + str(dbstatus['genomes']))
				print("\tnow:   " + str(new_dbstatus['genomes']))
				print("\tadded: " + str(new_dbstatus['genomes']-dbstatus['genomes']))
				print("number of unique sequences:")
				print("\twas:   " + str(dbstatus['seqs']))
				print("\tnow:   " + str(new_dbstatus['seqs']))
				print("\tadded: " + str(new_dbstatus['seqs']-dbstatus['seqs']))

	def match_genomes(self, include_profiles, exclude_profiles, accessions, lineages, zips, dates, labs, sources, collections, technologies, platforms, chemistries, software, software_version, materials, min_ct, max_ct, ambig, count=False, frameshifts=0, tsv=False):
		rows = self.db.match(include_profiles=include_profiles, exclude_profiles=exclude_profiles, accessions=accessions, lineages=lineages, zips=zips, dates=dates, labs=labs, sources=sources, collections=collections, technologies=technologies, chemistries=chemistries, software=software, software_version=software_version, materials=materials, min_ct=min_ct, max_ct=max_ct, ambig=ambig, count=count, frameshifts=frameshifts, debug=debug)
		if count:
			print(rows)
		else:
			self.rows_to_csv(rows, na="*** no match ***", tsv=tsv)

	def update_metadata(self, fname, accCol=None, lineageCol=None, zipCol=None, dateCol=None, gisaidCol=None, enaCol=None, labCol=None, sourceCol=None, collectionCol=None, technologyCol=None, platformCol=None, chemistryCol=None, softwareCol = None, versionCol = None, materialCol=None, ctCol=None, sep=",", pangolin=False, compressed=False):
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
		return self.db.restore_genome_using_dnavars(acc)

	def var2vcf(self, acc, date, output):

		sonartogo.export2VCF(self.dbfile,acc, date, output,self.db.refdescr)
		return 

	def view(self, acc):
		with sonardb.sonarDBManager(self.dbfile, readonly=True) as dbm:
			self.rows_to_csv(self.db.get_dna_vars(acc, dbm=dbm))

	def show_system_info(self):
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
	if hasattr(args, 'debug') and args.debug:
		debug = True
	else:
		debug = False

	if not args.db is None and args.tool != "add" and not os.path.isfile(args.db):
		sys.exit("input error: database does not exist.")

	snr = sonar(args.db, debug=debug)

	if not args.db is None:
		with sonardb.sonarDBManager(args.db, readonly=True) as dbm:
			dbm.check_db_compatibility()

	# add
	if args.tool == "add":

		# sanity check
		if not args.file and not args.cache:
			sys.exit("nothing to add.")

		snr.add(args.file, cachedir=args.cache, cpus=args.cpus, force=args.force, timeout=args.timeout, quiet=args.quiet, noprogress=args.noprogress, compressed=args.compressed, source=args.source, collection=args.collection, lab=args.lab)

	# match
	if args.tool == "match":
		# sanity check
		if args.date:
			regex = re.compile("^[0-9]{4}-[0-9]{2}-[0-9]{2}(?::[0-9]{4}-[0-9]{2}-[0-9]{2})?$")
			for d in args.date:
				if not regex.match(d):
					sys.exit("input error: " + d + " is not a valid date (YYYY-MM-DD) or time span (YYYY-MM-DD:YYYY-MM-DD).")
		if args.no_frameshifts:
			frameshifts = -1
		elif args.only_frameshifts:
			frameshifts = 1
		else:
			frameshifts = 0

		if args.lineage:
			args.lineage = [x.upper() for x in args.lineage]
		if args.lab:
			args.lab = [x.upper() for x in args.lab]
		if args.source:
			args.source = [x.upper() for x in args.source]
		if args.collection:
			args.collection = [x.upper() for x in args.collection]
		if args.technology:
			args.technology = [x.upper() for x in args.technology]
		if args.platform:
			args.platform = [x.upper() for x in args.platform]
		if args.chemistry:
			args.chemistry = [x.upper() for x in args.chemistry]
		if args.software:
			args.software = args.software.upper()
		if args.version:
			args.version = args.version.upper()
		if args.material:
			args.material = [x.upper() for x in args.material]

		snr.match_genomes(include_profiles=args.include, exclude_profiles=args.exclude, accessions=args.acc, lineages=args.lineage, zips=args.zip, dates=args.date, labs=args.lab, sources=args.source, collections=args.collection, technologies=args.technology, platforms=args.platform, chemistries=args.chemistry, software=args.software, software_version=args.version, materials=args.material, min_ct=args.min_ct, max_ct=args.max_ct, ambig=args.ambig, count=args.count, frameshifts=frameshifts, tsv=args.tsv)

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

	# Var2Vcf (export variants to  VCF)
	if args.tool == "var2vcf":
		if args.date:
			regex = re.compile("^[0-9]{4}-[0-9]{2}-[0-9]{2}(?::[0-9]{4}-[0-9]{2}-[0-9]{2})?$")
			for d in args.date:
				if not regex.match(d):
					sys.exit("input error: " + d + " is not a valid date (YYYY-MM-DD) or time span (YYYY-MM-DD:YYYY-MM-DD).")

		args.acc = set([x.strip() for x in args.acc])
		if args.file:
			if not os.path.isfile(args.file):
				sys.exit("input error: file " + args.file + " does not exist.")
			with snr.open_file(args.file, compressed='auto') as handle:
				for line in handle:
					args.acc.add(line.strip())
		if len(args.acc) == 0:
			print("Warning: export all variants, it will take a few minutes")

		snr.var2vcf(args.acc, args.date, args.output)

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
