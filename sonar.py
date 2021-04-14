#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

VERSION = "0.0.9"
import os
import sys
import csv
import argparse
from lib import sonardb
from Bio import SeqIO
import numpy as np
import tempfile
from collections import defaultdict
import itertools
import re
from tqdm import tqdm
import difflib
import glob
from time import sleep
import shutil
import traceback
from multiprocessing import Pool
import types

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
	#parser_add.add_argument('--paranoid', help="activate checks on seguid sequence collisions", action="store_true")

	# create the parser for the "match" command
	parser_match = subparsers.add_parser('match', parents=[general_parser], help='get mutations profiles for given accessions.')
	parser_match.add_argument('--include', '-i', metavar="STR", help="match genomes sharing the given mutation profile", type=str, action='append', nargs="+", default=None)
	parser_match.add_argument('--exclude', '-e', metavar="STR", help="match genomes not containing the mutation profile", type=str, action='append', nargs="+", default=None)
	parser_match.add_argument('--lineage', metavar="STR", help="match genomes of the given pangolin lineage(s) only", type=str, nargs="+", default=None)
	parser_match.add_argument('--acc', metavar="STR", help="match specific genomes defined by acession(s) only", type=str, nargs="+", default=None)
	parser_match.add_argument('--zip', metavar="INT", help="only match genomes of a given region(s) defined by zip code(s)", type=int,  nargs="+", default=None)
	parser_match.add_argument('--date', help="only match genomes sampled at a certain sampling date or time frame. Accepts single dates (YYYY-MM-DD) or time spans (YYYY-MM-DD:YYYY-MM-DD).", nargs="+", type=str)
	parser_match.add_argument('--exclusive', help="do not allow additional mutations", action="store_true")
	parser_match.add_argument('--count', help="count instead of listing matching genomes", action="store_true")
	parser_match.add_argument('--ambig', help="include ambiguos sites when reporting profiles (no effect when --count is used)", action="store_true")

	#create the parser for the "restore" command
	parser_restore = subparsers.add_parser('restore', parents=[general_parser], help='restore sequence(s) from the database.')
	parser_restore.add_argument('--acc', metavar="STR", help="acession(s) whose sequences are to be restored", type=str, nargs = "+", required=True)
	# parser_match.add_argument('--align', help="show aligned to reference sequence (if used, only a single accession can be processed)", action='store_true')

	#create the parser for the "view" command
	parser_view = subparsers.add_parser('view', parents=[general_parser], help='show dna profile.')
	parser_view.add_argument('--acc', metavar="STR", help="accession to consider", type=str, required=True)
	# parser_match.add_argument('--align', help="show aligned to reference sequence (if used, only a single accession can be processed)", action='store_true')

	# create the parser for the "update" command
	parser_update = subparsers.add_parser('update', parents=[general_parser], help='add or update meta information.')
	parser_update_input = parser_update.add_mutually_exclusive_group()
	parser_update_input.add_argument('--pangolin', metavar="FILE", help="import linegae information from csv file created by pangolin", type=str, default=None)
	parser_update_input.add_argument('--csv', metavar="FILE", help="import metadata from a csv file", type=str, default=None)
	parser_update_input.add_argument('--tsv', metavar="FILE", help="import metadata from a tsv file", type=str, default=None)
	parser_update.add_argument('--fields', metavar="STR", help="if --csv or --tsv is used, define relevant columns like \"pango={colname_in_cs} zip={colname_in_cs} date={colname_in_csv}\"", type=str, nargs="+", default=None)

	# create the parser for the "frameshift" command
	parser_frameshift = subparsers.add_parser('frameshift', parents=[general_parser], help='show all genomes with frameshifts')
	parser_frameshift_out = parser_frameshift.add_mutually_exclusive_group()
	parser_frameshift_out.add_argument('--count', help="count instead of listing matching genomes", action="store_true")
	parser_frameshift_out.add_argument('-o', '--out', help="write output to file(please note: the given file will be overwritten).", type=str, default=None)

	#create the parser for the "optimize" command
	parser_opt = subparsers.add_parser('optimize', parents=[general_parser], help='optimizes the database.')

	# version
	parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)

	return parser.parse_args()

class sonar():
	def __init__(self, db, gff=None):
		self.dbfile = db
		self.db = sonardb.sonarDB(self.dbfile)
		self.gff = gff

	def add(self, fnames=None, cachedir=None, cpus=1, timeout=600):
		'''
		Adds genome sequence(s) from the given FASTA file(s) to the database.
		If dir is not defined, a temporary directory will be used as cache.
		'''

		step = 0
		if cachedir and os.path.isdir(cachedir):
			step += 1
			print("[step", step, "] restoring ... ", file=sys.stderr)

		# create db if necessary
		if not os.path.isfile(self.dbfile):
			with sonardb.sonarDBManager(self.dbfile) as dbm:
				pass

		step = 0
		if cachedir and os.path.isdir(cachedir):
			step += 1
			print("[step", step, "] restoring ... ")

		with sonardb.sonarCache(cachedir) as cache:

			# add fasta files to cache
			if fnames:
				step += 1
				msg = "[step " + str(step) + "] caching ...   "
				for i in tqdm(range(len(fnames)), desc = msg):
					cache.add_fasta(fnames[i])
			step += 1
			msg = "[step " + str(step) + "] processing ..."
			new = []
			for seqhash in cache.cache:
				this = cache.get_cache_files(seqhash)
				if not os.path.isfile(this[1]):
					new.append(list(this) + [seqhash, timeout])

			pool = Pool(processes=cpus)
			failed = set()
			for status, seqhash in tqdm(pool.imap_unordered(self.db.multi_process_fasta_wrapper, new), total=len(new), desc = msg):
				if not status:
					failed.update([x[1] for x in cache.cache[seqhash]])
					del cache.cache[seqhash]
			if failed:
				print("timeout warning: following genomes were not added to the database since the respective sequence produced an timeout while aligning:", file=sys.stderr)
				for f in failed:
					print(f, file=sys.stderr)

			cache.write_idx(backup=True)

			step += 1
			msg = "[step " + str(step) + "] importing ... "
			data = []
			seqhashes = list(cache.cache.keys())
			self.db.import_genome_from_cache(cache.dirname, msg=msg)

	def match(self, include_profiles, exclude_profiles, accessions, lineages, zips, dates, exclusive, ambig, count=False, show_frame_shifts_only=False):
		rows = self.db.match(include_profiles, exclude_profiles, accessions, lineages, zips, dates, exclusive, ambig)
		if count:
			print(len(rows))
		else:
			self.rows_to_csv(rows, na="*** no match ***")

	def update_metadata(self, fname, accCol=None, lineageCol=None, zipCol=None, dateCol=None, gisaidCol=None, enaCol=None, sep=",", pangolin=False):
		updates = defaultdict(dict)
		if pangolin:
			with open(fname, "r", encoding='utf-8-sig') as handle:
				lines = csv.DictReader(handle, delimiter = ',', quoting=csv.QUOTE_MINIMAL)
				for line in lines:
					acc = line['Sequence name']
					updates[acc]['lineage'] = line['Lineage']
		elif accCol:
			with open(fname, "r") as handle:
				lines = csv.DictReader(handle, delimiter = sep)
				for line in lines:
					acc = line[accCol]
					if lineageCol and (acc not in updates or 'lineage' not in updates[acc]) :
						updates[acc]['lineage'] = line[lineageCol]
					if zipCol:
						updates[acc]['zip'] = line[zipCol]
					if dateCol:
						updates[acc]['date'] = line[dateCol]
					if gisaidCol:
						updates[acc]['gisaid'] = line[gisaidCol]
					if enaCol:
						updates[acc]['ena'] = line[enaCol]

		with sonardb.sonarDBManager(self.dbfile) as dbm:
			for acc, update in updates.items():
				elems = update.items()
				fieldList = [x[0] for x in elems]
				valList = [x[1] for x in elems]
				dbm.update("genome", fieldList, valList, "accession = ?", [acc])

	def restore(self, acc):
		with sonardb.sonarDBManager(self.dbfile, readonly=True) as dbm:
			return self.db.restore_genome(acc, dbm=dbm)

	def view(self, acc):
		with sonardb.sonarDBManager(self.dbfile, readonly=True) as dbm:
			self.rows_to_csv(self.db.get_dna_vars(acc, dbm=dbm))

	def show_frameshift(self, count=False, file=None):
		with sonardb.sonarDBManager(self.dbfile, readonly=True) as dbm:
			rows = [x for x in self.db.iter_frameshifts(dbm)]
		if count:
			print(len(rows))
		else:
			self.rows_to_csv(rows, file=file, na="*** no match ***")

	def rows_to_csv(self, rows, file=None, na="*** no data ***"):
		if len(rows) == 0:
			print(na, file=sys.stderr)
		else:
			file = sys.stdout if file is None else open(file, "w")
			writer = csv.DictWriter(file, rows[0].keys(), lineterminator=os.linesep)
			writer.writeheader()
			writer.writerows(rows)


def process_update_expressions(expr):
	allowed = {"accession": "accCol", "lineage": "lineageCol", "date": "dateCol", "zip": "zipCol", "gisaid": "gisaidCol", "ena": "enaCol"}
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
	snr = sonar(args.db)

	# add
	if args.tool == "add":

		# sanity check
		if not args.file and not args.cache:
			sys.exit("nothing to add.")

		snr.add(args.file, cachedir=args.cache, cpus=args.cpus, timeout=args.timeout)

	# match
	if args.tool == "match":
		# sanity check
		if args.date:
			regex = re.compile("^[0-9]{4}-[0-9]{2}-[0-9]{2}(?::[0-9]{4}-[0-9]{2}-[0-9]{2})?$")
			for d in args.date:
				if not regex.match(d):
					sys.exit("input error: " + d + " is not a valid date (YYYY-MM-DD) or time span (YYYY-MM-DD:YYYY-MM-DD).")
		snr.match(args.include, args.exclude, args.acc, args.lineage, args.zip, args.date, args.exclusive, args.ambig, args.count)

	# update
	if args.tool == "update":
		fields={}
		if args.csv:
			cols = process_update_expressions(args.fields)
			snr.update_metadata(args.csv, **cols, sep=",", pangolin=False)
		elif args.tsv:
			cols = process_update_expressions(args.fields)
			snr.update_metadata(args.tsv, **cols, sep="\t", pangolin=False)
		elif args.pangolin:
			snr.update_metadata(args.pangolin, pangolin=True)
		else:
			print("nothing to update.")

	# restore
	if args.tool == "restore":
		for acc in args.acc:
			print("\n".join(snr.restore(acc)))

	# view
	if args.tool == "view":
		snr.view(args.acc)

	# frameshift
	if args.tool == "frameshift":
		# sanity check
		snr.show_frameshift(args.count, args.out)

	# optimize
	if args.tool == "optimize":
		sonardb.sonarDBManager.optimize(args.db)
