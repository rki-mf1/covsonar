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
from joblib import Parallel, delayed
import itertools
import re
from tqdm import tqdm
import difflib
import glob
from time import sleep
import shutil
import traceback

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
	parser_add_input.add_argument('-o', '--outdir', metavar="DIR", help="use this directory for intermediate files (if set, temporary files are not deleted after import)", type=str, default=None)
	parser_add.add_argument('--cache', metavar="DIR", help="import data from a pre-processed sonar cache directory", type=str, default=None)
	parser_add.add_argument('--paranoid', help="activate checks on seguid sequence collisions", action="store_true")

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

	# create the parser for the "restore" command
	# parser_match = subparsers.add_parser('restore', parents=[general_parser], help='restore sequence (alignment) for a given accession.')
	# parser_match.add_argument('--acc', metavar="STR", help="match specific genomes defined by acession(s) only", type=str, nargs = "+", required=True)
	# parser_match.add_argument('--align', help="show aligned to reference sequence (if used, only a single accession can be processed)", action='store_true')

	# create the parser for the "update" command
	parser_update = subparsers.add_parser('update', parents=[general_parser], help='add or update meta information.')
	parser_update_input = parser_update.add_mutually_exclusive_group()
	parser_update_input.add_argument('--pangolin', metavar="FILE", help="import linegae information from csv file created by pangolin", type=str, default=None)
	parser_update_input.add_argument('--csv', metavar="FILE", help="import metadata from a csv file", type=str, default=None)
	parser_update_input.add_argument('--tsv', metavar="FILE", help="import metadata from a tsv file", type=str, default=None)
	parser_update.add_argument('--fields', metavar="STR", help="if --csv or --tsv is used, define relevant columns like \"pango={colname_in_cs} zip={colname_in_cs} date={colname_in_csv}\"", type=str, nargs="+", default=None)


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

	def add_fasta(self, *fnames, cachedir=None, cpus=1, paranoid=True):
		'''
		Adds genome sequence(s) from the given FASTA file(s) to the database.
		If dir is not defined, a temporary directory will be used as cache.
		'''
		cache = sonardb.sonarCache(cachedir)
		with sonardb.sonarDBManager(self.dbfile) as dbm:
			pass

		# add fasta files to cache
		msg = "[step 1 of 3] caching ...   "
		for i in tqdm(range(len(fnames)), desc = msg):
			cache.add_fasta(fnames[i])

		# execute adding method of the database module on the generated files
		seqhashes = cache.get_cached_seqhashes()
		msg = "[step 2 of 3] processing ..."
		r = Parallel(n_jobs=cpus)(delayed(self.db.process_fasta)(cache.cached_fasta_name(seqhashes[x]), cache.cached_pickle_name(seqhashes[x])) for x in tqdm(range(len(seqhashes)), desc = msg))

		msg = "[step 3 of 3] importing ... "
		data = []
		for seqhash in seqhashes:
			fname = cache.cached_pickle_name(seqhash)
			for acc, descr in cache.get_acc_descr(seqhash):
				data.append((fname, acc, descr))
		with sonardb.sonarDBManager(self.dbfile) as dbm:
			for i in tqdm(range(len(data)), desc = msg):
				dat = [data[i][1], data[i][2]] + cache.load_pickle(data[i][0])[2:]
				self.db.import_genome(*dat, paranoid = True, dbm=dbm)

	def add_cache(self, cachedir=None, paranoid=True):
		'''
		Adds genome sequence(s) from a preprocessed sonar cache.
		'''
		msg = "[step 1 of 2] restoring ... "
		fastas = [x for x in glob.glob(os.path.join(cachedir, "*"  + ".fasta"))]
		if not len(fastas):
			sys.exit("cache error: not a valid cache")
		data = []
		for i in tqdm(range(len(fastas)), desc = msg):
			record = SeqIO.readfieldList(fastas[i], "fasta")
			acc, descr, seq  = record.id, record.description, str(record.seq)
			seqhash = self.db.hash(seq)
			picklefile = os.path.join(cachedir, sonardb.sonarCache.slugify(seqhash) + ".pickle")
			if not os.path.isfile(picklefile):
				sys.exit("cache error: cache seems to be corrupt")
			data.append((picklefile, acc, descr))

		msg = "[step 2 of 2] importing ... "
		with sonardb.sonarDBManager(self.dbfile) as dbm:
			for i in tqdm(range(len(data)), desc = msg):
				dat = [data[i][1], data[i][2]] + sonardb.sonarCache.load_pickle(data[i][0])[2:]
				self.db.import_genome(*dat, paranoid = True, dbm=dbm)

	def match(self, include_profiles, exclude_profiles, accessions, lineages, zips, dates, exclusive, ambig, count=False):
		rows = self.db.match(include_profiles, exclude_profiles, accessions, lineages, zips, dates, exclusive, ambig)
		if count:
			print(len(rows))
		else:
			self.rows_to_csv(rows, na="*** no match ***")

	def update_metadata(self, fname, accCol=None, lineageCol=None, zipCol=None, dateCol=None, gisaidCol=None, enaCol=None, sep=",", pangolin=False):
		updates = defaultdict(dict)
		if pangolin:
			with open(fname, "r") as handle:
				lines = csv.DictReader(handle, delimiter = ',',  quoting=csv.QUOTE_MINIMAL)
				for line in lines:
					acc = line['﻿Sequence name']
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

	def rows_to_csv(self, rows, na="*** no data ***"):
		if len(rows) == 0:
			print(na, file=sys.stderr)
		else:
			writer = csv.DictWriter(sys.stdout, rows[0].keys(), lineterminator=os.linesep)
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

		# fasta file input
		if args.file:
			files = args.file

		#dir input
		elif args.dir:
			files = []
			for d in args.dir:
				if not os.path.isdir(dir):
					sys.exit("input error: " + dir + " is not a valid directory")
				files += [x for x in glob.glob(os.path.join(dir, "*.fasta *.fna"))]

		# cache input
		elif args.cache:
			if not os.path.isdir(args.cache):
				sys.exit("input error: " + args.cache + " does not exist")

		# sanity check
		if not files and not args.cache:
			sys.exit("nothing to add.")

		if args.cache:
			snr.add_cache(cachedir=args.cache)
		else:
			snr.add_fasta(*files, cachedir=args.outdir, cpus=args.cpus)


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


	# optimize
	if args.tool == "optimize":
		sonardb.sonarDBManager.optimize(args.db)
