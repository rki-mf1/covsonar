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
	parser_add.add_argument('-f', '--file', metavar="FILE", help="fasta file(s) containing DNA sequences to add", type=str, nargs="+", default=[])
	parser_add.add_argument('-o', '--outdir', metavar="DIR", help="use this directory for intermediate files (if set, temporary files are not deleted after import)", type=str, default=None)
	#parser_add.add_argument('-d', '--dir', metavar="DIR", help="add all files with fasta/pickle extension from the given input directory", type=str, default=None)
	parser_add.add_argument('--cache', metavar="DIR", help="import data from a pre-processed sonar cache directory", type=str, default="None")
	parser_add.add_argument('--paranoid', help="activate checks on seguid sequence collisions", action="store_true")

	# create the parser for the "match" command
	parser_match = subparsers.add_parser('match', parents=[general_parser], help='get mutations profiles for given accessions.')
	parser_match.add_argument('--include', '-i', metavar="STR", help="match genomes sharing the given mutation profile", type=str, action='append', nargs="+", default=None)
	parser_match.add_argument('--exclude', '-e', metavar="STR", help="match genomes not containing the mutation profile", type=str, action='append', nargs="+", default=None)
	parser_match.add_argument('--lineage', metavar="STR", help="match genomes of the given pangolin lineage(s) only", type=str, nargs="+", default=None)
	parser_match.add_argument('--acc', metavar="STR", help="match specific genomes defined by acession(s) only", type=str, nargs="+", default=None)
	parser_match.add_argument('--exclusive', help="do not allow additional mutations", action="store_true")
	parser_match.add_argument('--count', help="count only matching genomes", action="store_true")
	parser_match.add_argument('--logic', help="decide between OR or AND  logic when matching profile mutations (default: AND logic)", choices=["OR", "or", "AND", "and"], default="AND")
	parser_match.add_argument('--ambig', help="include ambiguos sites when reporting profiles (no effect when --count is used)", action="store_true")

	# create the parser for the "match" command
	parser_match = subparsers.add_parser('restore', parents=[general_parser], help='restore sequence (alignment) for a given accession.')
	parser_match.add_argument('--acc', metavar="STR", help="match specific genomes defined by acession(s) only", type=str, nargs = "+", required=True)
	parser_match.add_argument('--align', help="show aligned to reference sequence (if used, only a single accession can be processed)", action='store_true')

	# create the parser for the "view" command
	parser_add = subparsers.add_parser('view', parents=[general_parser], help='view database content.')
	parser_add.add_argument('-v', metavar="FILE", help="complete data view (default: dna)", choices=['dna', 'prot', 'profiles'], default="profiles")

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
			record = SeqIO.read(fastas[i], "fasta")
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

	def match(self, include_profiles, exclude_profiles, accessions, lineages, exclusive, ambig):
		self.rows_to_csv(self.db.match(include_profiles, exclude_profiles, None, None, exclusive, ambig), na="*** no match ***")

	def rows_to_csv(self, rows, na="*** no data ***"):
		if len(rows) == 0:
			print(na, file=sys.stderr)
		else:
			writer = csv.DictWriter(sys.stdout, rows[0].keys(), lineterminator=os.linesep)
			writer.writeheader()
			writer.writerows(rows)


if __name__ == "__main__":
	args = parse_args()
	snr = sonar(args.db)

	#add sequences
	if args.tool == "add":
		files = []
		if args.file:
			files += args.file

		if not files and args.cache is None:
			sys.exit("nothing to add.")

		if files:
			snr.add_fasta(*files, cachedir=args.outdir, cpus=args.cpus)
		else:
			snr.add_cache(cachedir=args.cache)

	#show
	if args.tool == "match":
		snr.match(args.include, args.exclude, args.acc, args.lineage, args.exclusive, args.ambig)

	#restore alignment
	if args.tool == "restore":
		for acc in args.accs:
			snr.restore(acc, args.align)

	#view data
	if args.tool == "view":
		rows = [x for x in snr.dbobj.iter_rows(args.v + "_view")]
		snr.rows_to_csv(rows)
