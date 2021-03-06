#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

VERSION = "0.0.9"
import os
import sys
import argparse
from lib import sonardb
from Bio import SeqIO

def parse_args():
	parser = argparse.ArgumentParser(prog="sonar.py", description="")
	subparsers = parser.add_subparsers(help='detect, store, and screen for mutations in SARS-CoV-2 genomic sequences')
	subparsers.dest = 'tool'
	subparsers.required = True

	#parent parser: db input
	db_parser = argparse.ArgumentParser(add_help=False)
	db_parser.add_argument('--db', metavar="DB_DIR", help="sonar database directory", type=str, required=True)

	# create the parser for the "count" command
	parser_add = subparsers.add_parser('add', parents=[db_parser], help='add genome sequences to the database.')
	parser_add.add_argument('-f', '--fasta', metavar="FILE", help="fasta file(s) containing DNA sequences to add", type=str, nargs="+")
	parser_add.add_argument('--paranoid', help="activate checks on seguid sequence collisions", action="store_true")

	parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)

	return parser.parse_args()

class sonar():
	def __init__(self, db):
		self.db = sonardb.sonarDB(db)

	def add(self, *fnames, paranoid=False):
		ref = os.path.join(os.path.dirname(os.path.realpath(__file__)), "resources", "reference.fna")
		for fname in fnames:
			for record in SeqIO.parse(fname, "fasta"):
				self.db.add(record.id, record.description, str(record.seq), ref, paranoid)

	def iter_genomes(self, view="dna"):
		for row in self.db.iter_rows(view + "_view"):
			yield row

if __name__ == "__main__":
	args = parse_args()
	snr = sonar(args.db)

	#add sequences
	if args.tool == "add":
		snr.add(*args.fasta, paranoid=args.paranoid)
		#for row in snr.iter_genomes():
		#	print(row)
