#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

VERSION = "0.0.9"
import os
import sys
import argparse
from lib import sonardb
from Bio import SeqIO
import csv

def parse_args():
	parser = argparse.ArgumentParser(prog="sonar.py", description="")
	subparsers = parser.add_subparsers(help='detect, store, and screen for mutations in SARS-CoV-2 genomic sequences')
	subparsers.dest = 'tool'
	subparsers.required = True

	#parent parser: db input
	db_parser = argparse.ArgumentParser(add_help=False)
	db_parser.add_argument('--db', metavar="DB_DIR", help="sonar database directory", type=str, required=True)

	# create the parser for the "add" command
	parser_add = subparsers.add_parser('add', parents=[db_parser], help='add genome sequences to the database.')
	parser_add.add_argument('-f', '--fasta', metavar="FILE", help="fasta file(s) containing DNA sequences to add", type=str, nargs="+")
	parser_add.add_argument('--paranoid', help="activate checks on seguid sequence collisions", action="store_true")

	# create the parser for the "view" command
	parser_add = subparsers.add_parser('view', parents=[db_parser], help='view database content.')
	parser_add.add_argument('-b', '--branch', metavar="FILE", help="data branch (default: dna)", choices=['dna', 'prot'], default="dna")

	# version
	parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)

	return parser.parse_args()

class sonar():
	def __init__(self, db):
		self.db = sonardb.sonarDB(db)

	def add_genome(self, *fnames, paranoid=False):
		for fname in fnames:
			for record in SeqIO.parse(fname, "fasta"):
				self.db.add_genome(record.id, record.description, str(record.seq), paranoid)
		snr.db.commit()


if __name__ == "__main__":
	args = parse_args()
	snr = sonar(args.db)

	#add sequences
	if args.tool == "add":
		snr.add_genome(*args.fasta, paranoid=args.paranoid)

	#view data
	if args.tool == "view":
		rows = [x for x in snr.db.iter_rows(args.branch + "_view")]
		if len(rows) == 0:
			print("*** no data ***")
		else:
			writer = csv.DictWriter(sys.stdout, rows[0].keys(), lineterminator=os.linesep)
			writer.writeheader()
			writer.writerows(snr.db.iter_rows(args.branch + "_view"))
