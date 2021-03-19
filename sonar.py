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
	parser_add.add_argument('-f', '--fasta', metavar="FILE", help="fasta file(s) containing DNA sequences to add", type=str, nargs="+")
	parser_add.add_argument('--paranoid', help="activate checks on seguid sequence collisions", action="store_true")

	# create the parser for the "show" command
	parser_show = subparsers.add_parser('show', parents=[general_parser], help='get mutation profiles for given accessions.')
	parser_show.add_argument('--profile', metavar="STR", help="genomes to show have to contain variants", type=str, nargs="+", default=None)
	parser_show.add_argument('--lineage', metavar="STR", help="genomes to show have to be assigned to the given pangolin lineage", type=str, nargs="+", default=None)
	parser_show.add_argument('--acc', metavar="STR", help="only consider the given genome acessions", type=str, nargs="+", default=None)
	parser_show.add_argument('--exclusive', help="do not allow additional variantions", action="store_true")
	parser_show.add_argument('--exact', help="exact deletion handling", action="store_true")

	# create the parser for the "view" command
	parser_add = subparsers.add_parser('view', parents=[general_parser], help='view database content.')
	parser_add.add_argument('-v', metavar="FILE", help="complete data view (default: dna)", choices=['dna', 'prot', 'profiles'], default="profiles")

	# version
	parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)

	return parser.parse_args()

class sonar():
	def __init__(self, db, gff=None):
		self.db = db
		self.dbobj = sonardb.sonarDB(self.db)
		self.gff = None

	def writefile(self, fname, *content):
		with open(fname, "w") as handle:
			handle.write("".join(content))

	def add(self, *fnames):

		# split fasta files to single entry files

		with tempfile.TemporaryDirectory(dir=os.getcwd()) as tmpdirname:
			i = 0
			entry = []
			for fname in fnames:
				with open(fname, "r") as inhandle:
					for line in inhandle:
						if line.startswith(">"):
							if entry:
								self.writefile(os.path.join(tmpdirname, str(i) + ".fasta"), *entry)
								i += 1
							entry = [line]
						else:
							entry.append(line)
			if entry:
				self.writefile(os.path.join(tmpdirname, str(i) + ".fasta"), *entry)

			db = sonardb.sonarDB(self.db)

			fnames = [ os.path.join(tmpdirname, x) for x in os.listdir(tmpdirname) if x.endswith(".fasta") ]
			for fname in fnames:
				db.add_genome_from_fasta(fname)

	def show(self, profiles=None, accessions=None, lineages=None, exclusive=False):
		where = []
		vals =[]
		if profiles:
			dna_profiles = [x for x in profiles if ":" not in x ]
			aa_profiles = [x for x in profiles if ":" in x ]
			if exclusive and dna_profiles and aa_profiles:
				sys.exit("input error: to excusively screen for mutation profiles profile has to be defined on nucleotide OR protein level.")
			for dna_profile in dna_profiles:
				where.append("dna_profile LIKE '% " + dna_profile + " %'")
			if dna_profiles and exclusive:
				where.append("length(dna_profile) = ?")
				vals.append(len(" " + " ".join(dna_profiles) + " "))
			for aa_profile in aa_profiles:
				where.append("aa_profile LIKE '% " + aa_profile + " %'")
			if aa_profiles and exclusive:
				where.append("length(aa_profile) = ?")
				vals.append(len(" " + " ".join(aa_profiles) + " "))
		if accessions:
			where.append("accession IN '(" + " ,".join(accessions) + ")'")
			vals.extend(accessions)
		if lineages:
			where.append("lineages IN '(" + " ,".join(lineages) + ")'")
			vals.extend(lineages)
		where = " AND ".join(where)
		rows = [x for x in self.dbobj.select("profile", whereClause=where, valList=vals)]
		self.rows_to_csv(rows)

	def rows_to_csv(self, rows):
		if len(rows) == 0:
			print("*** no data ***")
		else:
			writer = csv.DictWriter(sys.stdout, rows[0].keys(), lineterminator=os.linesep)
			writer.writeheader()
			writer.writerows(rows)

if __name__ == "__main__":
	args = parse_args()
	snr = sonar(args.db)
	# ray.init(num_cpus=args.cpus, include_dashboard=False)

	#add sequences
	if args.tool == "add":
		snr.add(*args.fasta)

	#show
	if args.tool == "show":
		snr.show(profiles=args.profile, lineages=args.lineage, accessions=args.acc, exclusive=args.exclusive)

	#view data
	if args.tool == "view":
		rows = [x for x in snr.dbobj.iter_rows(args.branch + "_view")]
		snr.rows_to_csv(rows)
