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

	# create the parser for the "match" command
	parser_match = subparsers.add_parser('match', parents=[general_parser], help='get mutations profiles for given accessions.')
	parser_match.add_argument('--profile', metavar="STR", help="match genomes sharing a given mutation profile", type=str, nargs="+", default=None)
	parser_match.add_argument('--lineage', metavar="STR", help="match genomes of a given pangolin lineage only", type=str, nargs="+", default=None)
	parser_match.add_argument('--acc', metavar="STR", help="match specific genomes defined by acessions only", type=str, nargs="+", default=None)
	parser_match.add_argument('--exclusive', help="do not allow additional mutations", action="store_true")
	parser_match.add_argument('--count', help="count only matching genomes", action="store_true")
	parser_match.add_argument('--logic', help="decide between OR or AND  logic when matching profile mutations (default: AND logic)", choices=["OR", "or", "AND", "and"], default="AND")
	parser_match.add_argument('--negate', help="show genomes NOT matching the mutation profile", action="store_true")
	parser_match.add_argument('--ambig', help="include ambiguos sites when reporting profiles (no effect when --count is used)", action="store_true")

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
		self.__iupac_nt_code = None
		self.__iupac_aa_code = None
		self.__iupac_explicit_nt_code = None
		self.__iupac_explicit_aa_code = None
		self.__terminal_letters_regex = re.compile("[A-Z]+$")

	def writefile(self, fname, *content):
		with open(fname, "w") as handle:
			handle.write("".join(content))

	@property
	def iupac_nt_code(self):
		if self.__iupac_nt_code is None:
 			self.__iupac_nt_code = { "A": set("A"), "C": set("C"), "G": set("G"), "T": set("T"), "R": set("AG"), "Y": set("CT"), "S": set("GC"), "W": set("AT"), "K": set("GT"), "M": set("AC"), "B": set("CGT"), "D": set("AGT"), "H": set("ACT"), "V": set("ACG"), "N": set("ATGC") }
		return self.__iupac_nt_code

	@property
	def iupac_explicit_nt_code(self):
		if self.__iupac_explicit_nt_code is None:
 			self.__iupac_explicit_nt_code = set([ x for x in self.iupac_nt_code if len(self.iupac_nt_code[x]) == 1 ])
		return self.__iupac_explicit_nt_code

	@property
	def iupac_aa_code(self):
		if self.__iupac_aa_code is None:
			self.__iupac_aa_code = { "A": set("A"), "R": set("R"), "N": set("N"), "D": set("D"), "C": set("C"), "Q": set("Q"), "E": set("E"), "G": set("G"), "H": set("H"), "I": set("I"), "L": set("L"), "K": set("K"), "M": set("M"), "F": set("F"), "P": set("P"), "S": set("S"), "T": set("T"), "W": set("W"), "Y": set("Y"), "V": set("V"), "U": set("U"), "O": set("O") }
			self.__iupac_aa_code['X'] = set(self.__iupac_aa_code.keys())
			self.__iupac_aa_code.update({"B": set("DN"), "Z": set("EQ"), "J": set("IL"), "Φ": set("VILFWYM"), "Ω": set("FWYH"), "Ψ": set("VILM"), "π": set("PGAS"), "ζ": set("STHNQEDKR"), "+": set("KRH"), "-": set("DE") })
		return self.__iupac_aa_code

	@property
	def iupac_explicit_aa_code(self):
		if self.__iupac_explicit_aa_code is None:
 			self.__iupac_explicit_aa_code = set([ x for x in self.iupac_aa_code if len(self.iupac_aa_code[x]) == 1 ])
		return self.__iupac_explicit_aa_code

	def pinpoint_mutation(self, mutation, code):
		match = self.__terminal_letters_regex.search(mutation)
		if not match:
			return mutation
		match = match.group(0)
		options = []
		for m in match:
			options.append(code[m])
		orig_stat = mutation[:-len(match)]
		return set([ orig_stat + "".join(x) for x in itertools.product(*options) ])

	def filter_ambig(self, profile, explicit_code):
		keep = []
		for mutation in list(filter(None, profile.split(" "))):
			match = self.__terminal_letters_regex.search(mutation)
			if match is None or match.group(0) in explicit_code:
				keep.append(mutation)
		return " ".join(keep)

	def add(self, *fnames, cpus=1):

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
			r = Parallel(n_jobs=cpus, verbose=10)(delayed(db.add_genome_from_fasta)(fname) for fname in fnames)
			#for fname in fnames:
			#	db.add_genome_from_fasta(fname)

	def match(self, profiles=None, accessions=None, lineages=None, exclusive=False, profile_logic="AND", negate=False, ambig=False, show_sql=True):
		profile_logic = " " + profile_logic.upper() + " "
		b = "()" if profile_logic != " OR " else ("", "")

		where = []
		vals =[]

		# adding profile condition to where clause
		if profiles:
			dna_profiles = [self.pinpoint_mutation(x, self.iupac_nt_code) for x in profiles if ":" not in x ]
			aa_profiles = [self.pinpoint_mutation(x, self.iupac_aa_code) for x in profiles if ":" in x ]

			## sanity checks of options set
			if exclusive and dna_profiles and aa_profiles:
				sys.exit("input error: --exclusive option can not applied when profile contains explicit nucleotide mutation information.")
			if exclusive and or_op:
				sys.exit("input error: --excusive and --or are mutually eclusive options")

			## assembling dna profiles
			for dna_profile in dna_profiles:
				if len(dna_profile) == 0:
					where.append("dna_profile LIKE '% " + dna_profile + " %'")
				else:
					where.append(b[0] + " OR ".join(["dna_profile LIKE '% " + x + " %'" for x in dna_profile]) + b[1])

			## assembling aa profiles
			for aa_profile in aa_profiles:
				if len(dna_profile) == 0:
					where.append("aa_profile LIKE '% " + aa_profile + " %'")
				else:
					where.append(b[0] + " OR ".join(["aa_profile LIKE '% " + x + " %'" for x in aa_profile]) + b[1])

			## generating logic-dependent syntax
			where = [profile_logic.join(where)]

			if profile_logic == " OR " and not where[0].endswith(")") and any((accessions, lineages)):
				where[0] = "(" + where[0] + ")"

			## adding condition for exclusive matching to where clause
			if dna_profiles and exclusive:
				where.append("length(dna_profile) = ?")
				vals.append(len(" " + " ".join([ x[0] for x in dna_profiles ]) + " "))

			if aa_profiles and exclusive:
				where.append("length(aa_profile) = ?")
				vals.append(len(" " + " ".join(aa_profiles) + " "))

		if accessions:
			where.append("accession IN '(" + " ,".join(accessions) + ")'")
		if lineages:
			where.append("lineages IN '(" + " ,".join(lineages) + ")'")

		where = " AND ".join(where)
		rows = [x for x in self.dbobj.select("essence", whereClause=where, valList=vals, show_sql=show_sql)]
		if negate:
			where = "accession NOT IN (" + " ,".join(['?'] * len(rows)) + ")"
			vals = [ x['accession'] for x in rows ]
			rows = [x for x in self.dbobj.select("essence", whereClause=where, valList=vals, show_sql=show_sql)]

		if not ambig:
			for i in range(len(rows)):
				rows[i]['dna_profile'] = self.filter_ambig(rows[i]['dna_profile'], self.iupac_explicit_nt_code)
				rows[i]['aa_profile'] = self.filter_ambig(rows[i]['aa_profile'], self.iupac_explicit_aa_code)

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
		snr.add(*args.fasta, cpus=args.cpus)

	#show
	if args.tool == "match":
		snr.match(profiles=args.profile, lineages=args.lineage, accessions=args.acc, profile_logic=args.logic, exclusive=args.exclusive, negate=args.negate, ambig=args.ambig)

	#view data
	if args.tool == "view":
		rows = [x for x in snr.dbobj.iter_rows(args.branch + "_view")]
		snr.rows_to_csv(rows)
