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
		self.db = db
		self.dbobj = sonardb.sonarDB(self.db)
		self.gff = None
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

	def writefile(self, fname, *content):
		with open(fname, "w") as handle:
			handle.write("".join(content))

	@property
	def dna_var_regex(self):
		'''
		Provides a regex matching to valid dna variant expressions.
		'''
		if self.__dna_var_regex is None:
			allowed_letters = "[" + "".join(self.iupac_nt_code.keys()) + "]"
			self.__dna_var_regex = re.compile("^(?:(?:del:[0-9]+:[0-9]+)|(?:" + allowed_letters + "[0-9]+" + allowed_letters + "+))$")
		return self.__dna_var_regex

	@property
	def aa_var_regex(self):
		'''
		Provides a regex matching to valid protein variant expressions.
		'''
		if self.__aa_var_regex is None:
			allowed_symbols = "(?:(?:" + ")|(?:".join(self.dbobj.refgffObj.symbols) + "))"
			allowed_letters = "[" + "".join(self.iupac_aa_code.keys()).replace("-", "") + "-" + "]"
			self.__aa_var_regex = re.compile("^" + allowed_symbols + ":(?:(?:del:[0-9]+:[0-9]+)|(?:" + allowed_letters + "[0-9]+" + allowed_letters + "+))$")
		return self.__aa_var_regex

	@property
	def del_regex(self):
		'''
		Provides a regex matching to valid protein variant expressions.
		'''
		if self.__del_regex is None:
			self.__del_regex = re.compile(":?del:[0-9]+:[0-9]+$")
		return self.__del_regex

	@property
	def iupac_nt_code(self):
		'''
		Provides a dict of all IUPAC nucleotide one letter codes
		(key: one letter code, value: set of assigned one letter explicit codes).
		'''
		if self.__iupac_nt_code is None:
 			self.__iupac_nt_code = { "A": set("A"), "C": set("C"), "G": set("G"), "T": set("T"), "R": set("AG"), "Y": set("CT"), "S": set("GC"), "W": set("AT"), "K": set("GT"), "M": set("AC"), "B": set("CGT"), "D": set("AGT"), "H": set("ACT"), "V": set("ACG"), "N": set("ATGC") }
		return self.__iupac_nt_code

	@property
	def iupac_explicit_nt_code(self):
		'''
		Provides a dict of the IUPAC nucleotide one letter codes only containing letters coding for exactly one nucleotide
		(key: one letter code, value: set of assigned one letter explicit codes).
		'''
		if self.__iupac_explicit_nt_code is None:
 			self.__iupac_explicit_nt_code = set([ x for x in self.iupac_nt_code if len(self.iupac_nt_code[x]) == 1 ])
		return self.__iupac_explicit_nt_code

	@property
	def iupac_ambig_nt_code(self):
		'''
		Provides a dict of the IUPAC nucleotide one letter codes only containing letters coding for more than one nucleotides
		(key: one letter code, value: set of assigned one letter explicit codes).
		'''
		if self.__iupac_ambig_nt_code is None:
 			self.__iupac_ambig_nt_code = set([ x for x in self.iupac_nt_code if len(self.iupac_nt_code[x]) > 1 ])
		return self.__iupac_ambig_nt_code

	@property
	def iupac_aa_code(self):
		'''
		Provides a dict of all IUPAC amino acid one letter codes
		(key: one letter code, value: set of assigned one letter explicit codes).
		'''
		if self.__iupac_aa_code is None:
			self.__iupac_aa_code = { "A": set("A"), "R": set("R"), "N": set("N"), "D": set("D"), "C": set("C"), "Q": set("Q"), "E": set("E"), "G": set("G"), "H": set("H"), "I": set("I"), "L": set("L"), "K": set("K"), "M": set("M"), "F": set("F"), "P": set("P"), "S": set("S"), "T": set("T"), "W": set("W"), "Y": set("Y"), "V": set("V"), "U": set("U"), "O": set("O") }
			self.__iupac_aa_code['X'] = set(self.__iupac_aa_code.keys())
			self.__iupac_aa_code.update({"B": set("DN"), "Z": set("EQ"), "J": set("IL"), "Φ": set("VILFWYM"), "Ω": set("FWYH"), "Ψ": set("VILM"), "π": set("PGAS"), "ζ": set("STHNQEDKR"), "+": set("KRH"), "-": set("DE") })
		return self.__iupac_aa_code

	@property
	def iupac_explicit_aa_code(self):
		'''
		Provides a dict of the IUPAC amino acid one letter codes only containing letters coding for exactly one amino acid
		(key: one letter code, value: set of assigned one letter explicit codes).
		'''
		if self.__iupac_explicit_aa_code is None:
 			self.__iupac_explicit_aa_code = set([ x for x in self.iupac_aa_code if len(self.iupac_aa_code[x]) == 1 ])
		return self.__iupac_explicit_aa_code

	@property
	def iupac_ambig_aa_code(self):
		'''
		Provides a dict of the IUPAC amino acid one letter codes only containing letters coding for more than one amino acid
		(key: one letter code, value: set of assigned one letter explicit codes).
		'''
		if self.__iupac_ambig_aa_code is None:
 			self.__iupac_ambig_aa_code = set([ x for x in self.iupac_aa_code if len(self.iupac_aa_code[x]) > 1 ])
		return self.__iupac_ambig_aa_code

	def pinpoint_mutation(self, mutation, code):
		'''
		Returns a set of explicit mutations based on the given, possibly ambiguous mutation definition.
		The mutation definition must follow the covsonar nomenclature system.
		'''
		# extract ALT call from mutation profile
		match = self.__terminal_letters_regex.search(mutation)
		if not match:
			return mutation
		match = match.group(0)

		# resolve ambiguities
		options = []
		for m in match:
			options.append(code[m])

		# generate the set of explicit mutations
		orig_stat = mutation[:-len(match)]
		return set([mutation] + [ orig_stat + "".join(x) for x in itertools.product(*options) ])

	def filter_ambig(self, profile, explicit_code, keep=None):
		'''
		Returns a mutation profile that do not include any ambiguous SNP anymore.
		Mutations listed in keep will be not excluded.
		'''
		out = []
		keep = set(keep) if keep else set()
		for mutation in list(filter(None, profile.split(" "))):
			if mutation in keep or self.del_regex.search(mutation):
				out.append(mutation)
				continue
			match = self.__terminal_letters_regex.search(mutation)
			if match and len(match.group(0)) == 1 and  match.group(0) not in explicit_code:
				continue
			out.append(mutation)
		return " ".join(out)

	def add(self, *fnames, cpus=1, paranoid=True):
		'''
		Adds genome sequence(s) from the given FASTA file(s) to the database.
		'''
		# split fasta files to single entry temporary files
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

			# execute adding method of the database module on the generated files
			db = sonardb.sonarDB(self.db)
			fnames = [ os.path.join(tmpdirname, x) for x in os.listdir(tmpdirname) if x.endswith(".fasta") ]
			r = Parallel(n_jobs=cpus, verbose=10)(delayed(db.add_genome_from_fasta)(fname) for fname in fnames)
			if paranoid:
				r = Parallel(n_jobs=cpus, verbose=10)(delayed(self.paranoid_check)(fname) for fname in fnames)

	def paranoid_check(self, fname):
		record = SeqIO.read(fname, "fasta")
		acc = record.id
		descr = record.description
		seq = str(record.seq.upper())
		restored_seq = self.restore([acc], aligned=False, seq_return=True)
		if seq != restored_seq:
			sys.exit("'good that I am paranoid' error: " + acc + " origninal and those restored from the database do not match.")

	def create_profile_clause(self, profile, dna=True, exclusive=False, negate=False):
		'''
		Adds genome sequence(s) from the given FASTA file(s) to the database.
		'''
		# configuring
		condition = " " if not negate else " NOT "
		config = {
				   "dna": {
							"field": "dna_profile",
					  		"code": self.iupac_nt_code,
					  		"explicit_code": self.iupac_explicit_nt_code
						  },
					"aa": {
					        "field": "aa_profile",
					        "code": self.iupac_aa_code,
					        "explicit_code": self.iupac_explicit_aa_code
					}
				  }

		# generating clause
		clause = []
		vars = set(profile)
		dna_profile_length = 1
		aa_profile_length = 1
		for var in vars:
			if self.isdnavar(var):
				dna_profile_length += 1 + len(var)
				conf = config['dna']
			else:
				aa_profile_length += 1 + len(var)
				conf = config['aa']
			for v in self.pinpoint_mutation(var, conf['code']):
				clause.append(conf['field'] + condition + "LIKE '% " + v + " %'")

		clause = " AND ".join(clause)

		# add exclusiveness condition not allowing additional mutations
		if exclusive:
			if dna_profile_length > 1 and aa_profile_length > 1:
				sys.exit("input error: exclusive profiles must be defined on dna OR protein level only.")
			if dna_profile_length > 0:
				clause += " AND length(" + conf['field'] + ") = " + str(dna_profile_length)
			elif aa_profile_length > 0:
				clause += " AND length(" + conf['field'] + ") = " + str(aa_profile_length)
		return clause

	def isdnavar(self, var):
		'''
		Returns True if a variant definition is a valid dna variant expression.
		'''
		match = self.dna_var_regex.match(var)
		if match:
			return True
		return False

	def isaavar(self, var):
		'''
		Returns True if a variant definition is a valid aa variant expression.
		'''
		match = self.aa_var_regex.match(var)
		if match:
			return True
		return False

		def isdel(self, var):
			'''
			Returns True if a variant definition is a valid aa variant expression.
			'''
			match = self.aa_var_regex.match(var)
			if match:
				return True
			return False

	def isvalidvar(self):
		pass

	def match(self, include_profiles=None, exclude_profiles=None, accessions=None, lineages=None, exclusive=False, ambig=False, show_sql=False):
		'''
		Provides mutation profile matching against sequences in the database.
		'''

		# For developers:
		# Aim is to bring all profile definitions and filters to one and the same
		# sqllite query to let the database do the actual work

		clause = []
		vals =[]

		#sanity check:
		check = []
		if include_profiles:
			check += [item for sublist in include_profiles for item in sublist]
		if exclude_profiles:
			check += [item for sublist in exclude_profiles for item in sublist]
		nonvalid = [ x for x in check if not self.isdnavar(x) and not self.isaavar(x) ]
		if nonvalid:
			sys.exit("input error: Non-valid variant expression(s) entered: " + ", ".join(nonvalid))

		# adding conditions of profiles to include to where clause
		if include_profiles:
			includes = []
			for profile in include_profiles:
				includes.append(self.create_profile_clause(profile, exclusive=exclusive, negate=False))
			if len(includes) > 1:
				includes = [ "(" + x + ")" if " AND " in x else x for x in includes ]
			clause.append(" OR ".join(includes))
			if len(includes) > 1:
				clause[-1] = "(" + clause[-1] + ")"

		# adding conditions of profiles to exclude to where clause
		if exclude_profiles:
			excludes = []
			for profile in exclude_profiles:
				excludes.append(self.create_profile_clause(profile, exclusive=exclusive, negate=True))
			excludes = [ "(" + x + ")" if " AND " in x else x for x in excludes ]
			clause.append(" AND ".join(excludes))

		# adding accession and lineage based conditions
		if accessions:
			clause.append("accession IN (" + " , ".join(['?'] * len(accessions)) + ")")
			vals.extend(accessions)
		if lineages:
			clause.append("lineages IN (" + " ,".join(['?'] * len(lineages))  + ")")
			vals.extend(lineages)

		# executing the query and storing results
		clause = " AND ".join(clause)
		rows = [x for x in self.dbobj.select("essence", whereClause=clause, valList=vals, show_sql=show_sql)]

		# remove ambiguities from database profiles if wished
		if not ambig:
			keep = [item for sublist in include_profiles for item in sublist] if include_profiles else None
			for i in range(len(rows)):
				rows[i]['dna_profile'] = self.filter_ambig(rows[i]['dna_profile'], self.iupac_explicit_nt_code, keep)
				rows[i]['aa_profile'] = self.filter_ambig(rows[i]['aa_profile'], self.iupac_explicit_aa_code, keep)

		self.rows_to_csv(rows)

	def restore(self, accs, aligned=False, seq_return=False):
		rows = [x for x in self.dbobj.select("dna_view", whereClause=" ".join(["accession = ?"] * len(accs)), valList=accs, orderby="start DESC")] # it's crucial to sort in desc to safely insert potential insertions
		if rows:
			gap = "-" if aligned else ""
			refseq = list(self.dbobj.refseq)
			qryseq = refseq[:]
			for row in rows:
				s = row['start']
				if row['ref'] != refseq[s]:
					sys.exit("data error: data is inconsistence (" + row['ref']+ " expected at position " + str(s+1) + " of the reference sequence, got " + refseq[s] + ").")
				qryseq[s] = gap if not row['alt'] else row['alt']
				if aligned and len(row['alt']) > 1:
					refseq[s] +=  "-" * (len(row['alt'])-1)
			if seq_return:
				return "".join(qryseq)
			if aligned:
				print(">" + self.dbobj.refdescr)
				print("".join(refseq))
			print(">" + rows[0]['description'])
			print("".join(qryseq))

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

	#add sequences
	if args.tool == "add":
		snr.add(*args.fasta, cpus=args.cpus)

	#show
	if args.tool == "match":
		snr.match(include_profiles=args.include, exclude_profiles=args.exclude, lineages=args.lineage, accessions=args.acc, exclusive=args.exclusive, ambig=args.ambig)

	#restore alignment
	if args.tool == "restore":
		snr.restore(args.acc, args.align)

	#view data
	if args.tool == "view":
		rows = [x for x in snr.dbobj.iter_rows(args.branch + "_view")]
		snr.rows_to_csv(rows)
