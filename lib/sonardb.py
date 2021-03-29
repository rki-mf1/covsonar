#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

import os
import re
import sys
import argparse
import sqlite3
from sqlite3 import Error
from Bio.SeqUtils.CheckSum import seguid
from Bio.Seq import Seq
from Bio import Align
from Bio import SeqIO
from Bio.Align import substitution_matrices
from packaging import version
import shutil
import base64
from collections import OrderedDict
import pickle
from tqdm import tqdm
from urllib.parse import quote as urlquote
from math import floor
from tempfile import TemporaryDirectory, mkdtemp
import traceback
import difflib
import itertools
import gc

class sonarCDS(object):
	'''
	This class provides properties and functions for each coding sequence
	of the reference annotation
	'''
	def __init__(self, symbol, start, end, strand, seq, locus=None, translation_table=1):
		'''
		following parameters are needed to initialize the class:
		symbol = gene symbol
		start = CDS start coordinate on the reference genome (0-based, inclusive)
		end = CDS end coordinate on the reference genome (0-based, exclusive)
		strand = coding strand ("+" or "-")
		seq = coding nucleotide sequence
		locus = gene locus accession number
		translation_table = translation table to use (default: 1)
		'''
		self.symbol = symbol
		self.locus = locus
		self.start = min(start, end) # inclusive
		self.end = max(start, end) # exclusive
		self.strand = "+" if strand == "+" or int(strand) == 1 else "-"
		self.dna = seq
		self.translation_table = translation_table
		self.__aa = None

	@property
	def aa(self):
		'''
		Property providing the aminoc acid code of the CDS from start to the
		first stop codon.

		>>> cds=sonarCDS("Eno", 155, 170, "+", "ATGTTATGAATGGCC")
		>>> cds.aa
		'ML'
		'''
		if self.__aa is None:
			l = len(self.dna)
			if l%3 == 1:
				l = -1
			elif l%3 == 2:
				l = -2
			self.__aa = str(Seq.translate(self.dna[:l], table=self.translation_table, to_stop=True))
		return self.__aa

	@property
	def coords(self):
		'''
		Property providing a tuple of 0-based start (inclusive) and end (exclusive)
		coordinate  of the CDS in the reference genome. Start coordinate is alsways
		smaller than end coordinate.

		>>> cds=sonarCDS("Eno", 155, 170, "+", "ATGTTATGAATGGCC")
		>>> cds.coords
		(155, 170)
		'''
		return self.start, self.end

	@property
	def range(self):
		'''
		Property providing a range from start and end coordinate of the CDS in the reference genome.

		>>> cds=sonarCDS("Eno", 155, 170, "+", "ATGTTATGAATGGCC")
		>>> cds.range
		range(155, 170)
		'''
		return range(self.start, self.end)

class sonarGFF(object):
	'''
	This class provides properties and functions of CDS annotations extracted
	from a given file.
	'''
	def __init__(self, gff3, fna, translation_table=1):
		'''
		following parameters are needed to initialize the class:
		gff3 = gff3 file name containing the annotation
		fna = fasta file name containing the reference genome sequence (only one sequence allowed)
		translation_table = translation_table = translation table to use (default: 1)
		'''
		self.translation_table = translation_table
		self.cds = self.process_gff3(gff3, fna)
		self.coords = { x.symbol: (x.start, x.end) for x in self.cds }
		self.symbols = [x.symbol for x in self.cds ]

	def iscds(self, x, y=None):
		'''
		Checks if reference genome coordinate X [to y] is in an annotated CDS.
		Both coordinates have to be 0-based, x is inclusive, y exclusive.

		>>> os.chdir(os.path.dirname(os.path.realpath(__file__)))
		>>> gff=sonarGFF("doctest.gff3", "doctest.fna")
		>>> gff.iscds(2)
		False
		>>> gff.iscds(7)
		True

		'''
		x, y = sorted([x, y]) if y is not None else (x, x)
		for start, end in self.coords.values():
			if x >= start and y < end:
				return True
		return False

	def convert_pos_to_dna(self, protein_symbol, x):
		'''
		Returns a reference genome coordinate (0-based start of the codon)
		based on a given 0-based protein coordinate (x).

		>>> os.chdir(os.path.dirname(os.path.realpath(__file__)))
		>>> gff=sonarGFF("doctest.gff3", "doctest.fna")
		>>> gff.convert_pos_to_dna("ORF1", 2)
		8
		'''
		if x > 0:
			x = x*3 - 2
		return x + self.coords[protein_symbol][0]

	def convert_pos_to_aa(self, x, y=None):
		'''
		Based on a given 0-based reference genome coordinate x [and y] a tuple of
		respective protein symbol and 0-based protein coordinates are returned.
		If x is not in an annotated CDS None is returned.
		If multiple CDS cover the respective genome position(s), one is randomly
		selected.

		>>> os.chdir(os.path.dirname(os.path.realpath(__file__)))
		>>> gff=sonarGFF("doctest.gff3", "doctest.fna")
		>>> gff.convert_pos_to_aa(6,10)
		('ORF1', 0, 1)
		'''
		x, y = sorted([x, y]) if y is not None else (x, x)

		for symbol, coords in self.coords.items():
			if x >= coords[0] and y < coords[1]:
				x = floor((x - coords[0])/3)
				if y is None:
					return (symbol, x)
				y = floor((y - 1 - coords[0])/3)
				return (symbol, x, y)
		return None

	def process_gff3(self, gff, fna):
		'''
		Parse a GFF3 formatted file to extract all CDS. A list of sonarCDS objects
		is returned.

		>>> os.chdir(os.path.dirname(os.path.realpath(__file__)))
		>>> gff=sonarGFF("doctest.gff3", "doctest.fna")

		>>> gff.coords
		{'ORF1': (4, 13)}
		'''
		symbol_regex = re.compile("gene=([^;]+)(?:;|$)")
		locus_regex = re.compile("locus_tag=([^;]+)(?:;|$)")

		record = SeqIO.read(fna, "fasta")
		gseq = str(record.seq).upper()

		with open(gff, "r") as handle:
			cds = set()
			for line in handle:
				if line.startswith("#") or len(line.strip()) == 0:
					continue
				fields = line.rstrip("\r\n").split("\t")
				if fields[2] == "CDS":
					symbol = symbol_regex.search(fields[-1])
					symbol = symbol.groups(1)[0] if symbol else None
					locus = locus_regex.search(fields[-1])
					locus = locus.groups(1)[0] if locus else None
					strand = fields[6]
					s = int(fields[3])-1
					e = int(fields[4])
					dna = gseq[s:e]
					if strand == "-":
						dna = str(Seq.reverse_complement(dna))
					cds.add((symbol, s, e, strand, dna, locus, self.translation_table))
		return [sonarCDS(*x) for x in cds]

class sonarALIGN(object):
	'''
	This class provides alignment functions and properties to a given query and target
	sequence. Additional protein-based functionalities are provided if a sonarGFF object
	is provided.
	'''
	def __init__(self, query, target, sonarGFFObj = None, aligned = False, scoring=(-16,-4,16,-4,0,0)):
		'''
		following parameters are needed to initialize the class:
		query = query nucleotide sequence to be algned to the target sequence
		target = target nucleotide sequence
		sonarGFFObj = sonarGFF object to consider cds annotation and provide protein level based functions
		'''
        if aligned:
            self.aligned_query = query.upper()
            self.aligned_target = target.upper()
        else:   
            self.aligned_query, self.aligned_target = self.align_dna(query.upper(), target.upper(), *scoring)

		self.gff = sonarGFFObj if sonarGFFObj else None
		self.indel_regex = re.compile(".-+")
		self.codon_regex = re.compile(".-*.-*.-*")
		self.__dnadiff = None
		self.__aadiff = None
		self.__target_coords_matrix = None

	@property
	def dnadiff(self):
		if self.__dnadiff is None:
			self.__dnadiff = [ x for x in self.iter_dna_vars() ]
		return self.__dnadiff

	@property
	def aadiff(self):
		if self.__aadiff is None:
			self.__aadiff = [ x for x in self.iter_aa_vars() ]
		return self.__aadiff

	@property
	def target_coords_matrix(self):
		if self.__target_coords_matrix is None:
			self.__target_coords_matrix = [len(x.group()) for x in re.finditer(".-*", self.aligned_target)]
		return self.__target_coords_matrix
    
    @staticmethod
	def align_dna(self, query, target, open_gap_score= -16, extend_gap_score = -4, end_extend_gap_score = -16,
				  end_gap_score = -4, target_end_gap_score = 0, query_end_gap_score = 0,
				  matrixfile = os.path.join(os.path.dirname(os.path.realpath(__file__)), "EDNAFULL")):
		'''
		align two nucleotide sequences and return a tuple of the aligned sequences.
		EDNAFULL is used as default matrix.

		>>> algn = sonarALIGN("ATTGTTGTTATAATGGCCAGTT", "ATTGATGTTATGAATGGCCTTT", scoring=(-10, -0.5, -10, -0.5, 0, 0))
		>>> algn.align_dna("ATTGTTGTTATAATGGCCAGTT", "ATTGATGTTATGAATGGCCTTT", -10, -0.5, -10, -0.5, 0, 0)
		('ATTGTTGTTAT-AATGGCCAGTT-', 'ATTGATGTTATGAATGGCC--TTT')

		Instead of using this function, the object properties aligned_query and
		aligned_target should be used:

		>>> algn.aligned_query
		'ATTGTTGTTAT-AATGGCCAGTT-'
		>>> algn.aligned_target
		'ATTGATGTTATGAATGGCC--TTT'

		'''
		aligner = Align.PairwiseAligner()
		aligner.mode = 'global'
		aligner.substitution_matrix = substitution_matrices.read(matrixfile)
		aligner.open_gap_score = open_gap_score
		aligner.extend_gap_score = extend_gap_score
		aligner.end_extend_gap_score = end_extend_gap_score
		aligner.end_gap_score = end_gap_score
		aligner.target_end_gap_score = target_end_gap_score
		aligner.query_end_gap_score = query_end_gap_score
		alignments = aligner.align(query, target)
		best_alignment = None
		for alignment in alignments:
			if not best_alignment or best_alignment.score < alignment.score:
				best_alignment = alignment
		best_alignment = str(best_alignment).split("\n")
		return best_alignment[0], best_alignment[2]

	def align_aa(self, query, target):
		'''
		this function is not in use right now, and should be revised to mimick EMBOSS NEEDLE before.
		'''
		aligner = Align.PairwiseAligner()
		aligner.substitution_matrix = substitution_matrices.load("PAM250")
		aligner.mode = 'global'
		aligner.open_gap_score = -10
		aligner.extend_gap_score = -0.5
		aligner.target_end_gap_score = 0.0
		aligner.query_end_gap_score = 0.0
		alignment = str(sorted(aligner.align(query, target), key=lambda x: x.score, reverse=True)[0]).split("\n")
		return alignment[0], alignment[2]

	def real_pos(self, x):
		'''
		converts alignment position to the real position in the target sequence.
		Positions are meant to be 0-based.

		>>> algn = sonarALIGN("ATTGTTGTTATAATGGCCAGTT", "ATTGATGTTATGAATGGCCTTT", scoring=(-10, -0.5, -10, -0.5, 0, 0))
		>>> algn.real_pos(22)
		20
		'''
		return x - self.aligned_target[:x+1].count("-")

	def align_pos(self, x):
		'''
		converts a target sequence position to the respective position in the alignment.
		Positions are meant to be 0-based.

		>>> algn = sonarALIGN("ATTGTTGTTATAATGGCCAGTT", "ATTGATGTTATGAATGGCCTTT", scoring=(-10, -0.5, -10, -0.5, 0, 0))
		>>> algn.align_pos(20)
		22
		'''
		return sum(self.target_coords_matrix[:x])

	def iter_dna_vars(self):
		'''
		generating an iterator returning tuples for each variant alignment position
		consisting of target nucleotide, query nucleotide(s), target start position, target end position, None, None
		Last two None elements exist to harmonize output syntax between this and the iter_aa_vars function.

		>>> algn = sonarALIGN("ATTGTTGTTATATGGCCAGTT", "ATTGATGTTATGAATGGCCTTT", scoring=(-10, -0.5, -10, -0.5, 0, 0))
		>>> for x in algn.iter_dna_vars():
		... 	print(x)
		('C', 'CAG', 18, None, None, None)
		('A', 'T', 4, None, None, None)
		('G', '', 11, None, None, None)
		('A', '', 12, None, None, None)
		('T', '', 21, None, None, None)

		'''
		target = self.aligned_target
		query = self.aligned_query

		# insertions
		for match in self.indel_regex.finditer(target):
			s = self.real_pos(match.start())
			yield match.group()[0], query[match.start():match.end()], s, None, None, None

		# deletions and snps
		for i, pair in enumerate(zip(target, query)):
			if pair[0] != "-" and pair[0] != pair[1]:
				l = len(pair[1])
				s = self.real_pos(i)
				e = None if l == 1 else s + l
				yield pair[0], pair[1].replace("-", ""), s, e, None, None

	def iter_aa_vars(self):
		'''
		generating an iterator returning tuples for each variant alignment position that affects a protein sequence
		(based on given cds annotation provided by the sonarGFF object) consisting of target amino acid, query amino acid(s),
		target start position, target end position, protein symbol, locus accession.
		Last two None elements exist to harmonize output syntax between this and the iter_aa_vars function.

		>>> os.chdir(os.path.dirname(os.path.realpath(__file__)))
		>>> gff=sonarGFF("doctest.gff3", "doctest.fna")
		>>> algn = sonarALIGN("ATTGTTGTTATATGGCCAGTT", "ATTGATGTTATGAATGGCCTTT", gff)
		>>> for x in algn.iter_aa_vars():
		... 	print(x)
		('M', 'L', 0, None, 'ORF1', 'testseq_0001')
		('*', '', 2, None, 'ORF1', 'testseq_0001')

		'''
		if self.gff:
			for cds in self.gff.cds:
				if cds.strand == "+":
					s = self.align_pos(cds.start)
					e = self.align_pos(cds.end)
					query = self.aligned_query[s:e]
					target = self.aligned_target[s:e]
				else:
					s = self.align_pos(cds.start)
					e = self.align_pos(cds.end)
					query = str(Seq.reverse_complement(self.aligned_query[s:e]))
					target = str(Seq.reverse_complement(self.aligned_target[s:e]))

				for match in self.codon_regex.finditer(target):
					s = match.start()
					e = match.end()
					tcodon = match.group().replace("-", "")
					qcodon = query[match.start():match.end()]
					taa = self.translate(tcodon, cds.translation_table)
					if "-" in qcodon:
						yield taa, "", int(s/3), None, cds.symbol, cds.locus
						continue
					qaa = self.translate(qcodon, cds.translation_table)
					if qaa != taa:
						e = None if len(qaa) == 1 else int(e/3)
						yield taa, qaa, int(s/3), e, cds.symbol, cds.locus

	def translate(self, seq, translation_table=1):
		'''
		returns amino acid sequence based on a given nucleotide sequence and translation table.
		The given sequence is shortened if its length is not a multiple of 3.

		>>> algn = sonarALIGN("ATTGTTGTTATATGGCCAGTT", "ATTGATGTTATGAATGGCCTTT")
		>>> algn.translate("ATGTGAAA")
		'M*'

		'''
		l = len(seq)
		if l%3 == 1:
			l = -1
		elif l%3 == 2:
			l = -2
		return str(Seq.translate(seq[:l], table=translation_table))

class sonarDBManager():
	'''
	This class provides basic functions and properties to manage a given sonar database.
	At initialization, the database is created if necessary and a connection is
	established. All modyfying actions are grouped in a single transaction by default.
	The class should be always used in a context manager (with statement) to allow
	rollback at abnormal termination.
	'''
	def __init__(self, dbfile, timeout=-1, readonly=False):
		'''
		following parameters are needed to initialize the class:
		dbfile = database file (will be created if it does not exists)
		timeout = busy time out of the database
		readonly = allow read only access only
		'''
		self.dbfile = os.path.abspath(dbfile)
		self.connection = None
		self.cursor = None
		self.__timeout = timeout
		self.__mode = "ro" if readonly else "rwc"
		self.__uri = "file:" + urlquote(self.dbfile)

	def __enter__(self):
		if not os.path.isfile(self.dbfile) or os.stat(self.dbfile).st_size == 0:
			self.create_scheme()
		self.connection, self.cursor = self.connect()
		self.start_transaction()
		return self

	def __exit__(self, exc_type, exc_value, exc_traceback):
		if [exc_type, exc_value, exc_traceback].count(None) != 3:
			print("warning:", file=sys.stderr)
			print(traceback.format_exc(), file=sys.stderr)
			if self.__mode == "rwc":
				print("rollback", file=sys.stderr)
				self.rollback()
		elif self.__mode == "rwc":
			self.connection.commit()
		self.close()

	def __del__(self):
		if self.connection:
			self.close()

	def connect(self):
		'''
		Establishes a database connection and returns a tuple of connection and cursor
		'''
		con = sqlite3.connect(self.__uri + "?mode=" + self.__mode, self.__timeout, isolation_level = None, uri = True)
		con.row_factory = self.dict_factory
		cur = con.cursor()
		return con, cur

	def start_transaction(self):
		self.cursor.execute("BEGIN DEFERRED")

	def commit(self):
		self.connection.commit()

	def rollback(self):
		self.connection.rollback()

	def close(self):
		self.connection.close()

	def create_scheme(self):
		'''
		create database scheme.
		'''
		with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), "db.sqlite"), 'r') as handle:
			sql = handle.read()
		with sqlite3.connect(self.__uri + "?mode=rwc", uri = True) as con:
			con.executescript(sql)

	def select(self, table, fieldList=['*'], whereClause=[], whereVals=None, orderby=None, fetchone=False, print_sql=False):
		sql = "SELECT " + "".join(fieldList) + " FROM " + table
		if whereClause:
			sql += " WHERE " + whereClause
		else:
			whereVals = []
		if orderby:
			sql += " ORDER BY " + orderby
		if print_sql:
			print(sql)
		self.cursor.execute(sql, whereVals)
		if fetchone:
			return self.cursor.fetchone()
		return self.cursor.fetchall()

	def delete(self, table, whereClause=None, whereVals=[], print_sql=False):
		sql = "DELETE FROM " + table
		if whereClause:
			sql += " WHERE " + whereClause
		if print_sql:
			print(sql)
		self.cursor.execute(sql, whereVals)

	def insert(self, table, fieldList, valList, ignore=False, print_sql=False):
		ignore = "" if not ignore else "OR IGNORE "
		sql = "INSERT " + ignore + "INTO " + table + "(" + ", ".join(fieldList) + ") VALUES (" + ", ".join(['?']*len(valList)) + ")"
		if print_sql:
			print(sql)
		self.cursor.execute(sql, valList)
		return self.cursor.lastrowid

	def update(self, table, fieldList, valList, whereClause, whereVals, print_sql=False):
		sql = "UPDATE " + table + " SET " + ", ".join([ str(x) + "= ?" for x in fieldList ]) + " WHERE " + whereClause
		if print_sql:
			print(sql)
		self.cursor.execute(sql, valList + whereVals)

	@staticmethod
	def optimize(dbfile):
		with sqlite3.connect(dbfile) as con:
			con.executescript("VACUUM")

	@staticmethod
	def dict_factory(cursor, row):
		'''
		convert database result in list of dictionaries where keys are the respective column names.
		'''
		d = OrderedDict()
		for idx, col in enumerate(cursor.description):
			d[col[0]] = row[idx]
		return d


class sonarDB(object):
	'''
	This class provides sonarDB the actual functionalities.
	'''
	def __init__(self, dbfile, translation_table = 1):
		'''
		following parameters are needed to initialize the class:
		dbfile = database file (will be created if it does not exists)
		translation_table = translation database
		'''
		self.db = os.path.abspath(dbfile)
		self.moduledir = os.path.dirname(os.path.realpath(__file__))
		self.reffna = os.path.join(self.moduledir, "ref.fna")
		self.refgff = os.path.join(self.moduledir, "ref.gff3")
		self.translation_table = translation_table
		self.__refseq = None
		self.__refdescr = None
		self.__refgffObj = None
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

	# PROPERTIES ON DEMAND

	@property
	def refseq(self):
		'''
		Provides the built-in reference genome sequence.
		'''
		if not self.__refseq:
			with open(self.reffna, "r") as handle:
				self.__refseq = "".join([x.strip().upper() for x in handle.readlines()[1:]])
		return self.__refseq

	@property
	def refdescr(self):
		'''
		Provides the header information of the built-in reference genome.
		'''
		if not self.__refdescr:
			with open(self.reffna, "r") as handle:
				self.__refdescr = handle.readline().strip()[1:]
		return self.__refdescr

	@property
	def refgffObj(self):
		'''
		Provides the sonarGFF header based on the built-in reference genome
		annotation.
		'''
		if not self.__refgffObj:
			self.__refgffObj = sonarGFF(self.refgff, self.reffna, self.translation_table)
		return self.__refgffObj

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
			allowed_symbols = "(?:(?:" + ")|(?:".join(self.refgffObj.symbols) + "))"
			allowed_letters = "[" + "".join(self.iupac_aa_code.keys()).replace("-", "") + "-" + "]"
			self.__aa_var_regex = re.compile("^" + allowed_symbols + ":(?:(?:del:[0-9]+:[0-9]+)|(?:" + allowed_letters + "[0-9]+" + allowed_letters + "+))$")
		return self.__aa_var_regex

	@property
	def del_regex(self):
		'''
		Provides a regex matching to valid protein variant expressions.
		'''
		if self.__del_regex is None:
			allowed_symbols = "(?:(?:" + ")|(?:".join(self.refgffObj.symbols) + "))"
			self.__del_regex = re.compile("^(?:" + allowed_symbols + ":)?del:[0-9]+:[0-9]+$")
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

	# DATA IMPORT

	@staticmethod
	def hash(seq):
		'''
		Provides a SEGUID hash of a sequence (basically SHA-1 ot the sequence in upper-cases).
		'''
		return seguid(seq)

	def process_fasta(self, fname, cache=None, seqhashes_to_skip=None):
		'''
		Processes a given fasta and prepars the contained genome sequence for database import.
		Relevant sequence- and alignment-based information are provided in pickle format that
		is returned or, alternatively, cached in a file.

		Not using pickle cache:

		>>> os.chdir(os.path.dirname(os.path.realpath(__file__)))
		>>> db = sonarDB(DOCTESTDB + "_process_fasta_test")
		>>> data = db.process_fasta("doctest_b117.fna")
		>>> data[0]
		'b117'
		>>> data[1]
		'b117 Ideal severe acute respiratory syndrome coronavirus 2 lineage B.1.1.7, complete genome'
		>>> data[5]
		'C3267T C5388A T6954C del:11288:8 del:21765:5 del:21991:2 A23063T C23271A C23604A C23709T T24506G G24914C C27972T G28048T A28111G G28280C A28281T T28282A'
		>>> data[6]
		'N:D3L ORF8:Q27* ORF8:R52I S:del:68:2 ORF8:Y73C S:del:143:1 N:S235F S:N501Y S:A570D S:P681H S:T716I S:S982A ORF1ab:T1001I S:D1118H ORF1ab:A1708D ORF1ab:I2230T ORF1ab:del:3675:1'

		'''
		with sonarDBManager(self.db, readonly=True) as dbm:
			record = SeqIO.read(fname, "fasta")
			acc = record.id
			descr = record.description
			seq = str(record.seq.upper())
			seqhash = self.hash(seq)
			if not self.seq_exists(seqhash, dbm):
				alignment = sonarALIGN(seq, self.refseq, self.refgffObj)
				dnadiff = alignment.dnadiff
				aadiff = alignment.aadiff
				dna_profile = self.build_profile(*dnadiff)
				prot_profile = self.build_profile(*aadiff)
			else:
				alignment = None
				dnadiff = None
				aadiff = None
				dna_profile = None
				prot_profile = None

			if cache is None:
				return [acc, descr, seqhash, dnadiff, aadiff, dna_profile, prot_profile, seq]

			with open(cache, "wb") as handle:
				pickle.dump([acc, descr, seqhash, dnadiff, aadiff, dna_profile, prot_profile, seq], handle)

	def import_genome_from_fasta(self, *fnames, msg=None, paranoid=True):
		'''
		Imports genome from fasta file(s).

		>>> os.chdir(os.path.dirname(os.path.realpath(__file__)))
		>>> db = sonarDB(DOCTESTDB)
		>>> db.import_genome_from_fasta("doctest_b117.fna")

		'''
		if not msg is None:
			rng = tqdm(range(len(fnames)), desc = msg)
		else:
			rng = range(len(fnames))

		with sonarDBManager(self.db) as dbm:
			for i in rng:
				data = self.process_fasta(fnames[i])
				self.import_genome(*data, dbm=dbm, paranoid=paranoid)

	def import_genome_from_pickle(self, *fnames, msg=None, paranoid=True):
		'''
		Imports genome from cached pickle file(s).

		>>> os.chdir(os.path.dirname(os.path.realpath(__file__)))
		>>> db = sonarDB(DOCTESTDB)
		>>> db.import_genome_from_pickle("doctest_b117.pickle")

		'''
		if not msg is None:
			rng = tqdm(range(len(fnames)), desc = msg)
		else:
			rng = range(len(fnames))

		with sonarDBManager(self.db) as dbm:
			for i in rng:
				with open(fnames[i], 'rb') as handle:
					data = pickle.load(handle, encoding="bytes")
				self.import_genome(*data, dbm=dbm, paranoid=paranoid)

	def import_genome_from_cache(self, cacheObj, msg=None, paranoid=True):
		'''
		Imports genome from sonar cache object.

		>>> os.chdir(os.path.dirname(os.path.realpath(__file__)))
		>>> cache = sonarCache()
		>>> cachefile=cache.add_fasta("doctest_b117.fna")
		>>> db = sonarDB(DOCTESTDB + "_import_genome_from_cache")
		>>> db.import_genome_from_cache(cache)

		'''
		seqhashes = list(cacheObj.get_cached_seqhashes())
		if not msg is None:
			rng = tqdm(range(len(seqhashes)), desc = msg)
		else:
			rng = range(len(seqhashes))

		with sonarDBManager(self.db) as dbm:
			for i in rng:
				seqhash = seqhashes[i]
				processed_data = self.process_fasta(cacheObj.cached_fasta_name(seqhash))[2:]
				seq = cacheObj.get_cached_seq(seqhash)
				pfile = cacheObj.cached_pickle_name(seqhash)
				for acc, descr in cacheObj.get_acc_descr(seqhash):
					self.import_genome(acc, descr, *processed_data, dbm=dbm, paranoid=paranoid)

	def import_genome(self, acc, descr, seqhash, dnadiff, aadiff, dna_profile, prot_profile, seq, dbm, paranoid=True):
		'''
		imports processed genome data to database
		'''
		if not self.genome_exists(acc, descr, seqhash, dbm):
			self.insert_genome(acc, descr, seqhash, dbm)
		if not self.seq_exists(seqhash, dbm):
			self.insert_sequence(seqhash, dbm)
			self.insert_profile(seqhash, dna_profile, prot_profile, dbm)
			for ref, alt, s, e, _, __ in dnadiff:
				varid = self.get_dna_varid(ref, alt, s, e, dbm)
				if varid is None:
					varid = self.insert_dna_var(ref, alt, s, e, dbm)
					self.insert_sequence2dna(seqhash, varid, dbm)

			for ref, alt, s, e, protein, locus in aadiff:
				varid = self.get_prot_varid(protein, locus, ref, alt, s, e, dbm)
				if varid is None:
					varid = self.insert_prot_var(protein, locus, ref, alt, s, e, dbm)
				self.insert_sequence2prot(seqhash, varid, dbm)
		if True:
			self.be_paranoid(acc, seq, dbm, auto_delete=True)

	def insert_genome(self, acc, descr, seqhash, dbm):
		return dbm.insert('genome', ['accession', 'description', 'seqhash'], [acc, descr, seqhash])

	def insert_sequence(self, seqhash, dbm):
		return dbm.insert('sequence', ['seqhash'], [seqhash], ignore=True)

	def insert_profile(self, seqhash, dna_profile, aa_profile, dbm):
		dna_profile = " " + dna_profile.strip() + " "
		aa_profile = " " + aa_profile.strip() + " "
		return dbm.insert('profile', ['seqhash', 'dna_profile', 'aa_profile'], [seqhash, dna_profile, aa_profile], ignore=True)

	def insert_dna_var(self, ref, alt, start, end, dbm):
		return dbm.insert('dna', ['varid', 'ref', 'alt', 'start', 'end'], [None, ref, alt, start, end], ignore=True)

	def insert_sequence2dna(self, seqhash, varid, dbm):
		return dbm.insert('sequence2dna', ('seqhash', 'varid'), [seqhash, varid], ignore=True)

	def insert_prot_var(self, protein, locus, ref, alt, start, end, dbm):
		return dbm.insert('prot', ['varid', 'protein', 'locus', 'ref', 'alt', 'start', 'end'], [None, protein, locus, ref, alt, start, end], ignore=True)

	def insert_sequence2prot(self, seqhash, varid, dbm):
		return dbm.insert('sequence2prot', ['seqhash', 'varid'], [seqhash, varid], ignore=True)

	# DELETE DATA

	def delete_accession(self, acc, dbm):
		dbm.delete("genome", whereClause="accession = ?", whereVals=[acc], print_sql=False)

	# DATA SELECTS

	def genome_exists(self, acc, descr, seqhash, dbm):
		row = self.get_genome_data(acc, dbm=dbm)
		if row and row[0]['description'] == descr and row[0]['seqhash'] == seqhash:
			return True
		else:
			return False

	def seq_exists(self, seqhash, dbm):
		return dbm.select('sequence', fieldList=['COUNT(*)'], whereClause="seqhash = ?", whereVals=[seqhash], orderby=None, fetchone=True, print_sql=False)['COUNT(*)'] > 0

	def get_genome_data(self, *accs, dbm):
		where = " OR ".join(["accession = ?"] * len(accs))
		return dbm.select('genome', whereClause=where, whereVals=accs)

	def get_dna_varid(self, ref, alt, start, end, dbm):
		row = dbm.select('dna', whereClause="ref = ? AND alt = ? AND start = ? and end = ?", whereVals=[ref, alt, start, end], fetchone=True)
		if row:
			return row['varid']
		return None

	def get_prot_varid(self, protein, locus, ref, alt, start, end, dbm):
		row = dbm.select('prot', whereClause="protein = ? AND locus = ? AND ref = ? AND alt = ? AND start = ? AND end = ?", whereVals=[protein, locus, ref, alt, start, end], fetchone=True)
		if row:
			return row['varid']
		return None

	def get_dna_vars(self, acc, dbm):
		return dbm.select('dna_view', whereClause="accession = ?", whereVals=[acc], orderby="start DESC")
	def iter_table(self, table):
		sql = dbm.select('dna', whereClause="ref = ? AND alt = ? AND start = ? and end = ?", whereVals=[ref, alt, start, end], fetchone=True)
		for row in dbm.select(table):
			yield row

	# NOMENCLATURE

	def isdnavar(self, var):
		'''
		Returns True if a variant definition is a valid dna variant expression.

		>>> db = sonarDB(DOCTESTDB)
		>>> db.isdnavar("S:N501Y")
		False
		>>> db.isdnavar("A101T")
		True

		'''
		return bool(self.dna_var_regex.match(var))

	def isaavar(self, var):
		'''
		Returns True if a variant definition is a valid aa variant expression.

		>>> db = sonarDB(DOCTESTDB)
		>>> db.isaavar("S:N501Y")
		True
		>>> db.isaavar("A101T")
		False

		'''
		return bool(self.aa_var_regex.match(var))

	def isdel(self, var):
		'''
		Returns True if a variant definition is a valid aa or dna deletion expression.

		>>> db = sonarDB(DOCTESTDB)
		>>> db.isdel("del:100-118")
		False
		>>> db.isdel("del:100:18")
		True
		>>> db.isdel("ORF1ab:del:5:2")
		True

		'''
		return bool(self.del_regex.match(var))

	# PROFILE BUILDING

	def build_profile(self, *vars):
		'''
		build a profile based on given mutations according to sonarDB's profile syntax
		'''
		if len(vars) == 0:
			return ""
		profile = []
		if len(vars) == 1:
			this_ref, this_alt, this_start, this_end, this_protein, this_locus = vars[0]
		else:
			vars = sorted(vars, key=lambda x: x[2])
			for l in range(len(vars)-1):
				this_ref, this_alt, this_start, this_end, this_protein, this_locus = vars[l]
				next_ref, next_alt, next_start, next_end, next_protein, next_locus = vars[l+1]
				if this_alt != "":
					var = self.format_var(this_ref, this_alt, this_start, this_end, this_protein)
					if var not in profile:
						profile.append(var)
				elif this_alt == "" and this_start + len(this_ref) == next_start and this_protein == next_protein and this_locus == next_locus:
					vars[l+1] = (this_ref + next_ref, vars[l+1][1], this_start, next_start, this_protein, this_locus)
				else:
					var = self.format_var(this_ref, this_alt, this_start, this_end, this_protein)
					if var not in profile:
						profile.append(var)
		var = self.format_var(this_ref, this_alt, this_start, this_end, this_protein)
		if var not in profile:
			profile.append(var)

		return " ".join(profile)

	@staticmethod
	def format_var(ref, alt, start, end, protein=None, locus=None):
		'''
		format single mutation based on according to sonarDB's profile syntax
		'''
		if end is None:
			coord = str(start+1)
		else:
			ref = "del:"
			coord = str(start+1) + ":" + str(end-start)
		protein = protein + ":" if protein else ""
		return protein + ref + coord + alt

	# MATCHING

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


	def pinpoint_mutation(self, mutation, code):
		'''
		Returns a set of explicit mutations based on the given possibly ambiguous mutation definition.
		The mutation definition must follow the covsonar nomenclature system.

		>>>
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

	def match_builder(self, profile, dna=True, exclusive=False, negate=False):
		'''
		build where clause for profile matching.
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

	def match(self, include_profiles=None, exclude_profiles=None, accessions=None, lineages=None, zips=None, dates=None, exclusive=False, ambig=False, show_sql=False):
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
				includes.append(self.match_builder(profile, exclusive=exclusive, negate=False))
			if len(includes) > 1:
				includes = [ "(" + x + ")" if " AND " in x else x for x in includes ]
			clause.append(" OR ".join(includes))
			if len(includes) > 1:
				clause[-1] = "(" + clause[-1] + ")"

		# adding conditions of profiles to exclude to where clause
		if exclude_profiles:
			excludes = []
			for profile in exclude_profiles:
				excludes.append(self.match_builder(profile, exclusive=exclusive, negate=True))
			excludes = [ "(" + x + ")" if " AND " in x else x for x in excludes ]
			clause.append(" AND ".join(excludes))

		# adding accession, lineage, zips, and dates based conditions
		if accessions:
			clause.append("accession IN (" + " , ".join(['?'] * len(accessions)) + ")")
			vals.extend(accessions)
		if lineages:
			clause.append("lineages IN (" + " ,".join(['?'] * len(lineages))  + ")")
			vals.extend(lineages)
		if zips:
			z = []
			for zp in zips:
				z.append("zip LIKE '" + str(zp) + "%'")
			clause.append(" OR ".join(z))
			if len(z) > 1:
				clause[-1] = "(" + clause[-1] + ")"
		if dates:
			d = []
			for dt in dates:
				if ":" in dt:
					x, y = dt.split(":")
					d.append("(date BETWEEN '" + x + "' AND '" + y + "')")
				else:
					d.append("date = " + dt)
			clause.append(" OR ".join(d))
			if len(d) > 1:
				clause[-1] = "(" + clause[-1] + ")"

		# executing the query and storing results
		clause = " AND ".join(clause)
		with sonarDBManager(self.db, readonly=True) as dbm:
			rows = [x for x in dbm.select("essence", whereClause=clause, whereVals=vals)]

		# remove ambiguities from database profiles if wished
		if not ambig:
			keep = [item for sublist in include_profiles for item in sublist] if include_profiles else None
			for i in range(len(rows)):
				rows[i]['dna_profile'] = self.filter_ambig(rows[i]['dna_profile'], self.iupac_explicit_nt_code, keep)
				rows[i]['aa_profile'] = self.filter_ambig(rows[i]['aa_profile'], self.iupac_explicit_aa_code, keep)

		return rows

	# VALIDATION

	def restore(self, acc, dbm, aligned=False, seq_return=False):
		rows = self.get_dna_vars(acc, dbm)
		if rows:
			gap = "-" if aligned else ""
			refseq = list(self.refseq)
			qryseq = refseq[:]
			for row in rows:
				if not row['start'] is None:
					s = row['start']
					if row['ref'] != refseq[s]:
						sys.exit("data error: data inconsistency found for '" + acc + "' (" + row['ref']+ " expected at position " + str(s+1) + " of the reference sequence, got " + refseq[s] + ").")
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

	def be_paranoid(self, acc, orig_seq, dbm, auto_delete=True):
		orig_seq = orig_seq.upper()
		restored_seq = self.restore(acc, aligned=False, seq_return=True, dbm=dbm)
		if orig_seq != restored_seq:
			if auto_delete:
				self.delete_accession(acc, dbm)
			sys.exit("Good that you are paranoid: " + acc + " original and those restored from the database do not match.")

	def show_seq_diffs(self, seq1, seq2, stderr=False):
		target = sys.stderr if stderr else None
		for i,s in enumerate(difflib.ndiff(seq2, seq1)):
			if s[0]==' ': continue
			elif s[0]=='-':
				print(u'wrong "{}" at position {}'.format(s[-1],i), file = target)
			elif s[0]=='+':
				print(u'missing "{}" at position {}'.format(s[-1],i), file = target)

	# OTHER

	@staticmethod
	def get_version():
		with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), ".version"), "r") as handle:
			return handle.read().strip()

class sonarCache():
	def __init__(self, dir=None):
		'''
		following parameters are needed to initialize the class:
		dbfile = database file (will be created if it does not exists)
		timeout = busy time out of the database
		readonly = allow read only access only
		'''
		self.temp = not bool(dir)
		if self.temp:
			self.dirname = mkdtemp(prefix=".sonarCache_")
		else:
			self.dirname = os.path.abspath(dir)
			self.checkdir()
		self.cache = {}

	def __enter__(self):
		pass

	def __exit__(self, exc_type, exc_value, exc_traceback):
		if [exc_type, exc_value, exc_traceback].count(None) != 3:
			print("warning:", file=sys.stderr)
			print(traceback.format_exc(), file=sys.stderr)
		if os.path.isdir(self.dirname) and self.temp:
			shutil.rmtree(self.dirname)

	def checkdir(self):
		if not os.path.isdir(self.dirname):
			os.makedirs(self.dirname)
		elif os.listdir(self.dirname):
			self.restore_cache()

	def restore_cache(self):
		acclog = self.cached_acclog_name()
		if os.listdir(self.dirname):
			if not os.path.isfile(acclog):
				sys.exit("cache error: not a valid cache directory.")
			self.cache = self.load_pickle(acclog)
			for seqhash in self.get_cached_seqhashes():
				if not os.path.isfile(self.cached_fasta_name(seqhash)) or not os.path.isfile(self.cached_pickle_name(seqhash)) :
					sys.exit("cache error: cache directory is corrupted.")

	@staticmethod
	def slugify(string):
		'''
		Provides a filesystem safe representation of a string.
		'''
		return base64.urlsafe_b64encode(string.encode('UTF-8') ).decode('UTF-8')

	def iter_fasta(self, fname):
		'''
		generates an iterator returning the entries of a given FASTA file as tuple of header
		and sequence.
		'''
		for record in SeqIO.read(fname, "fasta"):
			yield record.description, str(record.seq).upper()

	def read_fasta(self, fname):
		'''
		reads single entry fasta file and returns a tuple of header and sequence.
		Exception raises if multiple entries are in the respective fasta file.
		'''
		record = SeqIO.read(fname, "fasta")
		return record.description[1:], str(record.seq).upper()

	def cached_fasta_name(self, seqhash):
		return os.path.join(self.dirname, self.slugify(seqhash) + ".fasta")

	def cached_pickle_name(self, seqhash):
		return os.path.join(self.dirname, self.slugify(seqhash) + ".pickle")

	def link_pickle_name(self, fastaname):
		with open(fastaname) as handle:
			seqhash = handle.readline().strip()[1:]
		return self.cached_pickle_name(seqhash)

	def cached_acclog_name(self):
		return os.path.join(self.dirname, "acclist")

	@staticmethod
	def load_pickle(fname):
		with open(fname, 'rb') as handle:
			return pickle.load(handle, encoding="bytes")

	def write_pickle(self, fname, data):
		with open(fname, 'wb') as handle:
			pickle.dump(data, handle)


	def add_fasta(self, fname):
		'''
		Adds entries of a given fasta file in a sequence-unredundant manner to
		the cache.
		'''
		fnames = []
		for record in SeqIO.parse(fname, "fasta"):
			acc = record.id
			descr = record.description
			seq = str(record.seq).upper()
			seqhash = sonarDB.hash(seq)
			fnames.append(self.cached_fasta_name(seqhash))

			# check for seqhash collision
			if os.path.isfile(fnames[-1]):
				_, s = self.read_fasta(fnames[-1])
				if s != seq:
					sys.exit("cache error: sequence collision for hash '" + seqhash + "'")
			else:
				with open(fnames[-1], "w") as handle:
					handle.write(">" + seqhash + os.linesep + seq)

			# check for accession collision
			if acc in self.cache and self.cache[acc] != (descr, seqhash):
				sys.exit("cache error: sequence collision for accession '" + acc + "'")

			if not acc in self.cache:
				self.cache[acc] = (descr, seqhash)

		#self.write_pickle(self.cached_acclog_name(), self.cache[acc])
		return fnames

	def get_cached_seqhashes(self):
		return sorted(set([x[1] for x in self.cache.values()]))

	def get_cached_fasta_files(self):
		return [ self.cached_fasta_name(x) for x in self.get_cached_seqhashes ]

	def get_cached_seq(self, seqhash):
		return self.read_fasta(self.cached_fasta_name(seqhash))[1]

	def get_acc_descr(self, seqhash):
		return [(x, self.cache[x][0]) for x in self.cache if self.cache[x][1] == seqhash]

if __name__ == "__main__":
	import doctest
	global DOCTESTDIR, DOCTESTDB
	print("sonarDB", sonarDB.get_version())
	print("performing unit tests ...")
	with TemporaryDirectory() as tmpdirname:
		DOCTESTDIR = tmpdirname
		DOCTESTDB = os.path.join(DOCTESTDIR, "testdb")
		print(doctest.testmod(verbose=False))
