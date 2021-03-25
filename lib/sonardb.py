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
from urllib.parse import quote
from math import floor

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
	def __init__(self, query, target, sonarGFFObj = None):
		'''
		following parameters are needed to initialize the class:
		query = query nucleotide sequence to be algned to the target sequence
		target = target nucleotide sequence
		sonarGFFObj = sonarGFF object to consider cds annotation and provide protein level based functions
		'''
		self.aligned_query, self.aligned_target = self.align_dna(query.upper(), target.upper())
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

	def align_dna(self, query, target, matrixfile = os.path.join(os.path.dirname(os.path.realpath(__file__)), "EDNAFULL")):
		'''
		align two nucleotide sequences and return a tuple of the aligned sequences.
		EDNAFULL is used as default matrix.

		>>> algn = sonarALIGN("ATTGTTGTTATAATGGCCAGTT", "ATTGATGTTATGAATGGCCTTT")
		>>> algn.align_dna("ATTGTTGTTATAATGGCCAGTT", "ATTGATGTTATGAATGGCCTTT")
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
		aligner.open_gap_score = -10
		aligner.extend_gap_score = -0.5
		aligner.end_extend_gap_score = -10
		aligner.end_gap_score = -0.5
		aligner.target_end_gap_score = 0.0
		aligner.query_end_gap_score = 0.0
		alignment = str(sorted(aligner.align(query, target), key=lambda x: x.score, reverse=True)[0]).split("\n")
		return alignment[0], alignment[2]

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

		>>> algn = sonarALIGN("ATTGTTGTTATAATGGCCAGTT", "ATTGATGTTATGAATGGCCTTT")
		>>> algn.real_pos(22)
		20
		'''
		return x - self.aligned_target[:x+1].count("-")

	def align_pos(self, x):
		'''
		converts a target sequence position to the respective position in the alignment.
		Positions are meant to be 0-based.

		>>> algn = sonarALIGN("ATTGTTGTTATAATGGCCAGTT", "ATTGATGTTATGAATGGCCTTT")
		>>> algn.align_pos(20)
		22
		'''
		return sum(self.target_coords_matrix[:x])

	def iter_dna_vars(self):
		'''
		generating an iterator returning tuples for each variant alignment position
		consisting of target nucleotide, query nucleotide(s), target start position, target end position, None, None
		Last two None elements exist to harmonize output syntax between this and the iter_aa_vars function.

		>>> algn = sonarALIGN("ATTGTTGTTATATGGCCAGTT", "ATTGATGTTATGAATGGCCTTT")
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
	def __init__(self, dbfile, timeout=-1, readonly=True):
		self.dbfile = os.path.abspath(dbfile)
		self.connection = None
		self.cursor = None
		self.__timeout = timeout
		self.__mode = "ro" if readonly else "rwc"
		self.dburi = "file:" + urlquote(self.dbfile) + "?mode=" + self.__mode

	def __enter__(self):
		if not os.path.isfile(self.dbfile):
			self.create_scheme()
		self.connection = self.connect()
		self.connection.row_factory = self.dict_factory
		self.cursor = self.cursor()
		if self.__mode != "ro":
			self.cursor.execute("BEGIN TRANSACTION")
		return self

	def __exit__(self, exc_type, exc_value, exc_traceback):
		if [exc_type, exc_value, exc_traceback].count(None) != 3:
			print("warning:", file=sys.stderr)
			print(e, file=sys.stderr)
			if self.__mode != "ro":
				print("rollback", file=sys.stderr)
				self.rollback()
		else:
			self.cursor.execute("END TRANSACTION")
			self.connection.commit()
		self.close()

	def connect(self):
		return sqlite3.connect(self.dbfile, self.__timeout, isolation_level = None, uri = self.__uri)

	def commit(self):
		self.connection.commit()

	def rollback(self):
		self.connection.rollback()

	def close(self):
		self.connection.close()

	def create_scheme(self):
		with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), "db.sqlite"), 'r') as handle:
			sql = handle.read()
		con.sqlite3.connect(self.uri, uri=True)
		con.executescript(sql)
		con.close()

	def select(self, table, fieldList=['*'], whereClause=None, whereVals=None, orderby=None, print_sql=False):
		sql = "SELECT " + "".join(fieldList) + " FROM " + table
		if whereClause:
			sql += " WHERE " + whereClause
		if orderby:
			sql += " ORDER BY " + orderby
		if print_sql:
			print(sql)
		return self.cursor.execute(sql, valList).fetchall()

	def delete(self, table, fieldList=['*'], whereClause=None, whereVals=None, print_sql=False):
		sql = "DELETE FROM " + table
		if whereClause:
			sql += " WHERE " + whereClause
		if print_sql:
			print(sql)
		return self.cursor.execute(sql, valList).fetchall()

	def insert(self, table, fieldList, valList, ignore=True, vals=None, print_sql=False):
		ignore = "" if not ignore else "OR IGNORE "
		sql = "INSERT " + ignore + "INTO " + table + "(" + ", ".join(fieldList) + ") VALUES (" + ", ".join(['?']*len(valList)) + ")"
		return self.cursor.execute(sql, valList).rowid

	@staticmethod
	def dict_factory(cursor, row):
		d = OrderedDict()
		for idx, col in enumerate(cursor.description):
			d[col[0]] = row[idx]
		return d


class sonarDB(object):
	def __init__(self, dbfile, translation_table = 1):
		self.moduledir = os.path.dirname(os.path.realpath(__file__))
		self.reffna = os.path.join(self.moduledir, "ref.fna")
		self.refgff = os.path.join(self.moduledir, "ref.gff3")
		self.translation_table = translation_table
		self.__refseq = None
		self.__refdescr = None
		self.__refgffObj = None

	@property
	def refseq(self):
		if not self.__refseq:
			with open(self.reffna, "r") as handle:
				self.__refseq = "".join([x.strip().upper() for x in handle.readlines()[1:]])
		return self.__refseq

	@property
	def refdescr(self):
		if not self.__refdescr:
			with open(self.reffna, "r") as handle:
				self.__refdescr = handle.readline().strip()[1:]
		return self.__refdescr

	@property
	def refgffObj(self):
		if not self.__refgffObj:
			self.__refgffObj = sonarGFF(self.refgff, self.refseq, self.translation_table)
		return self.__refgffObj

	def slugify_fname(self, fname):
		return os.path.join(self.seqdir, base64.urlsafe_b64encode(fname.encode('UTF-8') ).decode('UTF-8') + ".fna")

	@staticmethod
	def hash(seq):
		return seguid(seq)

	def process_fasta(self, fname, cache=None):
		record = SeqIO.read(fname, "fasta")
		acc = record.id
		descr = record.description
		seq = str(record.seq.upper())
		seqhash = self.hash(seq)
		seq_exists = self.sequence_exists(seqhash)
		if not seq_exists:
			alignment = sonarALIGN(seq, self.refseq, self.refgffObj)
			dnadiff = alignment.dnadiff
			aadiff = alignment.aadiff
			dna_profile = self.create_profile(*dnadiff)
			prot_profile = self.create_profile(*aadiff)
		else:
			alignment = None
			dnadiff = None
			aadiff = None
			dna_profile = None
			prot_profile = None

		if cache is None:
			return [acc, descr, seqhash, dnadiff, aadiff, dna_profile, prot_profile]

		with open(cache, "wb") as handle:
			pickle.dump([acc, descr, seqhash, dnadiff, aadiff, dna_profile, prot_profile], handle)

	def import_genome_from_fasta(self, *fnames, msg=None):
		if not msg is None:
			rng = tqdm(range(len(fnames)), desc = msg)
		else:
			rng = range(len(fnames))
		for i in rng:
			data = self.process_fasta(fnames[i])
			with sonarDBManager(self.dbfile) as dbm:
				self.import_genome(*data, dbcursor=dbm.cursor)

	def import_genome_from_pickle(self, *fnames, msg=None):
		if not msg is None:
			rng = tqdm(range(len(fnames)), desc = msg)
		else:
			rng = range(len(fnames))

		with sonarDBManager(self.dbfile) as dbm:
			for i in rng:
				with open(fnames[i], 'rb') as handle:
					self.import_genome(*pickle.load(handle, encoding="bytes"), dbcursor=dbm.cursor)

	def import_genome(self, acc, descr, seqhash, dnadiff, aadiff, dna_profile, prot_profile, dbcursor):
		if not self.genome_exists(acc, descr, seqhash, dbcursor):
			self.insert_genome(acc, descr, seqhash, dbcursor)
		if not self.sequence_exists(seqhash, dbcursor):
			self.insert_sequence(seqhash, dbcursor)
			self.insert_profile(seqhash, dna_profile, prot_profile, dbcursor)
			for ref, alt, s, e, _, __ in dnadiff:
				varid = self.get_dna_varid(ref, alt, s, e)
				if varid is None:
					varid = self.insert_dna_var(ref, alt, s, e, dbcursor)
					self.insert_sequence2dna(seqhash, varid, dbcursor)

			for ref, alt, s, e, protein, locus in aadiff:
				varid = self.get_prot_varid(protein, locus, ref, alt, s, e, dbcursor)
				if varid is None:
					varid = self.insert_prot_var(protein, locus, ref, alt, s, e, dbcursor)
				self.insert_sequence2prot(seqhash, varid, dbcursor)

	def create_profile(self, *vars):
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
					profile.append(self.format_var(this_ref, this_alt, this_start, this_end, this_protein))
				elif this_alt == "" and this_start + len(this_ref) == next_start and this_protein == next_protein and this_locus == next_locus:
					vars[l+1] = (this_ref + next_ref, vars[l+1][1], this_start, next_start, this_protein, this_locus)
				else:
					profile.append(self.format_var(this_ref, this_alt, this_start, this_end, this_protein, this_locus))
		profile.append(self.format_var(this_ref, this_alt, this_start, this_end, this_protein, this_locus))
		return " ".join(profile)

	@staticmethod
	def format_var(ref, alt, start, end, protein=None, locus=None):
		if end is None:
			coord = str(start+1)
		else:
			ref = "del:"
			coord = str(start+1) + ":" + str(end-start)
		protein = protein + ":" if protein else ""
		return protein + ref + coord + alt

	def delete_accession(self, *accs):
		with sqlite3.connect(self.dbfile) as con:
			sql = "DELETE FROM genome WHERE accession in (" + ", ".join(["?"] * len(accs)) + ")"
			con.execute(sql, accs)
			con.commit()


	def get_dna_profiles(self, *acc):
		sql = "SELECT accession, start, end, ref, alt from dna_view  WHERE " + " OR ".join(["accession = ?" for x in acc])
		rows = self.rocon.execute(sql, acc).fetchall()
		return rows

	def insert_genome(self, acc, descr, seqhash, dbcursor):
		return dbcursor.execute('INSERT INTO genome(accession, description, seqhash) VALUES (?, ?, ?)', (acc, descr, seqhash)).lastrowid

	def insert_sequence(self, seqhash, dbcursor):
		return dbcursor.execute('INSERT OR IGNORE INTO sequence(seqhash) VALUES(?)', (seqhash, )).lastrowid

	def insert_profile(self, seqhash, dna_profile, aa_profile, dbcursor):
		dna_profile = " " + dna_profile.strip() + dna_profile
		aa_profile = " " + aa_profile.strip() + aa_profile
		return dbcursor.execute('INSERT OR IGNORE INTO profile(seqhash, dna_profile, aa_profile) VALUES(?, ?, ?)', (seqhash, dna_profile, aa_profile)).lastrowid

	def insert_dna_var(self, ref, alt, start, end, dbcursor):
		return dbcursor.execute('INSERT INTO dna(varid, ref, alt, start, end) VALUES(?, ?, ?, ?, ?)', (None, ref, alt, start, end)).lastrowid

	def insert_sequence2dna(self, seqhash, varid, dbcursor):
		dbcursor.execute('INSERT INTO sequence2dna(seqhash, varid) VALUES(?, ?)', (seqhash, varid))

	def insert_prot_var(self, protein, locus, ref, alt, start, end, dbcursor):
		return dbcursor.execute('INSERT INTO prot(varid, protein, locus, ref, alt, start, end) VALUES(?, ?, ?, ?, ?, ?, ?)', (None, protein, locus, ref, alt, start, end)).lastrowid

	def insert_sequence2prot(self, seqhash, varid,dbcursor):
		return dbcursor.execute('INSERT INTO sequence2prot(seqhash, varid) VALUES(?, ?)', (seqhash, varid)).lastrowid

	def get_genome_data(self, *acc, dbcursor):
		sql = "SELECT * from genome WHERE " + " OR ".join(["accession = ?" for x in acc])
		return dbcursor.execute(sql, acc).fetchall()

	def genome_exists(self, acc, descr, seqhash, dbcursor):
		row = self.get_genome_data(acc, dbcursor=dbcursor)
		if row and row[0]['description'] == descr and row[0]['seqhash'] == seqhash:
			return True
		else:
			return False

	def sequence_exists(self, seqhash, dbcursor):
		sql = "SELECT COUNT(*) from sequence WHERE seqhash = ?"
		return dbcursor.execute(sql, (seqhash, )).fetchone()['COUNT(*)'] > 0

	def get_dna_varid(self, ref, alt, start, end):
		sql = "SELECT varid from dna WHERE ref = ? AND alt = ? AND start = ? and end = ?"
		row = self.rocon.execute(sql, (ref, alt, start, end)).fetchone()
		if row:
			return row['varid']
		return None

	def get_prot_varid(self, protein, locus, ref, alt, start, end, dbcursor):
		sql = "SELECT varid from prot WHERE protein = ? AND locus = ? AND ref = ? AND alt = ? AND start = ? and end = ?"
		row = dbcursor.execute(sql, (protein, locus, ref, alt, start, end)).fetchone()
		if row:
			return row['varid']
		return None

	def show_profiles(self, *acc, branch="both"):
		if branch == "dna":
			profile = "dna_profile"
		elif branch == "prot":
			profile = "aa_profile"
		elif branch == "both":
			profile = "dna_profile, aa_profile"
		sql = "SELECT accession, " + profile + " from essence WHERE " + " OR ".join(["accession = ?" for x in acc])
		for row in self.rocon.execute(sql, acc):
			yield row

	def iter_rows(self, table):
		sql = "SELECT * FROM " + table
		for row in self.rocon.execute(sql):
			yield row

	def create_tables(self):
		with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "resources", "db.sqlite"), 'r') as handle:
			sql = handle.read()
		with sqlite3.connect(self.dbfile) as con:
			con.executescript(sql)

	@staticmethod
	def get_version():
		with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), ".version"), "r") as handle:
			return handle.read().strip()

if __name__ == "__main__":
	import doctest
	print("sonarDB", sonarDB.get_version())
	doctest.testmod()
