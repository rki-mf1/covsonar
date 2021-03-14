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

class sonarCDS(object):
	def __init__(self, symbol, start, end, strand, seq, locus=None, translation_table=1):
		self.symbol = symbol
		self.locus = locus
		self.start = start
		self.end = end
		self.strand = "+" if strand == "+" or int(strand) == 1 else "-"
		self.dna = seq
		self.translation_table = translation_table
		self.__aa = None

	@property
	def aa(self):
		if self.__aa is None:
			self.__aa = str(Seq.translate(self.dna, table=self.translation_table, to_stop=True))
		return self.__aa

	@property
	def coords(self):
		return self.start, self.end

	@property
	def range(self):
		return range(self.start, self.end)

class sonarGFF(object):
	def __init__(self, gff3, genomeseq, translation_table=1):
		self.translation_table = translation_table
		self.cds = self.process_gff3(gff3, genomeseq)

	def process_gff3(self, gff_fname, genomeseq):
		symbol_regex = re.compile("gene=([^;]+)(?:;|$)")
		locus_regex = re.compile("locus_tag=([^;]+)(?:;|$)")

		with open(gff_fname, "r") as handle:
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
					dna = genomeseq[s:e]
					if strand == "-":
						dna = str(Seq.reverse_complement(dna))
					cds.add(sonarCDS(symbol, s, e, strand, dna, locus, self.translation_table))
		return cds

class sonarALIGN(object):
	def __init__(self, query, target, sonarGFFObj = None):
		self.aligned_query, self.aligned_target = self.align_dna(query, target)
		self.gff = sonarGFFObj if sonarGFFObj else None
		self.indel_regex = re.compile(".-+")
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

	def align_dna(self, query, target):
		aligner = Align.PairwiseAligner()
		aligner.mode = 'global'
		aligner.open_gap_score = -10
		aligner.extend_gap_score = -0.5
		aligner.target_end_gap_score = 0.0
		aligner.query_end_gap_score = 0.0
		alignment = str(sorted(aligner.align(query, target), key=lambda x: x.score, reverse=True)[0]).split("\n")
		return alignment[0], alignment[2]

	def align_aa(self, query, target):
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
		return x - self.aligned_target[:x+1].count("-")

	def align_pos(self, x):
		return sum(self.target_coords_matrix[:x])

	def compare_aligned_seqs(self, query, target):
		target = target.upper()
		query = query.upper()

		# insertions
		for match in self.indel_regex.finditer(target):
			s = match.start() - target[:match.start()].count("-")
			yield match.group()[0], query[match.start():match.end()], s, None

		# deletions
		for match in self.indel_regex.finditer(query):
			s = match.start() - query[:match.start()].count("-")
			e = s + len(match.group())
			yield target[match.start():match.end()], match.group()[0], s, e

		# sn/aps
		i = -1
		for pair in zip(target, query):
			if pair[1] == "-":
				i += 1
			elif pair[0] == pair[1]:
				i += 1
			elif pair[0] != "-":
				i += 1
				yield pair[0], pair[1], i, None

	def iter_dna_vars(self):
		for t, q, s, e in self.compare_aligned_seqs(self.aligned_query, self.aligned_target):
			s = self.real_pos(s)
			if not e is None:
				e = self.real_pos(e)
			yield t, q, s, e

	def iter_aa_vars(self):
		if self.gff:
			for cds in self.gff.cds:
				if cds.strand == "+":
					x = self.align_pos(cds.start)
					qdna = self.aligned_query[x:].replace("-", "")
				else:
					x = self.align_pos(cds.end)
					qdna = str(Seq.reverse_complement(self.aligned_query[:x].replace("-", "")))
				qaa = str(Seq.translate(qdna, table=self.gff.translation_table, to_stop=True))

				aligned_query_aa, aligned_target_aa = self.align_aa(qaa, cds.aa)

				for t, q, s, e in self.compare_aligned_seqs(aligned_query_aa, aligned_target_aa):
					yield cds.symbol, cds.locus, t, q, s, e


class sonarDB(object):
	def __init__(self, db, translation_table = 1, check_db=True):
		self.dbdir = os.path.abspath(db)
		self.vfile = os.path.join(self.dbdir, ".version")
		self.dbfile = os.path.join(self.dbdir, "sonar.db")
		self.datadir = os.path.join(self.dbdir, "data")
		self.seqdir = os.path.join(self.datadir, "seqs")
		self.algndir = os.path.join(self.datadir, "algn")
		self.lindir = os.path.join(self.datadir, "lin")
		self.pangodir = os.path.join(self.dbdir, "pangodir")
		self.pangoenv = os.path.join(self.dbdir, "pangoenv")
		self.modulebase = os.path.dirname(os.path.realpath(__file__))
		self.reffna = os.path.join(self.modulebase, "ref.fna")
		self.refgff = os.path.join(self.modulebase, "ref.gff3")
		self.translation_table = translation_table
		self.__max_supported_prev_version = "0.0.9"
		self.__refseq = None
		self.__iupac_nuc_code = None
		self.__annotation = None
		self.__rocon = None
		self.__refgffObj = None

		if check_db:
			self.check_dir()

	def __exit__(self, type, value, traceback):
		if self.__conn:
			self.__conn.close()

	@property
	def rocon(self):
		if self.__rocon is None:
			self.__rocon = sqlite3.connect(self.dbfile, uri=True)
			self.__rocon.row_factory = sonarDB.dict_factory
		return self.__rocon

	@property
	def refseq(self):
		if not self.__refseq:
			with open(self.reffna, "r") as handle:
				self.__refseq = "".join([x.strip().upper() for x in handle.readlines()[1:]])
		return self.__refseq

	@property
	def refgffObj(self):
		if not self.__refgffObj:
			self.__refgffObj = sonarGFF(self.refgff, self.refseq, self.translation_table)
		return self.__refgffObj

	@property
	def iupac_nuc_code(self):
		if not self.__iupac_nuc_code:
			self.__iupac_nuc_code = {
						"A": set("A"),
						"C": set("C"),
						"G": set("G"),
						"T": set("T"),
						"R": set("AG"),
						"Y": set("CT"),
						"S": set("GC"),
						"W": set("AT"),
						"K": set("GT"),
						"M": set("AC"),
						"B": set("CGT"),
						"D": set("AGT"),
						"H": set("ACT"),
						"V": set("ACG"),
						"N": set("ATGC")
					}
		return self.__iupac_nuc_code

	@staticmethod
	def hash(seq):
		return seguid(seq)

	# directory and file management

	def check_dir(self):
		if not os.path.isdir(self.dbdir) or len(os.listdir(self.dbdir)) == 0:
			self.create_dirs()
		elif not os.path.isdir(self.seqdir) or not os.path.isfile(self.vfile) or not os.path.isfile(self.dbfile):
			sys.exit("error: invalid sonar database directory")
		else:
			self.check_compatibility()

	def check_compatibility(self):
		with open(self.vfile, "r") as handle:
			dir_version = handle.read().strip()
		module_version = sonarDB.get_version()
		if version.parse(dir_version) < version.parse(self.__max_supported_prev_version):
			sys.exit("version conflict: sonar database version " + dir_version + " is not supported anymore by this sonar version (module_version).")
		if version.parse(dir_version) > version.parse(module_version):
			sys.exit("version conflict: please update sonar (current version: " + module_version + ") to fit your database (version: " + dir_version + ").")

	def create_dirs(self):
		try:
			# dirs
			os.makedirs(self.seqdir, exist_ok=True)

			# db
			self.create_tables()

			# version file
			with open(self.vfile, "w") as handle:
				handle.write(sonarDB.get_version())

		except Error as e:
			shutil.rmtree(self.dbdir)
			exit(e)

	def slugify_fname(self, fname):
		return os.path.join(self.seqdir, base64.urlsafe_b64encode(fname.encode('UTF-8') ).decode('UTF-8') + ".fna")

	# db management

	@staticmethod
	def dict_factory(cursor, row):
	    d = OrderedDict()
	    for idx, col in enumerate(cursor.description):
	        d[col[0]] = row[idx]
	    return d

	def add_genome_from_fasta(self, fname):
		record = SeqIO.read(fname, "fasta")
		acc = record.id
		descr = record.description
		seq = str(record.seq.upper())
		seqhash = sonarDB.hash(seq)

		if not self.sequence_exists(seqhash):
			alignment = sonarALIGN(seq, self.refseq, self.refgffObj)
		else:
			alignment = None

		with sqlite3.connect(self.dbfile) as con:
			if not self.genome_exists(acc, descr, seqhash):
				self.insert_genome(acc, descr, seqhash, con)
			if alignment:
				self.insert_sequence(seqhash, con)
				for ref, alt, s, e in alignment.iter_dna_vars():
					varid = self.get_dna_varid(ref, alt, s, e)
					if varid is None:
						varid = self.insert_dna_var(ref, alt, s, e, con)
					self.insert_sequence2dna(seqhash, varid, con)

				for protein, locus, ref, alt, s, e in alignment.iter_aa_vars():
					varid = self.get_prot_varid(protein, locus, ref, alt, s, e)
					if varid is None:
						varid = self.insert_prot_var(protein, locus, ref, alt, s, e, con)
					self.insert_sequence2prot(seqhash, varid, con)

	def insert_genome(self, acc, descr, seqhash, con):
		cursor = con.execute('INSERT INTO genome(accession, description, seqhash) VALUES (?, ?, ?)', (acc, descr, seqhash))
		return cursor.lastrowid

	def insert_sequence(self, seqhash, con):
		cursor = con.execute('INSERT INTO sequence(seqhash) VALUES(?)', (seqhash, ))
		return cursor.lastrowid

	def insert_dna_var(self, ref, alt, start, end, con):
		cursor = con.execute('INSERT INTO dna(varid, ref, alt, start, end) VALUES(?, ?, ?, ?, ?)', (None, ref, alt, start, end))
		return cursor.lastrowid

	def insert_sequence2dna(self, seqhash, varid, con):
		con.execute('INSERT INTO sequence2dna(seqhash, varid) VALUES(?, ?)', (seqhash, varid))

	def insert_prot_var(self, protein, locus, ref, alt, start, end, con):
		cursor = con.execute('INSERT INTO prot(varid, protein, locus, ref, alt, start, end) VALUES(?, ?, ?, ?, ?, ?, ?)', (None, protein, locus, ref, alt, start, end))
		return cursor.lastrowid

	def insert_sequence2prot(self, seqhash, varid, con):
		con.execute('INSERT INTO sequence2prot(seqhash, varid) VALUES(?, ?)', (seqhash, varid))

	def get_genome_data(self, *acc):
		sql = "SELECT * from genome WHERE " + " OR ".join(["accession = ?" for x in acc])
		return self.rocon.execute(sql, acc).fetchall()

	def genome_exists(self, acc, descr, seqhash):
		row = self.get_genome_data(acc)
		if row and row[0]['description'] == descr and row[0]['seqhash'] == seqhash:
			return True
		return False

	def sequence_exists(self, seqhash):
		sql = "SELECT COUNT(*) from sequence WHERE seqhash = ?"
		return self.rocon.execute(sql, (seqhash, )).fetchone()['COUNT(*)'] > 0

	def get_dna_varid(self, ref, alt, start, end):
		sql = "SELECT varid from dna WHERE ref = ? AND alt = ? AND start = ? and end = ?"
		row = self.rocon.execute(sql, (ref, alt, start, end)).fetchone()
		if row:
			return row['varid']
		return None

	def get_prot_varid(self, protein, locus, ref, alt, start, end):
		sql = "SELECT varid from prot WHERE protein = ? AND locus = ? AND ref = ? AND alt = ? AND start = ? and end = ?"
		row = self.rocon.execute(sql, (protein, locus, ref, alt, start, end)).fetchone()
		if row:
			return row['varid']
		return None

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
	print("sonarDB", sonarDB.get_version())
