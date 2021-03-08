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
from packaging import version
import shutil
import base64
from collections import OrderedDict

def dict_factory(cursor, row):
    d = OrderedDict()
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d

class sonarDB():
	def __init__(self, db, check_db=True):
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
		self.__max_supported_prev_version = "0.0.9"
		self.__conn = None
		self.__cursor = None
		self.__refseq = None
		self.__iupac_nuc_code = None
		self.__annotation = None
		self.__indel_regex = re.compile(".-+")
		if check_db:
			self.check_dir()

	def __exit__(self, type, value, traceback):
		if self.__conn:
			self.__conn.close()

	@property
	def conn(self):
		if not self.__conn:
			self.__conn = self.connect()
		return self.__conn

	@property
	def cursor(self):
		if not self.__cursor:
			self.__cursor = self.conn.cursor()
		return self.__cursor

	@property
	def refseq(self):
		if not self.__refseq:
			with open(self.reffna, "r") as handle:
				self.__refseq = "".join([x.strip().upper() for x in handle.readlines()[1:]])
		return self.__refseq

	@property
	def cds(self):
		if self.__annotation is None:
			gene_regex = re.compile("gene=([^;]+)(?:;|$)")
			with open(self.refgff, "r") as handle:
				self.__annotation = {}
				for line in handle:
					if not line.startswith("#"):
						continue
					fields = line.strip().split("\t")
					if fields[2] == "CDS":
						strand = fields[6]
						s = fields[3]
						e = fields[4]
						gene = gene_regex.search(fields[-1]).groups(1)
						self.__annotation[gene] = (gene, s, e, strand)
		return self.__annotation

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

	# directory management

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

	# data management

	def slugify_fname(self, fname):
		return os.path.join(self.seqdir, base64.urlsafe_b64encode(fname.encode('UTF-8') ).decode('UTF-8') + ".fna")

	def add_genome(self, acc, descr, seq, paranoid=False):
		seqhash = sonarDB.hash(seq)
		fname = self.slugify_fname(seqhash)
		if not os.path.isfile(fname):
			with open(fname, "w") as handle:
				handle.write(">" + seqhash + "\n" + seq.strip().upper() + "\n")
		elif paranoid:
			with open(fname, "r") as handle:
				handle.readline()
				if seq.upper() != handle.readline().strip():
					sys.exit("error: hash collision for " + hashval)
		self.insert_genome(acc, descr, seqhash)

	def compare_nucs(self, ref_nuc, qry_nuc):
		if qry_nuc in self.iupac_nuc_code[ref_nuc]:
			return True
		else:
			return False

	def add_alignment(self, aligned_target_seq, aligned_ref_seq, seqhash = None):
		if not seqhash:
			seqhash = sonarDB.hash(aligned_target_seq.replace("-", ""))
		vars = []

		# find insertions

		for match in self.__indel_regex.finditer(aligned_ref_seq):
			s = match.start()
			e = match.end()
			r = aligned_ref_seq[s]
			q = aligned_target_seq[s:e]
			gaps = aligned_ref_seq[:e].count("-")
			vars.append((s+1-gaps, None, r, q))

		## find deletions

		for match in self.__indel_regex.finditer(aligned_target_seq):
			s = match.start()
			e = match.end()
			r = aligned_ref_seq[s:e]
			q = aligned_target_seq[s]
			gaps = aligned_ref_seq[:e].count("-")
			vars.append((s+1-gaps, e-gaps, r, q))

		## find snps

		i = 0
		for r, q in zip(aligned_ref_seq, aligned_target_seq):
			if r == "-":
				continue
			elif self.compare_nucs(r, q) or q == "-":
				i += 1
				continue
			i += 1
			vars.append((i, None, r, q))


		# import to db
		for var in vars:
			self.insert_dna_var(seqhash, *var)

	# db management

	def connect(self):
		try:
			conn = sqlite3.connect(self.dbfile)
			conn.row_factory = dict_factory
		except Error as e:
			exit(e)
		return conn

	def close(self):
		if self.__conn:
			self.__conn.close()

	def commit(self):
		self.conn.commit()

	def insert_sequence(self, seqhash):
		self.cursor.execute('INSERT OR IGNORE INTO sequence(seqhash) VALUES(?)', (seqhash, ))

	def get_genome_data(self, acc):
		sql = "SELECT * from genome WHERE accession = '" + acc + "'"
		return self.cursor.execute(sql).fetchone()

	def insert_genome(self, acc, descr, seqhash):
		genome_data = self.get_genome_data(acc)
		if genome_data:
			if genome_data['description'] != descr or genome_data['seqhash'] != seqhash:
				sys.exit("error: data collision for " + acc + " use update to change data of an existing genome accession.")
		else:
			self.insert_sequence(seqhash)
			self.cursor.execute('INSERT INTO genome(accession, description, seqhash) VALUES (?, ?, ?)', (acc, descr, seqhash))

	def get_varid(self, start, end, ref, alt, protein=None):
		if protein is None:
			res = self.cursor.execute('SELECT varid from dna WHERE start = ? AND end = ? AND ref = ? AND alt = ?', (start, end, ref, alt))
		else:
			res = self.cursor.execute('SELECT varid from protein WHERE protein = ? AND start = ? AND end = ? AND ref = ? AND alt = ?', (protein, start, end, ref, alt))
		row = res.fetchone()
		if not row:
			return False
		else:
			return row['varid']

	def insert_dna_var(self, seqhash, start, end, ref, alt):
		varid = self.get_varid(start, end, ref, alt)
		if not varid:
			self.cursor.execute('INSERT INTO dna(start, end, ref, alt) VALUES (?, ?, ?, ?)', (start, end, ref, alt))
			varid = self.cursor.lastrowid
		self.cursor.execute('INSERT OR IGNORE INTO sequence2dna(seqhash, varid) VALUES (?, ?)', (seqhash, varid))

	def iter_rows(self, table):
		sql = "SELECT * FROM " + table
		for row in self.conn.cursor().execute(sql):
			yield row

	def exec(self, sql):
		try:
			self.conn.cursor().execute(sql)
		except Error as e:
			exit(e)

	def multiexec(self, multisql):
		for sql in multisql.split(";"):
			self.exec(sql)

	def create_tables(self):
		with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "resources", "db.sqlite"), 'r') as handle:
			sql = handle.read()
		self.multiexec(sql)
		self.commit()

	@staticmethod
	def get_version():
		with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), ".version"), "r") as handle:
			return handle.read().strip()

if __name__ == "__main__":
	print("sonarDB", sonarDB.get_version())
