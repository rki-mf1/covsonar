#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

import os
import sys
import argparse
import sqlite3
from sqlite3 import Error
from Bio.SeqUtils.CheckSum import seguid
from packaging import version
import shutil
import base64

def dict_factory(cursor, row):
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d

def parse_args():
	parser = argparse.ArgumentParser(prog="sonar.py", description="")
	parser.add_argument('--db', metavar="STR", help="existing or new database file (default: sonar.db)", type=str, default="sonar.db")
	parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
	return parser.parse_args()

class sonarDB():
	def __init__(self, db):
		self.dbdir = os.path.abspath(db)
		self.vfile = os.path.join(self.dbdir, ".version")
		self.dbfile = os.path.join(self.dbdir, "sonar.db")
		self.datadir = os.path.join(self.dbdir, "data")
		self.seqdir = os.path.join(self.datadir, "seqs")
		self.__max_supported_prev_version = "0.0.9"
		self.__conn = None

		self.check_dir()
		self.create_tables()

	def __exit__(self, type, value, traceback):
		if self.__conn:
			self.__conn.close()

	@property
	def conn(self):
		if not self.__conn:
			self.__conn = self.connect()
		return self.__conn

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
			os.makedirs(self.seqdir, exist_ok=True)
			open(self.dbfile, "w").close()
			with open(self.vfile, "w") as handle:
				handle.write(sonarDB.get_version())
		except Error as e:
			shutil.rmtree(self.dbdir)
			exit(e)

	# data management

	def slugify_seq_fname(self, seqhash):
		return os.path.join(self.seqdir, base64.urlsafe_b64encode(seqhash.encode('UTF-8') ).decode('UTF-8') + ".fna")

	def add(self, acc, descr, seq, ref, paranoid=False):
		seqhash = sonarDB.hash(seq)
		fname = self.slugify_seq_fname(seqhash)
		if not os.path.isfile(fname):
			with open(fname, "w") as handle:
				handle.write(">" + seqhash + "\n" + seq.strip().upper() + "\n" + ref)
		elif paranoid:
			with open(fname, "r") as handle:
				handle.readline()
				if seq.upper() != handle.readline().strip():
					sys.exit("error: hash collision for " + hashval)
		self.insert_genome(acc, descr, seqhash)

	# db management

	def connect(self):
		try:
			conn = sqlite3.connect(self.dbfile)
			conn.row_factory = dict_factory
		except Error as e:
			exit(e)
		return conn

	def insert_sequence(self, seguid):
		cursor = self.conn.cursor()
		cursor.execute('INSERT OR IGNORE INTO sequence(seguid) VALUES(?)', (seguid, ))

	def get_genome_data(self, acc):
		sql = "SELECT * from genome WHERE accession = '" + acc + "'"
		return self.conn.cursor().execute(sql).fetchone()

	def insert_genome(self, acc, descr, seqhash):
		genome_data = self.get_genome_data(acc)
		if genome_data:
			if genome_data['descr'] != descr or genome_data['seguid'] != seqhash:
				sys.exit("error: data collision for " + acc + " use update to change data of an existing genome accession.")
		else:
			self.insert_sequence(seqhash)
			cursor = self.conn.cursor()
			cursor.execute('INSERT INTO genome(accession, description, seguid) VALUES (?, ?, ?)', (acc, descr, seqhash))

	def iter_rows(self, table):
		sql = "SELECT * FROM " + table
		for row in self.conn.cursor().execute(sql):
			yield row

	def exec(self, sql):
		try:
			self.conn.cursor().execute(sql)
			self.conn.commit()
		except Error as e:
			exit(e)

	def multiexec(self, multisql):
		for sql in multisql.split(";"):
			self.exec(sql)

	def create_tables(self):
		with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "resources", "db.sqlite"), 'r') as handle:
			sql = handle.read()
		self.multiexec(sql)

	@staticmethod
	def get_version():
		with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), ".version"), "r") as handle:
			return handle.read().strip()

if __name__ == "__main__":
	print("sonarDB", sonarDB.get_version())
	#args = parse_args()
	#with sonarDB(args.db) as db:
	#	db.insert_genome("test1", "a test", "ATG")
	#	for row in db.iter_rows("dna_view"):
	#		print(row)
