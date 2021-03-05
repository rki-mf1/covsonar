#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

VERSION = "0.0.9"
import os
import sys
import argparse
import sqlite3
from sqlite3 import Error

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
	def __init__(self, db_file):
		self.db_file = db_file
		self.__conn = None
	def __enter__(self):
		self.create_tables()
		return self
	def __exit__(self, type, value, traceback):
		if self.__conn:
			self.__conn.close()


	def connect(self):
		try:
			conn = sqlite3.connect(self.db_file)
			conn.row_factory = dict_factory
		except Error as e:
			exit(e)
		return conn

	@property
	def conn(self):
		if not self.__conn:
			self.__conn = self.connect()
		return self.__conn

	def get_seqID(self, seq, acc=None):
		if seq:
			sql = "SELECT seqID from sequence WHERE seq = '" + seq + "'"
		elif acc:
			sql = "SELECT seqID from genome WHERE acc = '" + acc + "'"
		else:
			return None
		return self.conn.cursor().execute(sql).fetchone()

	def get_genome_data(self, acc):
		sql = "SELECT * from genome WHERE accession = '" + acc + "'"
		return self.conn.cursor().execute(sql).fetchone()

	def insert_sequence(self, seq):
		seqID = self.get_seqID(seq)
		if seqID is None:
			cursor = self.conn.cursor()
			cursor.execute('INSERT INTO sequence VALUES (?, ?)', (None, seq))
			seqID = cursor.lastrowid
		else:
			seqID = seqID['seqID']
		return seqID

	def insert_genome(self, acc, descr, seq, nuc_vars=None, aa_vars=None):
		genome_data = self.get_genome_data(acc)
		if genome_data:
			if genome_data['descr'] != descr or acc_seqID != self.get_seqID(seq=seq):
				sys.exit("error: data collision for " + acc + " use update to change data of an existing genome accession.")
			return genome_data['genomeID']
		else:
			seqID = self.insert_sequence(seq)
			vals = [acc, descr, seqID]
			cursor = self.conn.cursor()
			cursor.execute('INSERT INTO genome VALUES (?, ?, ?, ?)', (None, acc, descr, seqID))
			genomeID = cursor.lastrowid
			return genomeID

	def iter_genomes(self, table):
		sql = "SELECT * FROM " + table
		try:
			rows = self.conn.cursor().execute(sql)
		except Error as e:
			exit(e)
		for row in rows:
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


if __name__ == "__main__":
	args = parse_args()
	with sonarDB(args.db) as db:
		#print(db.insert_sequence("ATG"))
		#print(db.insert_sequence("ATC"))
		db.insert_genome("test1", "a test", "ATG")
		for row in db.iter_genomes("global"):
			print(row)
		#db.exec("INSERT INTO sequence(seq) VALUES('ATG')")
		#print(db.show_all("sequence"))
