#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

VERSION = "0.0.9"
import os
import sys
import argparse
import sqlite3
from sqlite3 import Error
from Bio.SeqUtils.CheckSum import seguid as bioseguid

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

	@property
	def conn(self):
		if not self.__conn:
			self.__conn = self.connect()
		return self.__conn

	@staticmethod
	def seguid(seq):
		return bioseguid(seq)

	def connect(self):
		try:
			conn = sqlite3.connect(self.db_file)
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

	def insert_genome(self, acc, descr, seq):
		genome_data = self.get_genome_data(acc)
		seguid = sonarDB.seguid(seq)
		if genome_data:
			if genome_data['descr'] != descr or genome_data['seguid'] != seguid:
				sys.exit("error: data collision for " + acc + " use update to change data of an existing genome accession.")
		else:
			self.insert_sequence(seguid)
			cursor = self.conn.cursor()
			cursor.execute('INSERT INTO genome(accession, description, seguid) VALUES (?, ?, ?)', (acc, descr, seguid))

	def iter_rows(self, table):
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
		for row in db.iter_rows("dna_view"):
			print(row)
		#db.exec("INSERT INTO sequence(seq) VALUES('ATG')")
		#print(db.show_all("sequence"))
