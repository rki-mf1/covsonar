#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

VERSION = "0.0.9"
import os
import sys
import argparse
import sqlite3
from sqlite3 import Error

def parse_args():
	parser = argparse.ArgumentParser(prog="sonar.py", description="")
    parser.add_argument('--add', metavar="STR", help="sequence(s) from this file will be added to the splitted fasta files", type=str)
    parser.add_argument('--outdir', metavar="STR", help="output_directory", type=str)
    parser.add_argument('--check', help="activate consistency check of existing files (default: overwrite existing files)", action="store_true")
	parser.add_argument('fasta', metavar="STR", help="fasta file containing one or more entries to split", type=str, default="sonar.db")
	parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
	return parser.parse_args()

def fasta2dict():


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
