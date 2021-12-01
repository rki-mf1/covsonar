#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

VERSION = "0.0.9"
import os
import sys
import csv
import argparse
import gzip
import lzma
from .....sonar import sonarActions
from .....sonar.dbm import sonarDBManager
from .....sonar.cache import sonarCache
from Bio import SeqIO
from tempfile import mkstemp
from collections import defaultdict
import re
from tqdm import tqdm
from multiprocessing import Pool
import base64

def parse_args():
	parser = argparse.ArgumentParser(prog="add_fasta.py", description="prepares fasta entries for import")
	parser.add_argument('fasta', help="fasta file", action="str", nargs="+")
	parser.add_argument('--db', metavar="FILE", help="sonar database", type=str, required=True)
	parser.add_argument('--dir', metavar="DIR", help="cache directory", type=str, default=None)
	parser.add_argument('--cpus', metavar="int", help="number of cpus to use (default: 1)", type=int, default=1)
	parser.add_argument('--update', help="allow updates on existing entries", action="store_true")
	parser.add_argument('--debug', help="activate debugging mode showing all sqllite queries on screen", action="store_true")

	# version
	parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)

	return parser.parse_args()


if __name__ == "__main__":
	args = arse_args()
	temp = False if args.dir else True
	cache = sonarCache(args.db, outdir=args.dir, logfile="import.log", allow_updates=args.update, temp=temp, debug=args.debug)
	cache.add_fasta(*args.fasta)
