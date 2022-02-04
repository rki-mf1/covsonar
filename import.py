#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

VERSION = "0.0.9"
import os
import sys
import argparse
from sonar import sonarActions, sonarDBManager, sonarCache, sonarAligner
import subprocess
from mpire import WorkerPool

def parse_args():
	parser = argparse.ArgumentParser(prog="import.py", description="import sequences from fasta file to a sonar database")
	parser.add_argument('--fasta', help="fasta file", type=str, nargs="+")
	parser.add_argument('--db', metavar="FILE", help="sonar database", type=str, required=True)
	parser.add_argument('--dir', metavar="DIR", help="cache directory", type=str, default=None)
	parser.add_argument('--threads', metavar="int", help="number of cpus to use (default: 1)", type=int, default=1)
	parser.add_argument('--update', help="allow updates on existing entries", action="store_true")
	parser.add_argument('--force', help="if using an existing cache, force all calculations to be rerun", action="store_true")
	parser.add_argument('--debug', help="activate debugging mode showing all sqllite queries on screen", action="store_true")
	parser.add_argument('--memsave', help="memory-friendly but slightly slower mode", action="store_true")

	# version
	parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)

	return parser.parse_args()


if __name__ == "__main__":
	args = parse_args()
	script_base = os.path.dirname(os.path.realpath(__file__))
	temp = False if args.dir else True
	cache = sonarCache(args.db, outdir=args.dir, logfile="import.log", allow_updates=args.update, temp=temp, debug=args.debug)
	cache.add_fasta(*args.fasta)

	aligner = sonarAligner()
	l = len(cache._samplefiles)
	with WorkerPool(n_jobs=args.threads, start_method='fork') as pool:
		results = pool.map(aligner.process_cached_sample, cache._samplefiles, iterable_len=l, n_splits=args.threads, progress_bar=True)

	cache.import_samples()

