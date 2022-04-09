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
from tqdm import tqdm


def parse_args():
	parser = argparse.ArgumentParser(prog="import.py", description="import sequences from fasta file to a sonar database")
	parser.add_argument('--fasta', help="fasta file", type=str, nargs="+", default=None)
	parser.add_argument('--tsv', help="tab-delimited sample property file", type=str, default=None)
	parser.add_argument('--cols', help="sample property file", type=str, nargs="+", default=None)
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

def blocks(files, size=65536):
	while True:
		b = files.read(size)
		if not b:
			break
		yield b


if __name__ == "__main__":
	args = parse_args()
	if args.tsv is None and args.fasta is None:
		print("Nothing to import.")
		exit(0)

	if args.tsv:
		with open(args.tsv) as handle:
			headline = [x.strip() for x in handle.readline().split("\t")]

		with sonarDBManager(args.db, debug=args.debug) as dbm:
			stored_propnames = set(list(dbm.properties.keys()))

		cols = {x: x for x in stored_propnames if x in headline}
		cols['sample'] = 'sample'

		if not args.cols is None:
			for x in args.cols:
				if x.count("=") != 1:
					sys.exit("input error: " + x + " is not a valid sample property assignment.")
				k, v = x.split("=")
				if k not in stored_propnames and k != "sample":
					sys.exit("input error: sample property " + k + " does not exist in the given database.")
				cols[k] = v

		for v in cols.values():
			if headline.count(v) != 1:
				sys.exit("input error: missing or ambiguous " + v + " column in tsv file provided.")

		cols = {x: headline.index(cols[x]) for x in cols}

		if len(cols) == 1:
			sys.exit("input error: tsv does not provide any informative column.")

		#with open("file", "r",encoding="utf-8", errors='ignore') as f:
		#	total_proplines = sum(bl.count("\n") for bl in blocks(f))

	temp = False if args.dir else True

	cache = sonarCache(args.db, outdir=args.dir, logfile="import.log", allow_updates=args.update, temp=temp, debug=args.debug)

	if not args.fasta is None:
		print("STAGE 1: caching data")
		cache.add_fasta(*args.fasta)
		print("STAGE 2: profiling genomes")
		aligner = sonarAligner()
		l = len(cache._samplefiles)
		with WorkerPool(n_jobs=args.threads, start_method='fork') as pool:
			results = pool.map(aligner.process_cached_sample, cache._samplefiles, iterable_len=l, n_splits=args.threads, progress_bar=True)
		print("STAGE 3: importing profiles")
		cache.import_samples()

	if args.tsv:
		print("STAGE 4: importing sample properties")
		with sonarDBManager(args.db, debug=args.debug) as dbm:
				with open(args.tsv) as handle:
					headline = handle.readline()
					for i, line in enumerate(tqdm(handle)):
						fields = [x.strip() for x in line.split("\t")]
						sample_name = fields[cols['sample']]
						sample_id = dbm.get_sample_id(sample_name)
						if not sample_id:
							continue
						for property_name, col_index in cols.items():
							if property_name == "sample":
								continue
							dbm.insert_property(sample_id, property_name, fields[col_index])
