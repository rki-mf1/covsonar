#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

VERSION = "0.0.9"
import os
import sys
import csv
import re
import argparse

def parse_args(cmd=None):
	parser = argparse.ArgumentParser(prog="markx3_to_diff.py", description="extracts variants from markx3 formatted pairwise alignment data")
	parser.add_argument('input_file', metavar="FILE", help="markx3 formatted pairwise alignment file", type=str)
	parser.add_argument('output_file', metavar="FILE", help="output file", type=str)
	parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
	return parser.parse_args(cmd)

def read_markx3(fname):
	with open(fname, "r") as handle:
		seqs = []
		for line in handle:
			if line.startswith("#") or len(line.strip()) == 0:
				continue
			elif line.startswith(">"):
				seqs.append([])
			else:
				seqs[-1].append(line.strip())
	return ["".join(x) for x in seqs]

def walk_algn(qry_seq, ref_seq):
	l = len(qry_seq)
	if l != len(ref_seq):
		sys.exit("error: sequences differ in length")
	qry_seq += " "
	ref_seq += " "
	i = 0
	offset = 0
	while i < l:
		#match
		if qry_seq[i] == ref_seq[i]:
			pass
		#deletion
		elif qry_seq[i] == "-":
			s = i
			while qry_seq[i+1] == "-":
				i += 1
			yield ref_seq[s:i+1], str(s-offset), str(i+1-offset), " "
		#insertion
		elif ref_seq[i] == "-":
			s = i-1
			while ref_seq[i+1] == "-":
				i += 1
			if s == -1:
				ref = " "
			else:
				ref = ref_seq[s]
			yield ref, str(s-offset), str(s-offset+1), qry_seq[s:i+1]
			offset += i-s
		#snps
		else:
			yield ref_seq[i], str(i-offset), str(i-offset+1), qry_seq[i]
		i += 1

def write_diff(aseq, bseq, fname):
	with open(fname, "w") as handle:
		handle.write("\n".join(["\t".join(x) for x in walk_algn(aseq.upper(), bseq.upper())]))

def main(cmd=None):
    args = parse_args(cmd)
    write_diff(*read_markx3(args.input_file), args.output_file)

if __name__ == "__main__":
	if "snakemake" in globals():
		cmd = [snakemake.input[0], snakemake.output[0]]
	else:
		cmd=None
	main(cmd=cmd)
