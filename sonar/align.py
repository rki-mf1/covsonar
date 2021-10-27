#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

# DEPENDENCIES
import os
import re
import sys
import mappy
import subprocess
import tempfile

# CLASS
class sonarAligner(object):
	"""
	this object performs a pairwise sequence alignment and provides/stores selected
	alignment functionalities/statistics.
	"""

	def __init__(self):
		self.__cs_regex = re.compile(r':[0-9]+|\*[a-z][a-z]|[+-][A-Za-z]+')

	def map(self, ref_seq, qry_seq, cpus=1):
		#aligner = mappy.Aligner(seq=ref_seq, preset="asm20", n_threads=threads)
		#scoring list contains in the respective order:
		# -A Matching score [2]
		# -B Mismatching penalty [4]
		# -q Gap open penalty [4]
		# -e Gap extension penalty [2]
		# -q2 Long gap open penalty [24]
		# -e2 Long gap extension penalty [1]
		# -sc_ambi score involving ambiguous bases [1]

		aligner = mappy.Aligner(seq=ref_seq, scoring=[2, 4, 4, 2, 24, 1, 0], best_n=1, n_threads=cpus)
		if not aligner:
			sys.exit("aligner error: minimap2 failed to build index")

		hits = aligner.map(qry_seq, cs=True, MD=False)
		try:
			hit = next(hits)
		except:
			return None
		return [hit.r_st, hit.r_en, hit.q_st, hit.q_en, hit.cs]

	def itermap(self, ref_seq, qry_seq, cpus=1):
		mapping = self.map(ref_seq, qry_seq, cpus)
		if not mapping:
			return None

		ref_start, ref_end, qry_start, qry_end, cs = mapping

		ref_last_pos = len(ref_seq)-1
		qry_last_pos = len(qry_seq)-1

		# iterations
		while ref_end < ref_last_pos and qry_end < qry_last_pos:
			rseq = ref_seq[ref_end:]
			qseq = qry_seq[qry_end:]
			iteration = self.map(rseq, qseq, cpus)
			if not iteration:
				qry_end = qry_last_pos
				cs += "+" + qseq
				break
			if iteration[2] != 0:
				if iteration[0] == iteration[2]:
					for i in range(iteration[0]):
						#print("*"+ rseq[i].lower() + qseq[i].lower())
						cs += "*"+ rseq[i].lower() + qseq[i].lower()
				else:
					sys.exit("alignment error: iterative alignment exception")
			ref_end += iteration[1]
			qry_end += iteration[3]
			cs += iteration[4]

		return ref_start, ref_end, qry_start, qry_end,cs


	def minimap2(self, ref_file, qry_file):
		with open(ref_file) as handle:
			ref_seq = "".join([x.strip() for x in handle.readlines()[1:]])

		with open(qry_file) as handle:
			qry_seq = "".join([x.strip() for x in handle.readlines()[1:]])

		cmd = ['minimap2', '-t', str(cpus), '--cs', '-x', 'asm10', '--score-N=0', "--end-bonus=30000",  "-z", "10000" , ref_file, qry_file]
		p = subprocess.Popen(cmd, encoding='utf8', stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
		output, outerr = p.communicate()

		return output.split('\n')[0].split('\t')


	def iter_diffs(self, ref_seq, qry_seq, cpus=1, terminal_fill="."):
		alignment = self.itermap(ref_seq, qry_seq, cpus=cpus)
		if not alignment:
			sys.exit("alignment error: alignment failed.")

		ref_start, ref_end, qry_start, qry_end, cs = alignment

		# handling unaligned regions at the start of the query sequence
		if qry_start > 0:
			yield "", -1, qry_seq[:qry_start]
			if ref_start > 0:
				yield ref_seq[:ref_start], 0, ""
		elif ref_start > 0:
			for i in range(ref_start):
				yield ref_seq[i], i, terminal_fill

		# handling variations / snps, deletions or insertions within the query sequence
		pos = ref_start-1
		for match in self.__cs_regex.finditer(cs):
			match = match.group(0)
			if match[0] == ":":
				pos += int(match[1:])
			elif match[0] == "*":
				pos += 1
				yield match[1].upper(), pos, match[2].upper()
			elif match[0] == "+":
				yield "", pos, match[1:].upper()
			elif match[0] == "-":
				yield match[1:].upper(), pos+1, ""
				pos += len(match)-1

		#ending insert
		if qry_end != len(qry_seq):
			yield "", -1, qry_seq[:qry_start]
			if ref_start > 0:
				yield ref_seq[:ref_start], 0, ""

		#trailing N
		for i in range(ref_end, len(ref_seq)):
			yield ref_seq[i], i, terminal_fill
