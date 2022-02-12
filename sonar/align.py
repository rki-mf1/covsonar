#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

# DEPENDENCIES
import os
import re
import sys
import pickle
import subprocess
from Bio.Emboss.Applications import StretcherCommandline


# CLASS
class sonarAligner(object):
	"""
	this object performs a pairwise sequence alignment and provides/stores selected
	alignment functionalities/statistics.
	"""

	def __init__(self):
		pass

	def align(self, qry, ref, gapopen= 16, gapextend = 4):
		"""
		"""
		cline = StretcherCommandline(asequence=qry, bsequence=ref, gapopen=gapopen, gapextend=gapextend, outfile='stdout', aformat="fasta")
		stdout, stderr = cline()
		s1=stdout.find("\n")+1
		e=stdout[1:].find(">")+1
		s2=stdout[e:].find("\n") + e
		qry = stdout[s1:e].replace("\n", "")
		ref = stdout[s2:].replace("\n", "")
		return qry, ref

	def process_cached_sample(self, fname):
		"""
		"""
		with open(fname, 'rb') as handle:
			data = pickle.load(handle, encoding="bytes")
		if os.path.isfile(data['var_file']):
			line = ""
			with open(data['var_file'], "r") as handle:
				for line in handle:
					pass
			if line == "//":
				return
		alignment = self.align(data['seq_file'], data['ref_file'])
		vars = "\n".join(["\t".join(x) for x in self.extract_vars(*alignment)])
		try:
			with open(data['var_file'], "w") as handle:
				handle.write(vars + "\n//")
		except:
			os.makedirs(os.path.dirname(data['var_file']))
			with open(data['var_file'], "w") as handle:
				handle.write(vars + "\n//")

	def extract_vars(self, qry_seq, ref_seq):
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
