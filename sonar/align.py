#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

# DEPENDENCIES
import os
import re
import sys
import pickle
import subprocess
import pandas as pd
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

		sourceid = str(data['sourceid'])
		if os.path.isfile(data['var_file']):
			line = ""
			with open(data['var_file'], "r") as handle:
				for line in handle:
					pass
			if line == "//":
				return
		alignment = self.align(data['seq_file'], data['ref_file'])
		nuc_vars = [x for x in self.extract_vars(*alignment, sourceid)]
		vars = "\n".join(["\t".join(x) for x in nuc_vars])
		if nuc_vars:
			aa_vars = "\n".join(["\t".join(x) for x in self.lift_vars(nuc_vars, data['lift_file'], data['tt_file'])])
			if aa_vars:
				vars += "\n" + aa_vars
		try:
			with open(data['var_file'], "w") as handle:
				handle.write(vars + "\n//")
		except:
			os.makedirs(os.path.dirname(data['var_file']))
			with open(data['var_file'], "w") as handle:
				handle.write(vars + "\n//")

	def extract_vars(self, qry_seq, ref_seq, elemid):
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
				start = s-offset
				end = i+1-offset
				if end-start == 1:
					label = "del:" + str(start+1)
				else:
					label = "del:" + str(start+1) + "-" + str(end)
				yield ref_seq[s:i+1], str(start), str(end), " ", elemid, label

			#insertion
			elif ref_seq[i] == "-":
				s = i-1
				while ref_seq[i+1] == "-":
					i += 1
				if s == -1:
					ref = " "
					alt = qry_seq[:i+1]
				else:
					ref = ref_seq[s]
					alt = qry_seq[s:i+1]
				pos = s-offset+1
				yield ref, str(pos-1), str(pos), alt, elemid, ref + str(pos) + alt
				offset += i-s
			#snps
			else:
				ref = ref_seq[i]
				alt = qry_seq[i]
				pos = i-offset+1
				yield ref, str(pos-1), str(pos), alt, elemid, ref + str(pos) + alt
			i += 1

	def translate(self, seq, tt):
		aa = []
		while len(seq)%3 != 0:
			seq = seq[:len(seq)-1]
		for codon in [seq[i:i+3] for i in range(0, len(seq), 3)]:
			aa.append(tt[codon])
		return "".join(aa)

	def lift_vars(self, nuc_vars, lift_file, tt_file):
		df = pd.read_pickle(lift_file)

		with open(tt_file, 'rb') as handle:
			tt = pickle.load(handle, encoding="bytes")
		for nuc_var in nuc_vars:
			for i in range(int(nuc_var[1]), int(nuc_var[2])):
				alt = "-" if nuc_var[3] == " " else nuc_var[3]
				df.loc[df['nucPos1'] == i, 'alt1'] = alt
				df.loc[df['nucPos2'] == i, 'alt2'] = alt
				df.loc[df['nucPos3'] == i, 'alt3'] = alt

		df = df.loc[(df["ref1"] != df["alt1"]) | (df["ref2"] != df["alt2"]) | (df["ref3"] != df["alt3"])]
		df["altAa"] = df.apply(lambda x: self.translate(x["alt1"]+x["alt1"]+x["alt1"], tt), axis=1)
		df = df.loc[df["aa"] != df["altAa"]]

		# snps or inserts
		for index, row in df.loc[(df["altAa"] != "-") & (df["altAa"] != "")].iterrows():
			pos = row["aaPos"]+1
			label = row["aa"] + str(pos) + row['altAa']
			yield row["aa"], str(pos-1), str(pos), row['altAa'], str(row["elemid"]), label

		# deletions
		prev_row = None
		for index, row in df.loc[(df["altAa"] == "-")].sort_values(['elemid', 'aaPos']).iterrows():
			if prev_row is None:
				prev_row = row
			elif prev_row["elemid"] == row["elemid"] and abs(prev_row["aaPos"]-row['aaPos'])==1:
				prev_row["aa"] += row['aa']
			else:
				start = prev_row["aaPos"]
				end = prev_row["aaPos"]+len(prev_row["aa"])
				if end-start == 1:
					label = "del:" + str(start+1)
				else:
					label = "del:" + str(start+1) + "-" + str(end)
				yield prev_row["aa"], str(start), str(end), " ", str(prev_row["elemid"]), label
				prev_row = row
		if not prev_row is None:
			start = prev_row["aaPos"]
			end = prev_row["aaPos"]+len(prev_row["aa"])
			if end-start == 1:
				label = "del:" + str(start+1)
			else:
				label = "del:" + str(start+1) + "-" + str(end)
			yield prev_row["aa"], str(start), str(end), " ", str(prev_row["elemid"]), label
