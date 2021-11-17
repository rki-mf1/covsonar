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
from Bio.Emboss.Applications import StretcherCommandline, NeedleallCommandline
from Bio import SeqIO
from Bio.SeqIO import FastaIO

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
		#print(hit.r_st, hit.r_en, hit.q_st, hit.q_en, hit.cs)
		return [hit.r_st, hit.r_en, hit.q_st, hit.q_en, hit.cs]


	def itermap(self, ref_seq, qry_seq, cpus=1):
		#mapping = self.hisat2( qry_seq, cpus)
		#mapping = self.map(ref_seq, qry_seq, cpus)
		mapping = self.use_stretcher(qry_seq, ref_seq)
		# mapping = self.hisat2( qry_seq, cpus)

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

	def hisat2(self, qry_seq, cpus=1):
		tmp_dir = "/home/kongkitimanonk/SCRATCH_NOBAK/CovSonar1/workdir_covsonar/"
		ref_index_dir = "/scratch/kongkitimanonk/CovSonar1/workdir_covsonar/ref_index"

		# create an intermediate file for use as an input to hisat2
		n = 1
		qry_file = os.path.join(tmp_dir, "tmp_qry_file.fa")
		with open(qry_file, 'w') as out:
				out.write('>' + str(n) + '\n' + qry_seq.strip())  
		sam_qry_file = os.path.join(tmp_dir, "tmp_hist2.output.sam")
		cmd = ['/scratch/kongkitimanonk/CovSonar1/hisat2-2.2.1/hisat2', "--no-unal",  "-f", ref_index_dir, '--quiet', "--no-softclip", "--threads", str(cpus)
			   , "-U", qry_file, "-S", sam_qry_file, "-k", str(10)]
		
		p = subprocess.Popen(cmd, encoding='utf8')
		output, outerr = p.communicate()

		# Convert from SAM format to PAF format by using  paftools.js from minimap2
		cmd = ['paftools.js', 'sam2paf', sam_qry_file ]
		p = subprocess.Popen(cmd, encoding='utf8', stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		output, outerr = p.communicate()
		print('output',output)
		print('error', outerr)
		#os.remove(qry_file)
		#os.remove(sam_qry_file)
		list = output.split('\n')[0].split('\t')
		print(len(list))
		if len(list) <=1:
			return None
		# [5:] remove cs:Z:
		# index 7 8 2 3 17 = hit.r_st, hit.r_en, hit.q_st, hit.q_en, hit.cs
		# print(int(list[7]), int(list[8]), int(list[2]), int(list[3]), list[17][5:].lower())
	
		print(int(list[7]), int(list[8]), int(list[2]), int(list[3]), list[17][5:].lower())
		return  int(list[7]), int(list[8]), int(list[2]), int(list[3]), list[17][5:].lower()


	def use_stretcher(self, qry_seq, target_file, out_file = None, gapopen= 16, gapextend = 4, right_align = True):
		from tempfile import mkstemp, TemporaryDirectory, mkdtemp
		"""
		function to perform a pairwise aligment using EMBOSS Stretcher

		Parameters
		----------
		query_file : str
			define a path to a valid FASTA file storing the query sequence
		target_file : str
			define a path to a valid FASTA file storing the target sequence
			(= reference)
		out_file : str
			define a path to a file that will store the alignment. Please consider,
			that an existing file will be overwritten.
		gapopen : int [ 16 ]
			define penalty for gap opening
		gapextend : int [ 4 ]
			define penalty for gap extension

		Returns
		-------
		list
		  list of aligned query and target sequence, in that order
		"""
		tmp_dir = "/home/kongkitimanonk/SCRATCH_NOBAK/CovSonar1/workdir_covsonar/mycache"
		target_file=os.path.join(os.path.dirname(os.path.realpath(__file__)), "ref.fna")
		n=123
		fd, query_file = mkstemp(dir=tmp_dir, suffix=".fasta")
		with open(query_file, 'w') as out:
				out.write('>' + str(n) + '\n' + qry_seq.strip())  

		
		temp = True if not out_file else False
		if temp:
			handle, out_file = mkstemp(dir=tmp_dir, suffix=".1.sam")
		#cline = StretcherCommandline(asequence=query_file, bsequence=target_file,
		#							 gapopen=gapopen, gapextend=gapextend, outfile=out_file, aformat="fasta")
		#stdout, stderr = cline()
		#alignment = [str(x.seq) for x in SeqIO.parse(out_file, "fasta")]
		#if right_align:
		#	query_tmp, target_tmp = self.left_align_gaps(*alignment) # query, target
		#handle, out_file_query = mkstemp(dir=tmp_dir, suffix=".query.fasta")
		#with open(out_file_query, 'w') as out:
		#		out.write('>' + str(n) + '\n' + query_tmp)  
		#handle, out_target_query = mkstemp(dir=tmp_dir, suffix=".ref.fasta")
		#with open(out_target_query, 'w') as out:
		#		out.write('>' + "NC_045512.2" + '\n' + target_tmp)
		#print(query)
		#print(len(query))
		#print(target)
		#print(len(target))
		#cline = NeedleallCommandline(asequence=target_file, bsequence=query_file,
		#							 gapopen=gapopen, gapextend=gapextend, outfile=out_file, aformat="sam")
		#stdout, stderr = cline()
		#handle, out_file_2 = mkstemp(dir=tmp_dir, suffix=".2.sam")
		#cline = StretcherCommandline(asequence=out_target_query, bsequence=out_file_query,
		#							 gapopen=gapopen, gapextend=gapextend, outfile=out_file_2, aformat="sam")
		#stdout, stderr = cline()
		#print()


		cmd = ['stretcher',  "-asequence", target_file , '-bsequence',query_file , "-outfile", out_file
			   , "-gapopen", str(gapopen), "-gapextend", str(gapextend), '-aformat', "sam", "-aglobal" ]
		
		p = subprocess.Popen(cmd, encoding='utf8')
		output, outerr = p.communicate()
		
		
		## Fix header 
		with open(out_file, 'r') as original: data = original.read()
		with open(out_file, 'w') as modified: modified.write("@SQ	SN:NC_045512.2	LN:29903\n" + data)
		
		# Right align   maybe use bioconda and fix and write result back 


		# Convert from SAM format to PAF format by using  paftools.js from minimap2
		cmd = ['paftools.js', 'sam2paf', out_file ]
		p = subprocess.Popen(cmd, encoding='utf8', stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		output_paf, outerr_paf = p.communicate()
		print(output_paf)
		#if temp:
		#	os.remove(out_file)
		#	os.remove(query_file)

		root_list = output_paf.split('\n')
		
		print("list",root_list)
		if len(root_list) <=1:
			return None

		for i in root_list:
			list = i.split('\t')
			if len(list) <=1:
				continue
			else:
				print(list)
				print(int(list[7]), int(list[8]), int(list[2]), int(list[3]), list[16][5:].lower())
				return  int(list[7]), int(list[8]), int(list[2]), int(list[3]), list[16][5:].lower()

				

		# [5:] remove cs:Z:
		# index 7 8 2 3 17 = hit.r_st, hit.r_en, hit.q_st, hit.q_en, hit.cs
		# print(int(list[7]), int(list[8]), int(list[2]), int(list[3]), list[17][5:].lower())
		
	
	def left_align_gaps(self, query, target):
		"""
		function to align gaps to the left in two aligned sequences

		Parameters
		----------
		query : str
			define the query sequence in aligned form
		target : str
			define the target sequence (reference) in aligned form

		Returns
		-------
		list
		  aligned query and target sequence strings with left-aligned gaps,
		  in that order.
		"""
		l = len(query)-1
		for match in re.finditer("-+", query):
			s = match.start()-1
			e = match.end()-1
			g = "-" * (e-s)
			while s >= 0 and e < l and query[s] == target[e]:
				query = query[:s] + g + query[s] + query[e+1:]
				s -= 1
				e -= 1
		for match in re.finditer("-+", target):
			s = match.start()-1
			e = match.end()-1
			g = "-" * (e-s)
			while s >= 0 and e < l and target[s] == query[e]:
				target = target[:s] + g + target[s] + target[e+1:]
				s -= 1
				e -= 1
		return query, target


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
