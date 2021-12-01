#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

import os
import re
import sys
import argparse
import sqlite3
from sqlite3 import Error
from Bio.SeqUtils.CheckSum import seguid
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Emboss.Applications import StretcherCommandline
from packaging import version
import shutil
import base64
from collections import OrderedDict, defaultdict
import pickle
from tqdm import tqdm
from urllib.parse import quote as urlquote
from math import floor, ceil
from tempfile import mkstemp, TemporaryDirectory, mkdtemp
import traceback
import itertools
import signal
import csv
from time import sleep
from contextlib import ExitStack
from more_itertools import consecutive_groups, split_when
from sonar import sonarActions, sonarDBManager
import hashlib

class sonarCache():
	"""

	"""
	def __init__(self, db=None, outdir=None, refacc=None, logfile= None, allow_updates=True, ignore_errors=False, temp=False, debug=False):
		self.db = db
		self.allow_updates = allow_updates
		self.debug = debug
		self.refacc = refacc
		self.temp = temp
		self.ignore_errors = ignore_errors

		with sonarDBManager(self.db, debug=self.debug) as dbm:
			self.refmols = dbm.get_molecule_data("\"molecule.accession\"", "\"molecule.id\"", "\"molecule.standard\"" , reference_accession=self.refacc)
			self.default_refmol_acc = [x for x in self.refmols if self.refmols[x]['molecule.standard'] == 1 ][0]
			self.sources = {x['molecule.accession']: dbm.get_source(x['molecule.id']) for x in self.refmols.values()}
			self.properties = dbm.properties
		self._propregex = re.compile("\[(" + "|".join(self.properties.keys()) + ")=([^\[\]=]+)\]")
		self._molregex = re.compile("\[molecule=([^\[\]=]+)\]")

		self.logfile = open(logfile, "w") if logfile else None
		self.basedir = mkdtemp(prefix=".sonarCache_")if not outdir else outdir
		self.sample_dir = os.path.join(self.basedir, "samples")
		self.fasta_dir = os.path.join(self.basedir, "fasta")
		self.sample_subdirs = set()
		self.fasta_subdirs = set()
		os.makedirs(self.basedir, exist_ok=True)
		os.makedirs(self.sample_dir, exist_ok=True)
		os.makedirs(self.fasta_dir, exist_ok=True)

	def __enter__(self):
		return self

	def __exit__(self, exc_type, exc_value, exc_traceback):
		#if [exc_type, exc_value, exc_traceback].count(None) != 3:
		#	print("warning:", file=sys.stderr)
		#	print(traceback.format_exc(), file=sys.stderr)
		if os.path.isdir(self.basedir) and self.temp:
			shutil.rmtree(self.basedir)
		if self.logfile:
			self.logfile.close()
		if self.sample_yaml:
			self.sample_yaml.close()

	@staticmethod
	def slugify(string):
		return base64.urlsafe_b64encode(string.encode('UTF-8')).decode('UTF-8').rstrip("=")

	@staticmethod
	def deslugify(string):
		while len(string)%3 != 0:
			string += "="
		return base64.urlsafe_b64decode(string).decode("utf-8")

	def log(self, msg, die=False, errtype = "error"):
		if self.logfile:
			self.logfile.write(msg)
		elif not die:
			sys.stderr(msg)
		else:
			exit(errtype + ": " + msg)

	@staticmethod
	def write_pickle(fname, data):
		with open(fname, 'wb') as handle:
			pickle.dump(data, handle)

	@staticmethod
	def read_pickle(fname):
		with open(fname, 'rb') as handle:
			return pickle.load(handle, encoding="bytes")

	@staticmethod
	def pickle_collision(fname, data):
		if os.path.isfile(fname) and load_pickle(fname) != data:
			return True
		return False

	@staticmethod
	def load_pickle(fname):
		with open(fname, 'rb') as handle:
			return pickle.load(handle, encoding="bytes")

	@staticmethod
	def fasta_collision(fname, data):
		with open(fname, "r") as handle:
			if handle.read() != data:
				return True
		return False

	@staticmethod
	def sample_collision(fname, datadict):
		if sonarCache.read_pickle(fname) != datadict:
			return True
		return False

	def get_fasta_fname(self, seqhash):
		seqhash = self.slugify(seqhash)
		path = os.path.join(self.fasta_dir, seqhash[:2])
		os.makedirs(path, exist_ok=True)
		return os.path.join(path, seqhash + ".seq")

	def get_sample_fname(self, fasta_header):
		fbasename = self.slugify(hashlib.sha1(fasta_header.encode('utf-8')).hexdigest())
		path = os.path.join(self.sample_dir, fbasename[:2])
		os.makedirs(path, exist_ok=True)
		return os.path.join(path, hashlib.sha1(fasta_header.encode('utf-8')).hexdigest() + ".sample")

	def cache_sample(self, name, sampleid, seqhash, header, sequence, refmol, algnid, fasta, properties):
		fname = self.get_sample_fname(header)
		data = {
			"name": name,
			"sampleid": sampleid,
			"sequence": sequence,
			"refmol": refmol,
			"algnid": algnid,
			"header": header,
			"seqhash": seqhash,
			"fasta": os.path.abspath(fasta),
			"properties": properties,
			}
		if os.path.isfile(fname):
			if self.sample_collision(fname, data):
				sys.exit("sample data collision: data differs for sample " + name + " (" + header + ").")
		else:
			self.write_pickle(fname, data)
		return fname

	def cache_sequence(self, seqhash, sequence, refhash, refseq):
		fname = self.get_fasta_fname(seqhash + "-" + refhash)
		content = ">" + seqhash + "\n" + sequence + "\n>REF_" + refhash + "\n" + refseq
		if os.path.isfile(fname):
			if self.fasta_collision(fname, content):
				sys.exit("seqhash collision: sequences differ for seqhash " + seqhash + ".")
		else:
			with open(fname, "w") as handle:
				handle.write(content)
		return fname

	def iter_fasta(self, *fnames):
		for fname in fnames:
			for record in SeqIO.parse(fname, "fasta"):
				refmol = self.get_refmol(record.description)
				if not refmol:
					sys.exit("input error: " +  record.id + " refers to an unknown reference molecule (" + self._molregex.search(fasta_header) + ").")
				seq = sonarActions.harmonize(record.seq)
				seqhash = sonarActions.hash(seq)
				yield {
					   'name': record.id,
					   'header': record.description,
					   'seqhash': sonarActions.hash(seq),
					   'sequence': seq,
					   'refmol': refmol,
					   'properties': self.get_properties(record.description)
					   }

	def get_refmol(self, fasta_header):
		if mol := self._molregex.search(fasta_header):
			try:
				return self.refmols[mol]['accession']
			except:
				None
		return self.default_refmol_acc

	def get_refseq(self, refmol_acc):
		try:
			return self.sources[refmol_acc]['sequence']
		except:
			return None

	def get_refseq_id(self, refmol_acc):
		try:
			return self.sources[refmol_acc]['id']
		except:
			return None

	def get_refhash(self, refmol_acc):
		try:
			if 'seqhash' not in self.sources[refmol_acc]:
				self.sources[refmol_acc]['seqhash'] = sonarActions.hash(self.sources[refmol_acc]['sequence'])
			return self.sources[refmol_acc]['seqhash']
		except:
			return None

	def get_properties(self, fasta_header):
		return { x.group(1): x.group(2) for x in self._propregex.finditer(fasta_header) }

	def add_fasta(self, *fnames):
		with sonarDBManager(self.db, debug=False) as dbm:
			for data in self.iter_fasta(*fnames):
				# check sample
				data['sampleid'] = dbm.get_sample_id(data['name'])

				# check properties
				if data['sampleid'] is None:
					data['properties'].update({ x: self.properties[x]['standard']  for x in self.properties if not x in data['properties'] })
				elif not self.allow_updates:
					self.log(data['name'] + " skipped as it exists in the database and updating is disabled")
					continue
				else:
					stored_properties = dbm.get_properties(sample_name)
					data['properties'] = { x: data['properties'][x] for x in data['properties'].items() if stored_properties[x] != y }

				# check reference
				refseq_id = self.get_refseq_id(data['refmol'])
				if not refseq_id:
					if not self.ignore_errors:
						self.log("fasta header refers to an unknown refrence (" + data['header'] + ")", True, "input error")
					else:
						self.log("skipping " + data['name'] + " referring to an unknown reference (" + data['header'] + ")")

				# check alignment
				data['algnid'] = dbm.get_alignment_id(data['seqhash'], refseq_id)
				fasta = None if not data['algnid'] is None else self.cache_sequence(data['seqhash'], data['sequence'], self.get_refhash(data['refmol']), self.get_refseq(data['refmol']))
				self.cache_sample(**data, fasta=fasta)

if __name__ == "__main__":
	pass
