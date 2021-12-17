#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

import os
import re
import sys
import argparse
import base64
import pickle
from sonar import sonarActions, sonarDBManager
import hashlib
import yaml
from Bio import SeqIO
import difflib as dl
from math import ceil

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
			self.default_refmol_id = [x for x in self.refmols if self.refmols[x]['molecule.id'] == 1 ][0]
			self.sources = {x['molecule.accession']: dbm.get_source(x['molecule.id']) for x in self.refmols.values()}
			self.properties = dbm.properties
		self._propregex = re.compile("\[(" + "|".join(self.properties.keys()) + ")=([^\[\]=]+)\]")
		self._molregex = re.compile("\[molecule=([^\[\]=]+)\]")

		self.logfile = open(logfile, "w") if logfile else None
		self.basedir = os.path.abspath(mkdtemp(prefix=".sonarCache_")) if not outdir else os.path.abspath(outdir)
		self.smk_config = os.path.join(self.basedir, "config.yaml")
		self.sample_dir = os.path.join(self.basedir, "samples")
		self.seq_dir = os.path.join(self.basedir, "seq")
		self.algn_dir = os.path.join(self.basedir, "algn")
		self.var_dir = os.path.join(self.basedir, "var")
		self.ref_dir = os.path.join(self.basedir, "ref")

		os.makedirs(self.basedir, exist_ok=True)
		os.makedirs(self.seq_dir, exist_ok=True)
		os.makedirs(self.ref_dir, exist_ok=True)
		#os.makedirs(self.algn_dir, exist_ok=True)
		os.makedirs(self.var_dir, exist_ok=True)
		os.makedirs(self.sample_dir, exist_ok=True)
		self._samplefiles = set()
		self._refs = set()

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
	def file_collision(fname, data):
		with open(fname, "r") as handle:
			if handle.read() != data:
				return True
		return False

	def sample_collision(self, key, datadict):
		if sonarCache.read_pickle(key) != datadict:
			return True
		return False

	def write_smk_config(self):
		data = {
			'debug': self.debug,
			'sample_dir': self.sample_dir,
			'seq_dir': self.seq_dir,
			'algn_dir': self.algn_dir,
			'var_dir': self.var_dir
		}
		with open(self.smk_config, 'w') as handle:
		    yaml.dump(data, handle)

	def get_seq_fname(self, seqhash):
		fn = self.slugify(seqhash)
		return os.path.join(self.seq_dir, fn[:2], fn + ".seq")

	def get_ref_fname(self, refid):
		return os.path.join(self.ref_dir, str(refid) + ".seq")

	def get_algn_fname(self, seqhash):
		fn = self.slugify(seqhash)
		return os.path.join(self.algn_dir, fn[:2], fn + ".algn")

	def get_var_fname(self, seqhash):
		fn = self.slugify(seqhash)
		return os.path.join(self.var_dir, fn[:2], fn + ".var")

	def get_sample_fname(self, sample_name):
		fn = self.slugify(hashlib.sha1(sample_name.encode('utf-8')).hexdigest())
		return os.path.join(self.sample_dir, fn[:2], fn + ".sample")

	def cache_sample(self, name, sampleid, seqhash, header, refmol, refmolid, algnid, seqfile, reffile, algnfile, varfile, properties):
		data = {
			"name": name,
			"sampleid": sampleid,
			"refmol": refmol,
			"refmolid": refmolid,
			"algnid": algnid,
			"header": header,
			"seqhash": seqhash,
			"seq_file": seqfile,
			"ref_file": reffile,
			"algn_file": algnfile,
			"var_file": varfile
			}
		fname = self.get_sample_fname(name)
		if os.path.isfile(fname):
			if self.sample_collision(fname, data):
				sys.exit("sample data collision: data differs for sample " + name + " (" + header + ").")
		else:
			try:
				self.write_pickle(fname, data)
			except:
				os.makedirs(os.path.dirname(fname))
				self.write_pickle(fname, data)
		self._samplefiles.add(fname)
		return fname

	def iter_samples(self):
		for fname in self._samplefiles:
			yield self.read_pickle(fname)

	def cache_sequence(self, seqhash, sequence):
		fname = self.get_seq_fname(seqhash)
		if os.path.isfile(fname):
			if self.file_collision(fname, sequence):
				sys.exit("seqhash collision: sequences differ for seqhash " + seqhash + ".")
		else:
			try:
				with open(fname, "w") as handle:
					handle.write(sequence)
			except:
				os.makedirs(os.path.dirname(fname))
				with open(fname, "w") as handle:
					handle.write(sequence)
		return fname

	def cache_reference(self, refid, sequence):
		fname = self.get_ref_fname(refid)
		if refid not in self._refs:
			self._refs.add(refid)
			with open(fname, "w") as handle:
				handle.write(sequence)
		return fname

	def iter_fasta(self, *fnames):
		for fname in fnames:
			for record in SeqIO.parse(fname, "fasta"):
				refmol = self.get_refmol(record.description)
				if not refmol:
					sys.exit("input error: " +  record.id + " refers to an unknown reference molecule (" + self._molregex.search(fasta_header) + ").")
				refmolid = self.refmols[refmol]['molecule.id']
				seq = sonarActions.harmonize(record.seq)
				seqhash = sonarActions.hash(seq)
				yield {
					   'name': record.id,
					   'header': record.description,
					   'seqhash': sonarActions.hash(seq),
					   'sequence': seq,
					   'refmol': refmol,
					   'refmolid': refmolid,
					   'properties': self.get_properties(record.description)
					   }

	def get_refmol(self, fasta_header):
		mol = self._molregex.search(fasta_header)
		if not mol:
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
				if data['algnid'] is None:
					data['seqfile'] = self.cache_sequence(data['seqhash'], data['sequence'])
					data['reffile'] = self.cache_reference(refseq_id, self.get_refseq(data['refmol']))
					data['algnfile'] = self.get_algn_fname(data['seqhash'] + "@" + self.get_refhash(data['refmol']))
					data['varfile'] = self.get_var_fname(data['seqhash'] + "@" + self.get_refhash(data['refmol']))
				else:
					data['seqfile'] = None
					data['reffile'] = None
					data['algnfile'] = None
					data['varfile'] = None
				del(data['sequence'])
				self.cache_sample(**data)

	def import_samples(self):
		refseqs = {}
		with sonarDBManager(self.db, debug=False) as dbm:
			for sample_data in self.iter_samples():
				# nucleotide level import
				sampid = dbm.insert_sample(sample_data['name'])
				dbm.insert_sequence(sampid, sample_data['seqhash'])
				if sample_data['algnid'] is None:
					algnid = dbm.insert_alignment(sample_data['seqhash'], sample_data['refmolid'])
					with open(sample_data['var_file'], "r") as handle:
						for line in handle:
							vardat = line.strip("\r\n").split("\t")
							dbm.insert_variant(algnid, vardat[0], vardat[3], vardat[1], vardat[2])

				# paranoia test
				try:
					seq = list(refseqs[sample_data['refmolid']])
				except:
					refseqs[sample_data['refmolid']] = list(dbm.get_sequence(sample_data['refmolid']))
					seq = list(refseqs[sample_data['refmolid']])

				for vardata in dbm.iter_dna_variants(sample_data['name'], sample_data['refmolid']):
					#sys.stderr.write(str(vardata))
					#sys.stderr.write("\n")
					#if vardata['variant.start'] != -1 and vardata['variant.ref'] != seq[vardata['variant.start']]:
					#	sys.exit("error: paranoid1 (expected " + vardata['variant.ref'] + " got " + seq[vardata['variant.start']] + " at reference position " + str(vardata['variant.start']) + ")")
					if vardata['variant.alt'] == " ":
						for i in range(vardata['variant.start'], vardata['variant.end']):
							seq[i] = ""
					else:
						seq[vardata['variant.start']] = vardata['variant.alt']
				seq = "".join(seq)
				with open(sample_data['seq_file'], "r") as handle:
					orig_seq = handle.read()
				if seq != orig_seq:
					for x, y in sample_data.items():
						if x != "sequence":
							print(x + ":", y)
							print()
					print()
					for line in dl.ndiff([sample_data['sequence']], [seq]):
						print(line)
					sys.exit("error: paranoid2 caused " + sample_data['var_file'])

if __name__ == "__main__":
	pass
