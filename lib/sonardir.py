#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

import os
import sys
import argparse
from packaging import version
import sonardb

class sonarDIR():
	def __init__(self, dirname):
		#config
		self.dirname = os.path.abspath(dirname)
		self.vfile = os.path.join(self.dirname, ".version")
		self.dbfile = os.path.join(self.dirname, "sonar.db")
		self.datadir = os.path.join(self.dirname, "data")
		self.seqdir = os.path.join(self.datadir, "seqs")
		self.__max_supported_prev_version = "0.0.9"
		self.check_dir()

	def check_dir(self):
		if not os.path.isdir(self.dirname) or len(os.listdir(self.dirname)) == 0:
			self.create_dir()
		elif not os.path.isdir(self.seqdir) or not os.path.isfile(self.vfile) or not os.path.isfile(self.dbfile):
			sys.exit("error: invalid sonar database directory")
		else:
			self.check_compatibility()

	def check_compatibility(self):
		with open(self.vfile, "r") as handle:
			dir_version = handle.read().strip()
		module_version = sonarDIR.get_version()
		if version.parse(dir_version) < version.parse(self.__max_supported_prev_version):
			sys.exit("version conflict: sonar database version " + dir_version + " is not supported anymore by this sonar version (module_version).")
		if version.parse(dir_version) > version.parse(module_version):
			sys.exit("version conflict: please update sonar (current version: " + module_version + ") to fit your database (version: " + dir_version + ").")

	def create_dir(self):
		os.makedirs(self.seqdir, exist_ok=True)
		open(self.dbfile, "w").close()
		with open(self.vfile, "w") as handle:
			handle.write(sonarDIR.get_version())

	def add_seq(self, seq, ref, paranoid=False):
		hashval = sonarDB.hash(seq)
		fname = os.path.join(self.seqdir, hashval + ".fna")
		if not os.isfile(fname):
			with open(fname, "w") as handle:
				handle.write(">" + hashval + "\n" + seq.strip().upper() + "\n" + ref)
		elif paranoid:
			with open(fname, "r") as handle:
				handle.readline()
				if seq.upper() != handle.readline.strip():
					sys.exit("error: hash collision for " + hashval)

	@staticmethod
	def get_version():
		with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), ".version"), "r") as handle:
			return handle.read().strip()

if __name__ == "__main__":
	print("sonarDIR", sonarDIR.get_version())
