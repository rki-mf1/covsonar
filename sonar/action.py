#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

# DEPENDENCIES
import os
import re
import sys
from Bio.SeqUtils.CheckSum import seguid
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Emboss.Applications import StretcherCommandline
import pickle
from tqdm import tqdm
import tempfile
import itertools
from contextlib import ExitStack
from sonar.dbm import sonarDBManager
from sonar.align import sonarAligner

# CLASS
class sonarActions(object):
	"""
	this object provides sonarActions functionalities and intelligence

	Notes
	-----
		Please note, that genomic and protein coordinates are expected to be  and
		returned 0-based by this object, except for formatted profiles.
		While start or single coordinates are inclusive, end coordinates of
		ranges are exclusive, expressed in a mathematical notation: [start, end).
		Only in formatted profiles start and end coordinates are 1-based and both
		inclusive.

	Examples
	--------

	In this example the path to the database is stored in DOCTESTDB.

	>>> db = sonarActions(DOCTESTDB)

	Parameters
	----------
	dbfile : str
		define a path to a non-existent or valid SONAR database file. If the
		file does not exist, a SONAR database is created.
	translation_table : int
		define the genetic code table used for in silico translation (see
		https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) [ 1 ]

	Attributes
	----------
	db : str
		stores the absolute path to the used SONAR database file
	reffna : str
		stores the absolute path to the built-in FASTA file containing the reference
		genome sequence
	refgff : str
		stores the absolute path to the built-in GFF3 file containing the reference
		genome annotation
	translation_table : int
		stores the genetic code table used for in silico translation (see
		https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) [ 1 ]
	refseq : str
		stores the upper-case sequence of the built-in reference genome
	refdescr : str
		stores the FASTA header of the built-in reference genome
	refgffObj : object
		stores the sonarGFF object based on the built-in reference genome
		annotation
	iupac_nt_code : dict
		stores a dict with IUPAC one-letter nucleotide codes as keys and the
		respective set of matching explicit IUPAC one-letter nucleotide codes
		as values (e.g {"W": set('A', 'T')})
	iupac_explicit_nt_code : dict
		stores a set containing all non-ambiguous IUPAC one-letter nucleotide codes
	iupac_ambig_nt_code : set
		stores a set containing all ambiguous IUPAC one-letter nucleotide codes
	iupac_aa_code : dict
		stores a dict with IUPAC one-letter amino acid codes as keys and
		the respective set of matching IUPAC one-letter amino acids codes as values
	iupac_explicit_aa_code : dict
		stores a set containing all non-ambiguous IUPAC one-letter amino acid codes
	iupac_ambig_aa_code : dict
		stores a set containing all ambiguous IUPAC one-letter amino acid codes
	dna_var_regex : compiled re expression
		stores a compiled re expression that matches to nucleotide profiles but
		not to amino acid profiles
	aa_var_regex : compiled re expression
		stores a compiled re expression that matches to amino acid profiles but
		not to nucleotide profiles
	del_regex : compiled re expression
		stores a compiled re expression that matches to deletion profiles on
		nucleotide as well as on amino acid level.
	dnavar_grep_regex : compiled re expression
		stores a compiled re expression that matches to snp or dna insertion
		profiles with eference allele, genomic position and variant allele
		as groups.
	codedict : dict
		stores a dictionary with "dna" and "aa" containing the field name in the
		database that stores the profile data, the one letter code with and
		without ambiguities
	"""
	def __init__(self, dbfile, translation_table = 1, debug=False):
		self.db = os.path.abspath(dbfile)
		self.debug = debug
		self.__moduledir = self.get_module_file()
		self.reffna = os.path.join(self.__moduledir, "ref.fna")
		self.refgff = os.path.join(self.__moduledir, "ref.gff3")
		self.translation_table = translation_table
		self.__refseq = None
		self.__refdescr = None
		self.__refgffObj = None
		self.__iupac_nt_code = None
		self.__iupac_aa_code = None
		self.__iupac_explicit_nt_code = None
		self.__iupac_explicit_aa_code = None
		self.__iupac_ambig_nt_code = None
		self.__iupac_ambig_aa_code = None
		self.__terminal_letters_regex = re.compile("[A-Z]$")
		self.__dna_var_regex = None
		self.__aa_var_regex = None
		self.__del_regex = None
		self.__dnavar_grep_regex = None
		self.__codedict = None
		self.__fasta_tag_regex = None

	# DB MAINTENANCE

	@staticmethod
	def get_module_file(*join_with):
		return os.path.join(os.path.dirname(os.path.realpath(__file__)), *join_with)

	@staticmethod
	def setup_db(fname, default_setup=True, debug=False):
		if os.path.isfile(fname):
			sys.exit("setup error: " + fname + " does already exist.")
		sonarDBManager.setup(fname, debug=debug)
		if default_setup:
			with sonarDBManager(fname, debug=debug) as dbm:
				## adding pre-defined sample properties
				dbm.add_property("accession", "text", "text", "sample accession", "")
				dbm.add_property("imported", "date", "date", "date sample has been imported to the database", "")
				dbm.add_property("modified", "date", "date", "date when sample data has been modified lastly", "")

				## adding built-in reference (unsegmented genome)
				records = [x for x in sonarActions.iter_genbank(sonarActions.get_module_file("ref.gb"))]
				ref_id = dbm.add_reference(records[0]['accession'], records[0]['description'], records[0]['organism'], 1, 1)

				## adding reference molecule and elements
				for i, record in enumerate(records):
					gene_ids = {}
					s = 1 if i == 0 else 0

					mol_id = dbm.insert_molecule(ref_id, record['moltype'], record['accession'], record['symbol'], record['description'], i, record['length'], s)

					### source handling
					source_id = dbm.insert_element(mol_id, "source", record['source']['accession'], record['source']['symbol'], record['source']['description'], record['source']['start'], record['source']['end'], record['source']['strand'], record['source']['sequence'], parts=record['source']['parts'])
					if record['source']['sequence'] != dbm.get_sequence(source_id):
						sys.exit("genbank error: could not recover sequence of '" + record['source']['accession'] + "' (source)")

					### gene handling
					for elem in record['gene']:
						gene_ids[elem['accession']] = dbm.insert_element(mol_id, "gene", elem['accession'], elem['symbol'], elem['description'], elem['start'], elem['end'], elem['strand'], elem['sequence'], parent_id=source_id, parts=elem['parts'])
						if elem['sequence'] != dbm.extract_sequence(gene_ids[elem['accession']]):
							sys.exit("genbank error: could not recover sequence of '" + elem['accession'] + "' (gene)")

					### cds handling
					for elem in record['cds']:
						cid = dbm.insert_element(mol_id, "cds", elem['accession'], elem['symbol'], elem['description'], elem['start'], elem['end'], "", elem['sequence'], gene_ids[elem['gene']], elem['parts'])
						if elem['sequence'] != dbm.extract_sequence(cid, translation_table=1):
							sys.exit("genbank error: could not recover sequence of '" + elem['accession'] + "' (cds)")
	# DATA IMPORT
	@property
	def fasta_tag_regex(self):
		if self.__fasta_tag_regex == None:
			tags = set("molecule")
			with sonarDBManager(self.db, readonly=True, debug=self.debug) as dbm:
				if dbm.properties:
					tags.update(dbm.properties.keys())
			self.__fasta_tag_regex = re.compile("\[(" + "|".join(tags) + ")=([^\[\]=]+)\]")
		return self.__fasta_tag_regex

	## fasta and genbank handling handling
	@staticmethod
	def iter_genbank(fname):
		gb_data = {}
		for gb_record in SeqIO.parse(fname, "genbank"):

			## adding general annotation
			gb_data['accession'] = gb_record.name + "." + str(gb_record.annotations['sequence_version'])
			gb_data['symbol'] = gb_record.annotations['symbol'] if 'symbol' in gb_record.annotations else ""
			gb_data['organism'] = gb_record.annotations['organism']
			gb_data['moltype'] = None
			gb_data['description'] = gb_record.description
			gb_data['length'] = None
			gb_data['segment'] = ""
			gb_data['gene'] = []
			gb_data['cds'] = []
			gb_data['source'] = None

			## adding source annotation
			source = [x for x in gb_record.features if x.type == "source"]
			if len(source) !=  1:
				sys.exit("genbank error: expecting exactly one source feature (got: " + str(len(source)) + ")")
			feat = source[0]
			gb_data['moltype'] = feat.qualifiers["mol_type"][0] if "mol_type" in feat.qualifiers else ""
			gb_data['source'] = {
				"accession": gb_data['accession'],
				"symbol": gb_data['accession'],
				"start": int(feat.location.start),
				"end": int(feat.location.end),
				"strand": "",
				"sequence": sonarActions.harmonize(feat.extract(gb_record.seq)),
				"description": "",
				"parts": [[int(x.start), int(x.end), x.strand, i] for i, x in enumerate(feat.location.parts, 1)]
			}
			gb_data['length'] = len(gb_data['source']['sequence'])
			if "segment" in feat.qualifiers:
				gb_data['segment'] = feat.qualifiers["segment"][0]

			for feat in gb_record.features:
				## adding gene annotation
				if feat.type == "gene":
					gb_data['gene'].append({
						"accession": feat.id if feat.id != "<unknown id>" else feat.qualifiers['gene'][0],
						"symbol": feat.qualifiers['gene'][0],
						"start": int(feat.location.start),
						"end": int(feat.location.end),
						"strand": feat.strand,
						"sequence": sonarActions.harmonize(feat.extract(gb_data['source']['sequence'])),
						"description": "",
						"parts": [[int(x.start), int(x.end), x.strand, i] for i, x in enumerate(feat.location.parts, 1)]
					})
				## adding cds annotation
				elif feat.type == "CDS":
					gb_data['cds'].append({
						"accession": feat.qualifiers['protein_id'][0],
						"symbol": feat.qualifiers['gene'][0],
						"start": int(feat.location.start),
						"end": int(feat.location.end),
						"strand": feat.strand,
						"gene": feat.qualifiers['gene'][0],
						"sequence": feat.qualifiers['translation'][0],
						"description": feat.qualifiers['product'][0],
						"parts": [[int(x.start), int(x.end), x.strand, i] for i, x in enumerate(feat.location.parts, 1)]
					})
			yield gb_data

	@staticmethod
	def iter_fasta(fname):
		seq = []
		header = None
		with open(fname) as handle:
			for line in handle:
				line = line.strip()
				if line.startswith(">"):
					if header and seq:
						yield id, header, sonarActions.harmonize("".join(seq))
					header = line[1:]
					id = header.split()[0]
					seq = []
				elif line != "":
					seq.append(line)
			if header and seq:
				yield id, header, sonarActions.harmonize("".join(seq))

	## seq handling
	@staticmethod
	def hash(seq):
		"""
		static function to hash any sequence using SEGUID (SHA-1 hash of the
		upper-case sequence)

		Parameters
		----------
		seq : str
			define a sequence to hash

		Returns
		-------
		str
			seguid

		"""
		return seguid(seq)

	@staticmethod
	def harmonize(seq):
		"""
		static function to return a sequence in upper case format and with T instead of U

		Parameters
		----------
		seq : str
			define a sequence to harmonize

		Returns
		-------
		str
			sequence

		"""
		return str(seq).strip().upper().replace("U", "T")

	## alignment handling

	## adding sample
	def add_sample(self, name, seguid, property_dict, dbm=None):
		with ExitStack() as stack:
			if dbm is None:
				dbm = stack.enter_context(sonarDBManager(self.db, debug=self.debug))
			return dbm.insert_sample(name, seguid, **property_dict)

	def add_samples_from_fasta(self, *fnames, reference_accession=None, dbm=None, allow_update=True, ignore_errors=True, disable_progress=False):
		with ExitStack() as stack:
			if dbm is None:
				dbm = stack.enter_context(sonarDBManager(self.db, debug=self.debug))

			# get molecule and annotation data from database
			default_molecule = dbm.get_default_molecule_accession(reference_accession=reference_accession)
			molecules = dbm.get_molecule_data(reference_accession=reference_accession)
			sources = {x['reference.accession']: dbm.get_source(x['reference.id']) for x in molecules.values()}

			#prepare progressbar
			total_entries = 0
			if not disable_progress:
				for fname in fnames:
					with open(fname, "r") as handle:
						for line in handle:
							if line.startswith(">"):
								total_entries += 1

			# read fasta
			reference_files = {}
			try:
				with tqdm(total=total_entries, desc="importing data ...", disable=disable_progress) as pbar:
					for fname in fnames:
						for sample_name, header, sequence in self.iter_fasta(fname):
							pbar.update(1)

							## extract molecule accession and sample properties from fasta header
							properties = {}
							for match in self.fasta_tag_regex.finditer(header):
								properties[match.group(1)] = match.group(2)

							molecule_accession = molecules[molecules[properties]['molecule']['id']] if "molecule" in properties else default_molecule
							molecule_id = molecules[molecule_accession]['molecule.id']
							molecule_len = molecules[molecule_accession]['molecule.length']

							try:
								source_id = sources[molecule_accession]['id']
							except:
								sys.exit("reference error: no molecule '" + molecule + "' registered for reference '" + reference_accession + "'.")
							molecule_seq = sources[molecule_accession]['sequence']

							## check sample data
							sample_id = dbm.get_sample_id(sample_name)
							seguid = sonarActions.hash(sequence)

							### handle sample updates
							if sample_id:
								#### check provided sequence (molecule specific)
								stored_seqhash = dbm.get_alignment_data(sample_id, source_id, "\"sequence.hash\"")['sequence.hash']
								if stored_seqhash != seguid:
									if ignore_errors:
										seguid = stored_seqhash
									elif not allow_update:
										sys.exit("import error: cannot update sequence of molecule " + molecule_sample + " from sample " + sample_name)

								#### check sample properties
								stored_properties = dbm.get_properties(sample_name)
								for pname, pval in properties.items():
									if stored_properties[pname] != pval:
										if ignore_errors:
											properties[pname] = stored_properties[pname]
										elif not allow_update:
											sys.exit("import error: cannot update property '" + pname + "' of sample " + sample_name)

							## set property default values
							for pname in [x for x in dbm.properties if x not in properties]:
								properties[pname] = dbm.properties[pname]['standard']

							## insert data
							sid = dbm.insert_sample(sample_name)
							dbm.insert_sequence(sid, seguid)
							for pname, pval in properties.items():
								dbm.insert_property(sid, pname, pval)

							## align if needed and insert variants
							aid = dbm.get_alignment_id(seguid, source_id)
							if aid is None:
								aid = dbm.insert_alignment(seguid, source_id)
								aligner = sonarAligner()
								result = aligner.iter_diffs(molecule_seq, sequence, cpus=3)
								#print(list(result))
								for variant in result:
									dbm.insert_variant(aid, variant[0], variant[2], variant[1], variant[1] + max(len(variant[0]), len(variant[2])))

							## paranoia check (checking nucleotide level variants)
							self.paranoia_test(sample_name, source_id, molecule_seq, sequence, molecule_len, dbm=dbm)
			finally:
				for fname in reference_files.values():
					os.remove(fname)


	def paranoia_test(self, sample_name, source_id, reference_sequence, target_sequence, ref_len, dbm=None):
		errmsg = None
		with ExitStack() as stack:
			if dbm is None:
				dbm = stack.enter_context(sonarDBManager(self.db, debug=self.debug))
			restored_sequence = list(reference_sequence)
			restored_sequence_prefix = ""
			for v in dbm.iter_profile(sample_name, source_id):
				#print(v)
				#gap handling
				if v['variant.alt'] == "":
					for j, i in enumerate(range(v['variant.start'], v['variant.end'])):
						if reference_sequence[i] != v['variant.ref'][j]:
							errmsg = "paranoia alarm: failed reference deduction for sample '" + sample_name +"' (expected "+ reference_sequence[i] +" got " + v['variant.ref'] + " at reference position " + str(i+1) + ")"
						restored_sequence[i] = ""
				#insert handling
				elif v['variant.ref'] == "":
					pos = v['variant.start']
					if pos == -1:
						restored_sequence_prefix = v['variant.alt']
					else:
						restored_sequence[pos] += v['variant.alt']
				#snp handling
				else:
					if reference_sequence[v['variant.start']] != v['variant.ref']:
						errmsg = "paranoia alarm: mismatch in reference sequence (expected " + reference_sequence[v['variant.start']] + " got " + v['variant.ref'] + " at position " + str(v['variant.start']) + ")"
					restored_sequence[v['variant.start']] = v['variant.alt']

			restored_sequence = restored_sequence_prefix + "".join(restored_sequence)
			if restored_sequence.strip("N.") != target_sequence.strip("N."):
				errmsg = "paranoia alarm: sequence restore failed for " + sample_name

			if errmsg:
				with open("paranoia.fna", "w") as handle:
					handle.write(">original " + sample_name + "\n" + target_sequence + "\n>recovered\n" + "".join(restored_sequence) +"\n>ref\n" + reference_sequence)
				with open("paranoia.log", "w") as handle:
					i = 0
					handle.write("0-based pos\texpected\tgot\n")
					for nucs in zip(target_sequence.strip("N."), "".join(restored_sequence).strip("N.")):
						if nucs[0] != nucs[1]:
							handle.write(str(i) + "\t" + "\t".join(nucs) + "\n")
						i += 1
				sys.exit(errmsg)
