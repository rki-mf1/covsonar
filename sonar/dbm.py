# DEPENDENCIES
import os
import re
import sys
import sqlite3
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from collections import OrderedDict, defaultdict
from urllib.parse import quote as urlquote
import itertools

# COMPATIBILITY
SUPPORTED_DB_VERSION = 4

# CLASS
class sonarDBManager():
	"""
	This object provides a sonarDB SQLite manager handler managing connections and
	providing context-safe transaction control.

	Notes
	-----
		This object should be called using a context manager to ensure rollbacks
		on abnormal termination.

	Example
	-------

	In this example the DOCTESTDB variable store the path to a database file

	>>> with sonarDBManager(DOCTESTDB) as dbm:
	... 	pass

	Parameters
	----------

	dbfile : str
		define a path to a non-existent or valid SONAR database file. If the
		file does not exist, a SONAR database is created.
	timeout : int [ -1 ]
		define busy timeout. Use -1 for no timeout.
	readonly : bool [ False ]
		define if the connection should be read-only
	debug : bool [ False ]
		debug mode (print selected sql queries)

	Attributes
	----------
	dbfile : str
		stores the path to the used SONAR database file.
	connection : object
		stores the SQLite3 connection
	cursor : method
		stores the SQLite3 cursor

	Dev Note
	--------
	A database row is returned as dictionary with column name as keys. Multiple
	rows are returned as list of those dictionaries.

	"""

	def __init__(self, dbfile, timeout=-1, readonly=False, debug=False, autocreate=False):
		self.connection = None
		if not autocreate and not os.path.isfile(dbfile):
			sys.exit("database error: database does not exists")
		self.dbfile = os.path.abspath(dbfile)
		self.cursor = None
		self.__timeout = timeout
		self.__mode = "ro" if readonly else "rwc"
		self.__uri = self.get_uri(dbfile)
		self.debug = debug
		self.__properties = False
		self.__illegal_properties = {'molecule', }
		self.__default_reference_id = False
		self.operators = { "genuine": {
								"=": "=",
								">": ">",
								"<": "<",
								">=": ">=",
								"<=": "<=",
								"IN": "IN",
								"LIKE": "LIKE",
								"BETWEEN": "BETWEEN"
								},
						   "negated": {
						   		"=": "!=",
								">": "<=",
								"<": ">=",
								">=": "<",
								"<=": ">",
								"IN": "NOT IN",
								"LIKE": "NOT LIKE",
								"BETWEEN": "NOT BETWEEN"
								}
							}


	def __enter__(self):
		self.connection, self.cursor = self.connect()
		if self.debug:
			self.connection.set_trace_callback(print)
		self.start_transaction()
		return self

	def __exit__(self, exc_type, exc_value, exc_traceback):
		if [exc_type, exc_value, exc_traceback].count(None) != 3:
			if self.__mode == "rwc":
				print("rollback database", file=sys.stderr)
				self.rollback()
		elif self.__mode == "rwc":
			self.commit()
		self.close()

	def __del__(self):
		if self.connection:
			self.close()

	def connect(self):
		con = sqlite3.connect(self.__uri + "?mode=" + self.__mode, self.__timeout, isolation_level = None, uri = True)
		con.row_factory = self.dict_factory
		cur = con.cursor()
		return con, cur

	def start_transaction(self):
		self.cursor.execute("BEGIN DEFERRED")

	def commit(self):
		self.connection.commit()

	def rollback(self):
		self.connection.rollback()

	def close(self):
		self.connection.close()

	def get_db_version(self):
		return self.cursor.execute('pragma user_version').fetchone()['user_version']

	def check_db_compatibility(self):
		dbver = self.get_db_version()
		if dbver != SUPPORTED_DB_VERSION:
			sys.exit("compatibility error: the given database is not compatible with this version of sonar (database version: " + str(dbver) + "; supported database version: " + str(SUPPORTED_DB_VERSION) +")")

	@property
	def default_reference_id(self):
		if self.__default_reference == False:
			sql = "SELECT accession FROM reference WHERE standard = 1 LIMIT 1"
			self.__default_reference = self.cursor.execute(sql).fetchone()
		return self.__default_reference

	@property
	def properties(self):
		if self.__properties == False:
			sql = "SELECT * FROM property"
			rows = self.cursor.execute(sql).fetchall()
			self.__properties = {} if not rows else { x['name']: x for x in rows }
			for pname in self.__properties:
				if self.__properties[pname]['standard'] is not None:
					if self.__properties[pname]['datatype'] == 'integer':
						self.__properties[pname]['standard'] = int(self.__properties['standard'])
				elif self.__properties[pname]['datatype'] == 'float':
					self.__properties[pname]['standard'] = float(self.__properties['standard'])
		return self.__properties

	# SETUP

	@staticmethod
	def get_uri(dbfile):
		return "file:" + urlquote(dbfile)

	@staticmethod
	def setup(filename, debug=False):
		with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), "db.sqlite"), 'r') as handle:
			sql = handle.read()
		uri = sonarDBManager.get_uri(filename)
		with sqlite3.connect(uri + "?mode=rwc", uri = True) as con:
			if debug:
				con.set_trace_callback(print)
			con.executescript(sql)

	def add_codon(self, translation_table, codon, aa):
		sql = "INSERT OR IGNORE INTO translation (id, codon, aa) VALUES(?, ?, ?);"
		self.cursor.execute(sql, [translation_table, codon, aa])

	def add_property(self, name, datatype, querytype, description, standard):
		if name in self.__illegal_properties:
			sys.exit("property error: " + name + " is not allowed to be used as property name.")
		sql = "INSERT INTO property (name, datatype, querytype, description, standard) VALUES(?, ?, ?, ?, ?);"
		self.cursor.execute(sql, [name, datatype, querytype, description, standard])

	def add_translation_table(self, translation_table):
		sql = "SELECT COUNT(*) FROM translation  WHERE id = ?;"
		if self.cursor.execute(sql, [translation_table]).fetchone()['COUNT(*)'] != 4096:
			for codon in itertools.product("ATGCRYSWKMBDHVN-", repeat = 3):
				codon = "".join(codon)
				try:
					aa = str(Seq.translate(codon, table=translation_table))
				except:
					aa = ""
				self.add_codon(translation_table, codon, aa)

	def add_reference(self, accession, description, organism, translation_table, standard=0):
		self.add_translation_table(translation_table)
		if standard:
			sql = "UPDATE reference SET standard = 0 WHERE standard != 0"
			self.cursor.execute(sql)
		sql = "INSERT INTO reference (id, accession, description, organism, translation_id, standard) VALUES(?, ?, ?, ?, ?, ?);"
		self.cursor.execute(sql, [None, accession, description, organism, translation_table, standard])
		return self.cursor.lastrowid


	# INSERTING DATA

	def insert_property(self, sample_id, property_name, property_value):
		sql = "INSERT OR REPLACE INTO sample2property (sample_id, property_id, value_" + self.properties[property_name]['datatype'] + ") VALUES(?, ?, ?);"
		self.cursor.execute(sql, [sample_id, self.properties[property_name]['id'], property_value])

	def insert_sample(self, sample_name):
		sql = "INSERT OR IGNORE INTO sample (name) VALUES(?);"
		self.cursor.execute(sql, [sample_name])
		sql = "SELECT id FROM sample WHERE name = ?"
		return self.cursor.execute(sql, [sample_name]).fetchone()['id']

	def insert_sequence(self, sample_id, sequence_hash):
		sql = "INSERT OR IGNORE INTO sample2sequence (sample_id, sequence_hash) VALUES(?, ?);"
		self.cursor.execute(sql, [sample_id, sequence_hash])
		sql = "SELECT id FROM sample2sequence WHERE sample_id = ? AND sequence_hash = ?"

	def insert_alignment(self, seqhash, element_id):
		sql = "INSERT OR IGNORE INTO alignment (id, sequence_hash, element_id) VALUES(?, ?, ?);"
		self.cursor.execute(sql, [None, seqhash, element_id])
		sql = "SELECT id FROM alignment WHERE element_id = ? AND sequence_hash = ?;"
		aid = self.cursor.execute(sql, [element_id, seqhash]).fetchone()['id']
		return aid

	def insert_molecule(self, reference_id, type, accession, symbol, description, segment, length, standard=0):
		if standard:
			sql = "UPDATE molecule SET standard = 0 WHERE reference_id = ? AND standard != 0"
			self.cursor.execute(sql, [reference_id])
		sql = "INSERT INTO molecule (id, reference_id, type, accession, symbol, description, segment, length, standard) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?);"
		self.cursor.execute(sql, [None, reference_id, type, accession, symbol, description, segment, length, standard])
		sql = "SELECT id FROM molecule WHERE accession = ?"
		mid = self.cursor.execute(sql, [accession]).fetchone()['id']
		return mid

	def insert_element(self, molecule_id, type, accession, symbol, description, start, end, strand, sequence, parent_id="", parts=None):
		sql = "INSERT INTO element (id, molecule_id, type, accession, symbol, description, start, end, strand, sequence) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?);"
		self.cursor.execute(sql, [None, molecule_id, type, accession, symbol, description, start, end, strand, sequence])
		sql = "SELECT id FROM element WHERE accession = ?"
		mid = self.cursor.execute(sql, [accession]).fetchone()['id']
		if parent_id != "":
			sql = "INSERT OR IGNORE INTO element2element (parent_id, child_id) VALUES(?, ?);"
			self.cursor.execute(sql, [parent_id, mid])
		if parts is not None:
			for part in parts:
				sql = "INSERT OR IGNORE INTO elempart (element_id, start, end, strand, segment) VALUES(?, ?, ?, ?, ?);"
				self.cursor.execute(sql, [mid] + part)
		return mid

	def insert_variant(self, alignment_id, ref, alt, start, end, parent_id=None):
		vid = self.get_variant_id(ref, alt, start, end, parent_id)
		if not vid:
			sql = "INSERT OR IGNORE INTO variant (id, ref, alt, start, end, parent_id) VALUES(?, ?, ?, ?, ?, ?);"
			self.cursor.execute(sql, [None, ref, alt, start, end, parent_id])
			vid = self.cursor.lastrowid
		sql = "INSERT OR IGNORE INTO alignment2variant (alignment_id, variant_id) VALUES(?, ?);"
		self.cursor.execute(sql, [alignment_id, vid])
		return vid

	# DELETING DATA

	def delete_sample(self, sample_name):
		sql = "DELETE FROM sample WHERE name = ?;"
		self.cursor.execute(sql, [sample_name])

	# SELECTING DATA

	def sample_exists(self, sample_name):
		sql = "SELECT EXISTS(SELECT 1 FROM sample WHERE name=? LIMIT 1) as found";
		return bool(self.cursor.execute(sql, [sample_name]).fetchone()['found'])

	def get_sample_id(self, sample_name):
		sql = "SELECT id FROM sample WHERE name = ? LIMIT 1;";
		row = self.cursor.execute(sql, [sample_name]).fetchone()
		return row['id'] if row else None

	def seq_exists(self, seqhash):
		sql = "SELECT EXISTS(SELECT 1 FROM sample2sequence WHERE sequence_id=? LIMIT 1) as found;"
		return bool(self.cursor.execute(sql, [seqhash]).fetchone()['found'])

	def get_alignment_id(self, seqhash, element_id):
		sql = "SELECT id FROM alignment WHERE sequence_hash=? AND element_id = ? LIMIT 1;"
		row = self.cursor.execute(sql, [seqhash, element_id]).fetchone()
		return None if row is None else row['id']

	def get_molecule_ids(self, reference_accession=None):
		if reference_accession:
			 condition = "reference.accession = ?"
			 val = reference_accession
		else:
			 condition = "reference.standard = ?"
			 val = 1
		sql = "SELECT molecule.accession, molecule.id FROM referenceView WHERE " + condition
		return {x['molecule.accession']: x[molecule.id] for x in self.cursor.execute(sql, val).fetchall() if x is not None }

	def get_molecule_id(self, molecule_accession):
		sql = "SELECT id FROM molecule WHERE accession = ?;"
		row = self.cursor.execute(sql, [molecule_accession]).fetchone()
		if row:
			row = row['id']
		return row

	def get_molecule_data(self, *fields, reference_accession=None):
		if not fields:
			fields = "*"
		elif "\"molecule.accession\"" not in fields:
			fields = list(fields)+ ["\"molecule.accession\""]
		if reference_accession:
			 condition = "\"reference.accession\" = ?"
			 vals = [reference_accession]
		else:
			 condition = "\"reference.standard\" = ?"
			 vals = [1]
		sql = "SELECT " + ", ".join(fields) + " FROM referenceView WHERE " + condition + ";"
		row = self.cursor.execute(sql, vals).fetchall()
		if row:
			return { x['molecule.accession']: x for x in row }
		return {}

	def get_elements(self, molecule_id, *types):
		sql = "SELECT * FROM element WHERE molecule_id = ?"
		if types:
			sql += " AND type IN (" + ", ".join(["?"] * len(types)) + ");"
		row = self.cursor.execute(sql, [molecule_id] + list(types)).fetchall()
		if not row:
			return []
		return row

	def get_source(self, molecule_id):
		return self.get_elements(molecule_id, 'source')[0]

	def get_default_molecule(self, reference_accession=None):
		if reference_accession:
			 condition = "\"reference.accession\" = ?"
			 val = reference_accession
		else:
			 condition = "\"reference.standard\" = ?"
			 val = 1
		sql = "SELECT \"molecule.accession\" FROM referenceView WHERE \"molecule.standard\" = ? AND " + condition
		return self.cursor.execute(sql, [1, val]).fetchone()

	def get_alignment_data(self, sample_id, element_id, *fields, limit=1):
		if not fields:
			fields = ['*']
		sql = "SELECT " + ", ".join(fields) + " FROM alignmentView WHERE \"sample.id\" = ? AND \"element.id\" = ? LIMIT "+ LIMIT +";"
		return self.cursor.execute(sql, [sample_id, element_id]).fetchall()

	def get_alignment_id(self, seqhash, element_id):
		sql = "SELECT id FROM alignment WHERE \"sequence_hash\" = ? AND \"element.id\" = ? LIMIT 1;"
		row = self.cursor.execute(sql, [seqhash, element_id]).fetchone()
		if row:
			return row['id']
		return None

	def get_variant_id(self, ref, alt, start, end, parent_id):
		sql = "SELECT id FROM variant WHERE ref = ? AND alt = ? AND start = ? AND end = ? AND parent_id = ?;"
		row = self.cursor.execute(sql, [ref, alt, start, end, parent_id]).fetchone()
		return None if row is None else row['id']

	def iter_dna_variants(self, sample_name, *element_ids):
		sql = "SELECT element_id, element.symbol, start, end, ref, alt, frameshift FROM variant WHERE name = ? AND element.id IN (" + ", ".join(['?'] * len(element_ids)) + ") AND type != 'protein' LEFT JOIN element ON element.id = element_id;"
		for row in self.cursor.execute(sql, [sample_name] + element_ids):
			yield row

	def iter_profile(self, sample_name, element_id):
		sql = "SELECT * FROM variantView WHERE \"sample.name\" = ? AND \"element.id\" = ? ORDER BY \"variant.start\";"
		for row in self.cursor.execute(sql, [sample_name, element_id]):
			yield row

	def iter_protein_variants(self, sample_name, *element_ids):
		sql = "SELECT element_id, element.symbol, start, end, ref, alt FROM variant WHERE name = ? AND element.id IN (" + ", ".join(['?'] * len(element_ids)) + ") AND type = 'protein' LEFT JOIN element ON element.id = element_id;"
		for row in self.cursor.execute(sql, [sample_name] + element_ids):
			yield row

	def get_frameshift_vars(self, sample_name, *element_ids):
		sql = "SELECT element_id, element.symbol, start, end, ref, alt FROM variant WHERE name = ? AND element.id IN (" + ", ".join(['?'] * len(element_ids)) + ") AND frameshift = 1 LEFT JOIN element ON element.id = element_id;"
		for row in self.cursor.execute(sql, [sample_name] + element_ids):
			yield row

	def count_samples(self):
		sql = "SELECT COUNT(*) FROM sample;"
		return self.cursor.execute(sql).fetchone()['COUNT(*)']

	def count_sequences(self):
		sql = "SELECT COUNT(DISTINCT seqhash) FROM sample2sequence;"
		return self.cursor.execute(sql).fetchone()['COUNT(seqhash)']

	def count_property(self, property_name, distinct=False, not_null=True):
		d = "DISTINCT " if distinct else ""
		n = " WHERE " + self.properties[property_name] + " IS NOT NULL"
		sql = "SELECT COUNT(" + d + self.properties[property_name] + ") as count FROM sample2property" + n + ";"
		return self.cursor.execute(sql).fetchone()['count']

	def get_properties(self, sample_name, property_name=None):
		if property_name:
			sql = "SELECT * FROM propertyView WHERE \"sample.name\" = ? AND property_name = ?;"
			vals = [sample_name, property_name]
		else:
			sql = "SELECT * FROM propertyView WHERE \"sample.name\" = ?"
			vals = [sample_name]
		properties = {}
		for row in self.cursor.execute(sql, vals):
			if row == None:
				return {}
			properties[row['property.name']] = row["sample2property.value_" + row['property.datatype']]
		return properties

	def get_seqhashes(self, sample_name):
		sql = "SELECT \"sequence.hash\" FROM sequenceView WHERE sequence.name = ?"
		return [ x[sequence.hash] for x in self.cursor.execute(sql, [sample_name]).fetchall() if x is not None ]

	def get_sequence(self, sample_name):
		sql = "SELECT \"sequence.hash\" FROM sequenceView WHERE sequence.name = ?"
		return [ x[sequence.hash] for x in self.cursor.execute(sql, [sample_name]).fetchall() if x is not None ]


	def get_earliest_import(self):
		sql = "SELECT MIN(imported) as import FROM genome WHERE import IS NOT NULL;"
		return self.cursor.execute(sql).fetchone()['import']

	def get_latest_import(self):
		sql = "SELECT MAX(imported) as import FROM genome WHERE import IS NOT NULL;"
		return self.cursor.execute(sql).fetchone()['import']

	def get_earliest_date(self):
		sql = "SELECT MIN(date) as date FROM genome WHERE date IS NOT NULL;"
		return self.cursor.execute(sql).fetchone()['date']

	def get_latest_date(self):
		sql = "SELECT MAX(date) as date FROM genome WHERE date IS NOT NULL;"
		return self.cursor.execute(sql).fetchone()['date']

	def get_element_parts(self, element_id=None):
		sql = "SELECT start, end, strand FROM elempart WHERE element_id = ? ORDER BY segment"
		return self.cursor.execute(sql, [element_id]).fetchall()

	def get_sequence(self, element_id=None):
		sql = "SELECT sequence, type FROM element WHERE id = ?"
		row = self.cursor.execute(sql, [element_id]).fetchone()
		return None if row is None else row['sequence']

	def extract_sequence(self, element_id=None, translation_table=None):
		sql = "SELECT sequence, type, id FROM element WHERE id = ?"
		row = self.cursor.execute(sql, [element_id]).fetchone()
		if not row:
			return None
		element_id = row['id']
		while row and row['type'] != "source":
			sql = "SELECT parent_id FROM element2element WHERE child_id = ?"
			row = self.cursor.execute(sql, [row['id']]).fetchone()
			sql = "SELECT sequence, type, id FROM element WHERE id = ?"
			row = self.cursor.execute(sql, [row['parent_id']]).fetchone()
		sequence = row['sequence']
		parts = []
		for part in self.get_element_parts(element_id):
			parts.append(FeatureLocation(part['start'], part['end'], strand=part['strand']))
		feat = CompoundLocation(parts) if len(parts) > 1 else parts[0]
		if translation_table is None:
			return str(feat.extract(sequence))
		return str(Seq(feat.extract(sequence)).translate(table=translation_table, stop_symbol=""))

	# MATCHING PROFILES
	def get_operator(self, val):
		if val.startswith("^"):
			return val[1:], self.operators["negated"]
		return val, self.operators["genuine"]

	def get_conditional_expr(self, field, operator, *vals):
		condlist = []
		vallist = []

		#transforming operators
		if operator == "=" and len(vals) > 1:
			operator = "IN"
		elif operator == "!=" and len(vals) > 1:
			operator = "NOT IN"

		#creating conditions
		if operator == "IN" or operator == "NOT IN":
			condlist.append(field + " " + operator + "(" + ", ".join(['?'] * len(vals)) + ")" )
			vallist.extend(vals)

		elif operator == "BETWEEN" or operator == "NOT BETWEEN":
			for val in vals:
				condlist.append(field + " " + operator + " ? AND ?" )
				vallist.extend(val)

		else:
			condlist.extend([field + " " + operator + " ?"] * len(vals))
			vallist.extend(vals)

		return condlist, vallist

	def query_numeric(self, field, *vals, link="OR"):
		link = " " + link.strip() + " "
		defaultop = "="
		op1 = re.compile(r'^(\^*)((?:>|>=|<|<=|!=|=)?)(-?[1-9]+[0-9]*)$')
		op2 = re.compile(r'^(\^*)(-?[1-9]+[0-9]*):(-?[1-9]+[0-9]*)$')
		errmsg =  "query error: numeric value or range expected for field " + field + "(got: "
		data = defaultdict(set)

		for val in vals:
			val = str(val).strip()
			# single value
			if not ":" in val:
				match = op1.match(val)
				if not match:
					sys.exit(errmsg + val + ")")
				op = match.group(2) if match.group(2) else defaultop
				operator = self.operators['negated'][op] if match.group(1) else self.operators['genuine'][op]
				val = match.group(3)
				data[operator].add(val)

			# range
			else:
				match = op2.match(val)
				if not match:
					sys.exit(errmsg + val + ")")
				operator = self.operators['negated']["BETWEEN"] if match.group(1) else self.operators['genuine']["BETWEEN"]
				val = (match.group(2), match.group(3))
				data[operator].add(val)

		for operator, valset in data.items():
			condition, vals = self.get_conditional_expr(field, operator, *valset)
			conditions.extend(condition)
			vallist.extend(vals)

		return link.join(conditions), vallist

	def query_float(self, field, *vals, link="OR"):
		link = " " + link.strip() + " "
		defaultop = "="
		op1 = re.compile(r'^(\^*)((?:>|>=|<|<=|!=|=)?)(-?[1-9]+[0-9]*(?:.[0-9]+)*)$')
		op2 = re.compile(r'^(\^*)(-?[1-9]+[0-9]*(?:.[0-9]+)*):(-?[1-9]+[0-9]*(?:.[0-9]+)*)$')
		errmsg =  "query error: decimal value or range expected for field " + field + "(got: "
		data = defaultdict(set)

		for val in vals:
			val = str(val).strip()
			# single value
			if not ":" in val:
				match = op1.match(val)
				if not match:
					sys.exit(errmsg + val + ")")
				op = match.group(2) if match.group(2) else defaultop
				operator = self.operators['negated'][op] if match.group(1) else self.operators['genuine'][op]
				val = match.group(3)
				data[operator].add(val)

			# range
			else:
				match = op2.match(val)
				if not match:
					sys.exit(errmsg + val + ")")
				operator = self.operators['negated']["BETWEEN"] if match.group(1) else self.operators['genuine']["BETWEEN"]
				val = (match.group(2), match.group(3))
				data[operator].add(val)

		for operator, valset in data.items():
			condition, vals = self.get_conditional_expr(field, operator, *valset)
			conditions.extend(condition)
			vallist.extend(vals)

		return link.join(conditions), vallist

	def query_dates(self, field, *vals, link="OR"):
		link = " " + link.strip() + " "
		defaultop = "="
		op1 = re.compile(r'^(\^*)((?:>|>=|<|<=|!=|=)?)([0-9]{4}-[0-9]{2}-[0-9]{2})$')
		op2 = re.compile(r'^(\^*)([0-9]{4}-[0-9]{2}-[0-9]{2}):([0-9]{4}-[0-9]{2}-[0-9]{2})$')
		errmsg =  "query error: date or date range expected for field " + field + "(got: "
		data = defaultdict(set)

		for val in vals:
			val = str(val).strip()
			# single date
			if not ":" in val:
				match = op1.match(val)
				if not match:
					sys.exit(errmsg + val + ")")
				op = match.group(2) if match.group(2) else defaultop
				operator = self.operators['negated'][op] if match.group(1) else self.operators['genuine'][op]
				val = match.group(3)
				data[operator].add(val)

			# date range
			else:
				match = op2.match(val)
				if not match:
					sys.exit(errmsg + val + ")")
				operator = self.operators['negated']["BETWEEN"] if match.group(1) else self.operators['genuine']["BETWEEN"]
				val = (match.group(2), match.group(3))
				data[operator].add(val)

		for operator, valset in data.items():
			condition, vals = self.get_conditional_expr(field, operator, *valset)
			conditions.extend(condition)
			vallist.extend(vals)

		return link.join(conditions), vallist

	def query_string(self, field, *vals, link="OR"):
		link = " " + link.strip() + " "
		data = defaultdict(set)

		for val in vals:
			if val.startswith("^"):
				val = val[1:]
				opkey = 'negated'
			else:
				opkey = 'genuine'
			if val.startswith("%") or val.endswith("%"):
				operator = self.operators[opkey]["LIKE"]
			else:
				operator = self.operators[opkey]["="]

			data[operator].add(val)

		for operator, valset in data.items():
			condition, vals = self.get_conditional_expr(field, operator, *valset)
			conditions.extend(condition)
			vallist.extend(vals)

		return link.join(conditions), vallist

	def query_zip(self, field, *vals, link="OR"):
		link = " " + link.strip() + " "
		data = defaultdict(set)

		for val in vals:
			if val.startswith("^"):
				val = val[1:]
				opkey = 'negated'
			else:
				opkey = 'genuine'
			val = val.strip("%") + "%"
			operator = self.operators[opkey]["LIKE"]

			data[operator].add(val)

		for operator, valset in data.items():
			condition, vals = self.get_conditional_expr(field, operator, *valset)
			conditions.extend(condition)
			vallist.extend(vals)

		return link.join(conditions), vallist

	def query_metadata(self, name, *vals):
		conditions = []
		valueList = []

		# query dates
		if self.properties[name]['querytype'] == "date":
			a,b = self.query_dates(name, vals)
			conditions.append(a)
			valueList.append(b)

		# query numeric
		elif self.properties[name]['querytype'] == "numeric":
			a,b = self.query_numeric(name, vals)
			conditions.append(a)
			valueList.append(b)

		# query text
		elif self.properties[name]['querytype'] == "string":
			a,b = self.query_string(name, vals)
			conditions.append(a)
			valueList.append(b)

		# query zip
		elif self.properties[name]['querytype'] == "zip":
			a,b = self.query_zip(name, vals)
			conditions.append(a)
			valueList.append(b)

		return "SELECT sample.id FROM metadataView WHERE " + " AND ".join(conditions), valueList

	def query_profile(*profiles, reference_accession=None, molecule_symbol=None):
		profile_regex = re.compile(r'^(|[^:]+:)?([^:]+:)?(del:\*?[0-9]+(?:-[0-9]+)?\*?|[A-Z]+[0-9]+[A-Z]+)$')
		snv_regex = re.compile(r'^([A-Z]+)([0-9]+)([A-Z]+)$')
		conditions = []
		valueList = []

		# set reference
		if reference_accession is None:
			conditions.append("reference.standard = ?")
			valueList.append(1)
		else:
			conditions.append("reference.accession = ?")
			valueList.append(reference_accession)

		# set variants
		for profile in filter(None, profile.split(" ")):
			elems = profiles.split(":")
			condition = []

			match = profile_regex.match(profile)

			if not match:
				sys.exit("error: " + profile + " is not a valid variant definition.")

			## set molecule
			if match.group(1):
				condition.append("molecule.symbol = ?")
				valueList.append(match.group(1)[:-1])
			else:
				condition.append("molecule.standard = ?")
				valueList.append(1)

			## set element
			if match.group(2):
				condition.append("element.type = ?")
				valueList.append("protein")
				condition.append("element.symbol = ?")
				valueList.append(match.group(2)[:-1])
			else:
				condition.append("element.type = ?")
				valueList.append("source")

			## set variant:
			### snv
			if not match.group(3).startswith("del:"):
				mut = snv_regex.match()
				condition.append("variant.ref = ?")
				condition.append("variant.start = ?")
				condition.append("variant.end = ?")
				condition.append("variant.alt = ?")
				valueList.extend([mut.group(1), int(mut.group(2))-1, int(mut.group(2)), mut.group(3)])

			### del
			else:
				coords = match.group(3)[4:]
				xop = "<=" if coords.startswith("*") else "="
				yop = ">=" if coords.endswith("*") else "="
				coords = [int(x) for x in coords.strip("*").split("-")]
				x = coords[0]-1
				y = coords[1] if len(coords) == 2 else coords[0]
				condition.append(["variant.start = ? AND variant.end = ?"])
				valueList.extend([x, y])

			# deletion handling
			if elems[-1].startswith("del"):
				coords = elems[-1][3:].split("-")
				start = int(coords[0].strip("*"))
			if coords[0][0] == "*" and start > 0:
				condition.append("variant.start=" + str(start-1) + "variant.alt != ''")
			if len(coords) == 2:
				end = int(coords[1].strip("*"))
				if coords[0][0] == "*" and start > 0:
					pass

		return "SELECT sample.id FROM filterView WHERE " + " AND ".join(conditions), valueList

	def match(metadata, *profiles, reference=None):

		#sqls for metadata-based filtering
		msqls = []
		mvals = []
		for mname, vals in metadata.items():
			msql, mval = self.query_metadata(mname, *vals)
			msqls.append(msql)
			mvals.extend(mval)

		#sqls for genomic profile based filtering
		psqls = []
		pvals = []
		for profile in profiles.items():
			psql, pval = psqls.extend(self.query_profile(profile, reference=reference))
			psqls.append(psql)
			pvals.extend(pval)


	# UPDATE DATA


	# MISC

	@staticmethod
	def optimize(dbfile):
		with sqlite3.connect(dbfile) as con:
			con.executescript("VACUUM")

	@staticmethod
	def dict_factory(cursor, row):
		d = OrderedDict()
		for idx, col in enumerate(cursor.description):
			d[col[0]] = row[idx]
		return d
