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
		if self.debug:
			con.set_trace_callback(print)
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
			#for pname in self.__properties:
			#	if self.__properties[pname]['standard'] is not None:
			#		if self.__properties[pname]['datatype'] == 'integer':
			#			self.__properties[pname]['standard'] = int(self.__properties['standard'])
			#		elif self.__properties[pname]['datatype'] == 'float':
			#			self.__properties[pname]['standard'] = float(self.__properties['standard'])
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

	def add_property(self, name, datatype, querytype, description, standard=None):
		if not re.match("^[a-zA-Z0-9_]+$", name):
			sys.exit("error: invalid property name (property names can contain only letters, numbers and underscores)")
		if name in self.properties:
			sys.exit("error: a property named " + name + " already exists in the given database.")
		sql = "INSERT INTO property (name, datatype, querytype, description, standard) VALUES(?, ?, ?, ?, ?);"
		self.cursor.execute(sql, [name, datatype, querytype, description, standard])
		self.__properties = False
		pid = self.properties[name]['id']
		if not standard is None:
			sql = "INSERT INTO sample2property (property_id, value_" + self.properties[name]['datatype'] + ", sample_id) SELECT ?, ?, id FROM sample WHERE 1"
			vals = [pid, standard]
			self.cursor.execute(sql, vals)
		return pid

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

	def insert_sequence(self, seqhash):
		sql = "INSERT OR IGNORE INTO sequence (seqhash) VALUES(?);"
		self.cursor.execute(sql, [seqhash])
		return seqhash

	def insert_sample(self, sample_name, seqhash):
		self.insert_sequence(seqhash)
		sql = "INSERT OR IGNORE INTO sample (name, seqhash) VALUES(?, ?);"
		sid = self.cursor.execute(sql, [sample_name, seqhash])
		sql = "SELECT id FROM sample WHERE name = ?"
		sid = self.cursor.execute(sql, [sample_name]).fetchone()['id']
		for pname in self.properties:
			if not self.properties[pname]['standard'] is None:
				self.insert_property(sid, pname, self.properties[pname]['standard'])
		return

	def insert_alignment(self, seqhash, element_id):
		sql = "INSERT OR IGNORE INTO alignment (id, seqhash, element_id) VALUES(?, ?, ?);"
		self.cursor.execute(sql, [None, seqhash, element_id])
		sql = "SELECT id FROM alignment WHERE element_id = ? AND seqhash = ?;"
		aid = self.cursor.execute(sql, [element_id, seqhash]).fetchone()['id']
		return aid

	def insert_molecule(self, reference_id, type, accession, symbol, description, segment, length, standard=0):
		if symbol.strip() == "":
			symbol = accession
		if standard:
			sql = "UPDATE molecule SET standard = ? WHERE reference_id = ? AND standard != 0"
			self.cursor.execute(sql, [0, reference_id])
		sql = "INSERT INTO molecule (id, reference_id, type, accession, symbol, description, segment, length, standard) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?);"
		self.cursor.execute(sql, [None, reference_id, type, accession, symbol, description, segment, length, standard])
		sql = "SELECT id FROM molecule WHERE accession = ?"
		mid = self.cursor.execute(sql, [accession]).fetchone()['id']
		return mid

	def insert_element(self, molecule_id, type, accession, symbol, description, start, end, strand, sequence, standard=0, parent_id="", parts=None):
		if symbol.strip() == "":
			symbol = accession
		if standard:
			sql = "UPDATE element SET standard = ? WHERE molecule_id = ? AND standard != 0"
			self.cursor.execute(sql, [0, molecule_id])
		sql = "INSERT INTO element (id, molecule_id, type, accession, symbol, description, start, end, strand, sequence, standard, parent_id) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);"
		self.cursor.execute(sql, [None, molecule_id, type, accession, symbol, description, start, end, strand, sequence, standard, parent_id])
		sql = "SELECT id FROM element WHERE accession = ?"
		mid = self.cursor.execute(sql, [accession]).fetchone()['id']
		if parts is not None:
			for part in parts:
				sql = "INSERT OR IGNORE INTO elempart (element_id, start, end, strand, base, segment) VALUES(?, ?, ?, ?, ?, ?);"
				self.cursor.execute(sql, [mid] + part)
		return mid

	def insert_variant(self, alignment_id, ref, alt, start, end, parent_id=""):
		vid = self.get_variant_id(ref, alt, start, end, parent_id)
		if not vid:
			sql = "INSERT OR IGNORE INTO variant (id, start, end, ref, alt, parent_id) VALUES(?, ?, ?, ?, ?, ?);"
			self.cursor.execute(sql, [None, start, end, ref, alt, parent_id])
			vid = self.cursor.lastrowid
		sql = "INSERT OR IGNORE INTO alignment2variant (alignment_id, variant_id) VALUES(?, ?);"
		self.cursor.execute(sql, [alignment_id, vid])
		return vid

	# DELETING DATA

	def delete_sample(self, sample_name):
		sql = "DELETE FROM sample WHERE name = ?;"
		self.cursor.execute(sql, [sample_name])

	def delete_property(self, property_name):
		sql = "DELETE FROM sample2property WHERE property_id = ?;"
		self.cursor.execute(sql, [self.properties[property_name]['id']])
		sql = "DELETE FROM property WHERE name = ?;"
		self.cursor.execute(sql, [property_name])

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
		sql = "SELECT id FROM alignment WHERE seqhash = ? AND element_id = ? LIMIT 1;"
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

	def get_annotation(self, reference_accession=None, molecule_accession=None, element_accession=None, fields=['*']):
		conditions = []
		vals = []
		if reference_accession:
			 conditions.append("\"reference.accession\" = ?")
			 vals.append(reference_accession)
		else:
			 conditions.append("\"reference.standard\" = ?")
			 vals.append(1)
		if molecule_accession:
			 conditions.append("\"molecule.accession\" = ?")
			 vals.append(molecule_accession)
		else:
			 conditions.append("\"molecule.standard\" = ?")
			 vals.append(1)
		if element_accession:
			 conditions.append("\"element.accession\" = ?")
			 vals.append(element_accession)
		else:
			 conditions.append("\"element.type\" = ?")
			 vals.append('source')
		sql = "SELECT " + ", ".join(fields) + " FROM referenceView WHERE " + " AND ".join(conditions)
		return self.cursor.execute(sql, vals).fetchone()

	def get_alignment_data(self, sample_id, element_id, *fields, limit=1):
		if not fields:
			fields = ['*']
		sql = "SELECT " + ", ".join(fields) + " FROM alignmentView WHERE \"sample.id\" = ? AND \"element.id\" = ? LIMIT "+ LIMIT +";"
		return self.cursor.execute(sql, [sample_id, element_id]).fetchall()

	def get_alignment_id(self, seqhash, element_id):
		sql = "SELECT id FROM alignment WHERE \"seqhash\" = ? AND \"element.id\" = ? LIMIT 1;"
		row = self.cursor.execute(sql, [seqhash, element_id]).fetchone()
		if row:
			return row['id']
		return None

	def get_variant_id(self, ref, alt, start, end, parent_id):
		sql = "SELECT id FROM variant WHERE start = ? AND end = ? AND ref = ? AND alt = ? AND parent_id = ?;"
		row = self.cursor.execute(sql, [start, end, ref, alt, parent_id]).fetchone()
		return None if row is None else row['id']

	def iter_dna_variants(self, sample_name, *element_ids):
		sql = "SELECT \"element.id\", \"element.symbol\", \"variant.start\", \"variant.end\", \"variant.ref\", \"variant.alt\" FROM variantView WHERE \"sample.name\" = ? AND \"element.id\" IN (" + ", ".join(['?'] * len(element_ids)) + ");"
		for row in self.cursor.execute(sql, [sample_name] + [str(x) for x in element_ids]):
			if row["variant.start"] is not None:
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

	def count_property(self, property_name, distinct=False, ignore_standard=False):
		d = "DISTINCT " if distinct else ""
		c = "WHERE property_id = ?"
		v = [self.properties[property_name]['id']]
		if ignore_standard and self.properties[property_name]['standard'] is not None:
			c += " AND value_" + self.properties[property_name]['datatype'] + " != ?"
			v.append(self.properties[property_name]['standard'])
		sql = "SELECT COUNT(" + d + " value_" + self.properties[property_name]['datatype'] + ") as count FROM sample2property " + c + ";"
		return self.cursor.execute(sql, v).fetchone()['count']

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

	def get_seqhash(self, sample_name):
		sql = "SELECT \"sample.seqhash\" FROM sequenceView WHERE sample.name = ?"
		return [ x[sample.hash] for x in self.cursor.execute(sql, [sample_name]).fetchall() if x is not None ]

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
		sql = "SELECT sequence, type, id, parent_id FROM element WHERE id = ?"
		row = self.cursor.execute(sql, [element_id]).fetchone()
		if not row:
			return None
		element_id = row['id']
		while row and row['type'] not in {"source", "CDS"}:
			sql = "SELECT sequence, type, parent_id FROM element WHERE id = ?"
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

		intersect_sqls = []
		intersect_vals = []
		except_sqls = []

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
		errmsg =  "query error: date or date range expected for field " + field + " (got: "
		data = defaultdict(set)
		conditions = []
		vallist = []

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
		targetfield = "value_" + self.properties[name]['datatype']

		# query dates
		if self.properties[name]['querytype'] == "date":
			a, b = self.query_dates(targetfield, *vals)
			conditions.append(a)
			valueList.extend(b)

		# query numeric
		elif self.properties[name]['querytype'] == "numeric":
			a,b = self.query_numeric(targetfield, vals)
			conditions.append(a)
			valueList.extend(b)

		# query text
		elif self.properties[name]['querytype'] == "string":
			a,b = self.query_string(targetfield, vals)
			conditions.append(a)
			valueList.extend(b)

		# query zip
		elif self.properties[name]['querytype'] == "zip":
			a,b = self.query_zip(targetfield, vals)
			conditions.append(a)
			valueList.extend(b)

		else:
			sys.exit("error: unknown query type.")

		return "SELECT \"sample.id\" FROM propertyView WHERE " + " AND ".join(conditions), valueList

	def query_profile(self, *vars, reference_accession=None, molecule_accession=None, element_accession=None):
		iupac_nt_code = { "A": set("A"), "C": set("C"), "G": set("G"), "T": set("T"), "R": set("AGR"), "Y": set("CTY"), "S": set("GCS"), "W": set("ATW"), "K": set("GTK"), "M": set("ACM"), "B": set("CGTB"), "D": set("AGTD"), "H": set("ACTH"), "V": set("ACGV"), "N": set("ACGTRYSWKMBDHVN") }
		iupac_aa_code = { "A": set("A"), "R": set("R"), "N": set("N"), "D": set("D"), "C": set("C"), "Q": set("Q"), "E": set("E"), "G": set("G"), "H": set("H"), "I": set("I"), "L": set("L"), "K": set("K"), "M": set("M"), "F": set("F"), "P": set("P"), "S": set("S"), "T": set("T"), "W": set("W"), "Y": set("Y"), "V": set("V"), "U": set("U"), "O": set("O"), "B": set("DNB"), "Z": set("EQZ"), "J": set("ILJ"), "Φ": set("VILFWYMΦ"), "Ω": set("FWYHΩ"), "Ψ": set("VILMΨ"), "π": set("PGASπ"), "ζ": set("STHNQEDKRζ"), "+": set("KRH+"), "-": set("DE-"), "X": set("ARNDCQEGHILKMFPSTWYVUOBZJΦΩΨπζ+-X") }
		del_regex = re.compile(r'^(|[^:]+:)?([^:]+:)?del:(=?[0-9]+)(|-=?[0-9]+)?$')
		snv_regex = re.compile(r'^(|[^:]+:)?([^:]+:)?([A-Z]+)([0-9]+)(=?[A-Z]+)$')

		# set variants and generate sql
		base_sql = "SELECT \"sample.id\" FROM variantView WHERE "
		intersect_sqls = []
		intersect_vals = []
		except_sqls = []
		except_vals = []
		for var in vars:
			c = []
			v = []

			if var.startswith("^"):
				var = var[1:]
				negate = True
			else:
				negate = False

			## variant typing
			if match := snv_regex.match(var):
				snv = True
			elif match := del_regex.match(var):
				snv = False
			else:
				sys.exit("input error: " + var + " is not a valid variant definition.")

			## set molecule
			if match.group(1):
				c.append("\"molecule.symbol\" = ?")
				v.append(match.group(1)[:-1])
			else:
				c.append("\"molecule.standard\" = ?")
				v.append(1)

			## set element
			if match.group(2):
				c.append("\"element.type\" = ?")
				v.append("protein")
				c.append("\"element.symbol\" = ?")
				v.append(match.group(2)[:-1])
				code = iupac_aa_code
			else:
				c.append("\"element.standard\" = ?")
				v.append(1)
				code = iupac_nt_code

			## snp and insert handling
			if snv:
				c.append("\"variant.start\" = ?")
				v.append(int(match.group(4))-1)
				c.append("\"variant.end\" = ?")
				v.append(int(match.group(4)))
				c.append("\"variant.ref\" = ?")
				v.append(match.group(3))

				### explicit alternate allele
				if match.group(5).startswith("="):
					c.append("\"variant.alt\" = ?")
					v.append(match.group(5)[1:])

				### potentially ambiguous alternate snp
				elif len(match.group(5)) == 1:
					l = len(code[match.group(5)])
					if l == 1:
						c.append("\"variant.alt\" = ?")
						v.append(match.group(5))
					else:
						c.append("\"variant.alt\" IN (" + ", ".join(['?'] * l) + ")")
						v.extend(code[match.group(5)])

				### potentially ambiguous alternate insert
				else:
					a = ["".join(x) for x in itertools.product(*[code[x] for x in match.group(5)])]
					l = len(a)
					if l == 1:
						c.append("\"variant.alt\" = ?")
						v.extend(a)
					else:
						c.append("\"variant.alt\" IN (" + ", ".join(['?'] * l) + ")")
						v.extend(a)

			## deletion handling
			else:
				s = match.group(3)
				e = match.group(4)[1:]

				if s.startswith("="):
					s = s[1:]
					c.append("\"variant.start\" = ?")
				else:
					c.append("\"variant.start\" <= ?")

				if e.startswith("="):
					e = e[1:]
					c.append("\"variant.end\" = ?")
				else:
					c.append("\"variant.end\" >= ?")
				v.append(int(s) - 1)
				v.append(int(e))

				c.append("\"variant.alt\" = ?")
				v.append(" ")

			## assemble sub-sql
			if negate:
				except_sqls.append(base_sql + " AND ".join(c))
				except_vals.extend(v)
			else:
				intersect_sqls.append(base_sql + " AND ".join(c))
				intersect_vals.extend(v)

		# assemble final sql
		if not intersect_sqls:
			intersect_sqls = [base_sql + "1"]

		sql = " INTERSECT ".join(intersect_sqls)

		if except_sqls:
			sql += " EXCEPT " + " EXCEPT ".join(except_sqls)

		return sql, intersect_vals + except_vals

	def match(self, *profiles, properties=None, reference_accession=None, molecule_accession=None, element_accession = None):
		#collecting sqls for metadata-based filtering
		property_sqls = []
		property_vals = []
		if properties:
			for pname, vals in properties.items():
				sql, val = self.query_metadata(pname, *vals)
				property_sqls.append(sql)
				property_vals.extend(val)

		property_sqls = " INTERSECT ".join(property_sqls)
		#collecting sqls for genomic profile based filtering
		profile_sqls = []
		profile_vals = []
		for profile in profiles:
			sql, val = self.query_profile(*profile, reference_accession=reference_accession)
			profile_sqls.append(sql)
			profile_vals.extend(val)

		if len(profiles) == 1:
			profile_sqls = profile_sqls[0]
		elif len(profiles) > 1:
			profile_sqls = " UNION ".join(["SELECT * FROM (" + x + ")" for x in profile_sqls])
		else:
			profile_sqls = ""

		# assembling final compound query
		if property_sqls and profile_sqls:
			if len(profiles) > 1:
				sql = property_sqls + " INTERSECT SELECT * FROM (" + profile_sqls + ")"
			else:
				sql = property_sqls + " INTERSECT " + profile_sqls
		else:
			sql = property_sqls + profile_sqls

		print(sql)
		return self.cursor.execute(sql, property_vals + profile_vals).fetchall()

	# UPDATE DATA


	# MISCg

	@staticmethod
	def optimize(dbfile):
		with sqlite3.connect(dbfile) as con:
			con.executescript("PRAGMA analysis_limit=400;")
			con.executescript("PRAGMA optimize;")
			con.executescript("VACUUM;")

	@staticmethod
	def dict_factory(cursor, row):
		d = OrderedDict()
		for idx, col in enumerate(cursor.description):
			d[col[0]] = row[idx]
		return d
