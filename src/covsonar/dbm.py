from collections import defaultdict
import itertools
import logging
import os
import pkgutil
import re
import sqlite3
import sys
from urllib.parse import quote as urlquote

from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation
from Bio.SeqFeature import FeatureLocation
import pandas as pd

# COMPATIBILITY
SUPPORTED_DB_VERSION = 4


class sonarDBManager:
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

    def __init__(
        self, dbfile, timeout=-1, readonly=True, debug=False, autocreate=False
    ):
        logging.basicConfig(format="%(asctime)s %(message)s")
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
        self.__illegal_properties = {}
        self.__default_reference_id = False
        self.__lineage_sublineage_dict = None
        self.__rootdir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        self.lineagewithsublineages = os.path.join(
            self.__rootdir, "lib", "lineage.all.tsv"
        )

        self.operators = {
            "genuine": {
                "=": "=",
                ">": ">",
                "<": "<",
                ">=": ">=",
                "<=": "<=",
                "IN": "IN",
                "LIKE": "LIKE",
                "BETWEEN": "BETWEEN",
            },
            "negated": {
                "=": "!=",
                ">": "<=",
                "<": ">=",
                ">=": "<",
                "<=": ">",
                "IN": "NOT IN",
                "LIKE": "NOT LIKE",
                "BETWEEN": "NOT BETWEEN",
            },
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
        con = sqlite3.connect(
            self.__uri + "?mode=" + self.__mode,
            self.__timeout,
            isolation_level=None,
            uri=True,
        )
        if self.debug:
            con.set_trace_callback(logging.warning)
        con.row_factory = self.dict_factory
        cur = con.cursor()
        return con, cur

    def start_transaction(self):
        self.cursor.execute("BEGIN DEFERRED")

    def new_transaction(self):
        self.commit()
        self.start_transaction()

    def commit(self):
        self.connection.commit()

    def rollback(self):
        self.connection.rollback()

    def close(self):
        self.connection.close()

    def get_db_version(self):
        return self.cursor.execute("pragma user_version").fetchone()["user_version"]

    def check_db_compatibility(self):
        dbver = self.get_db_version()
        if dbver != SUPPORTED_DB_VERSION:
            sys.exit(
                "compatibility error: the given database is not compatible with this version of sonar (database version: "
                + str(dbver)
                + "; supported database version: "
                + str(SUPPORTED_DB_VERSION)
                + ")"
            )

    @property
    def lineage_sublineage_dict(self):
        if not self.__lineage_sublineage_dict:
            df = pd.read_sql(
                "SELECT * FROM lineages", self.connection
            )  # pd.read_csv(self.lineagewithsublineages, sep='\t')
            self.__lineage_sublineage_dict = dict(zip(df.lineage, df.sublineage))
        return self.__lineage_sublineage_dict

    @property
    def default_reference_id(self):
        if self.__default_reference is False:
            sql = "SELECT accession FROM reference WHERE standard = 1 LIMIT 1"
            self.__default_reference = self.cursor.execute(sql).fetchone()
        return self.__default_reference

    @property
    def properties(self):
        if self.__properties is False:
            sql = "SELECT * FROM property"
            rows = self.cursor.execute(sql).fetchall()
            self.__properties = {} if not rows else {x["name"]: x for x in rows}
        return self.__properties

    # SETUP

    @staticmethod
    def get_uri(dbfile):
        return "file:" + urlquote(dbfile)

    @staticmethod
    def setup(filename, debug=False):
        sql = pkgutil.get_data(__name__, "data/db.sqlite").decode()
        uri = sonarDBManager.get_uri(filename)
        with sqlite3.connect(uri + "?mode=rwc", uri=True) as con:
            if debug:
                con.set_trace_callback(logging.warning)
            con.executescript(sql)

    def add_codon(self, translation_table, codon, aa):
        sql = "INSERT OR IGNORE INTO translation (id, codon, aa) VALUES(?, ?, ?);"
        self.cursor.execute(sql, [translation_table, codon, aa])

    def add_property(self, name, datatype, querytype, description, standard=None):
        name = name.upper()
        if name in self.__illegal_properties:
            sys.exit(
                "error: '"
                + str(name)
                + "' is reserved and cannot be used as property name"
            )
        if not re.match("^[a-zA-Z0-9_]+$", name):
            sys.exit(
                "error: invalid property name (property names can contain only letters, numbers and underscores)"
            )
        if name in self.properties:
            sys.exit(
                "error: a property named "
                + name
                + " already exists in the given database."
            )
        try:
            sql = "INSERT INTO property (name, datatype, querytype, description, standard) VALUES(?, ?, ?, ?, ?);"
            self.cursor.execute(sql, [name, datatype, querytype, description, standard])
            self.__properties = False
            pid = self.properties[name]["id"]
            if standard is not None:
                sql = (
                    "INSERT INTO sample2property (property_id, value_"
                    + self.properties[name]["datatype"]
                    + ", sample_id) SELECT ?, ?, id FROM sample WHERE 1"
                )
                vals = [pid, standard]
                self.cursor.execute(sql, vals)
        except sqlite3.Error as error:
            print("Failed to insert data into sqlite table", error)
        return pid

    def add_translation_table(self, translation_table):
        sql = "SELECT COUNT(*) FROM translation  WHERE id = ?;"
        if self.cursor.execute(sql, [translation_table]).fetchone()["COUNT(*)"] != 4096:
            for codon in itertools.product("ATGCRYSWKMBDHVN-", repeat=3):
                codon = "".join(codon)
                try:
                    aa = str(Seq.translate(codon, table=translation_table))
                except Exception:
                    aa = ""
                self.add_codon(translation_table, codon, aa)

    def add_reference(
        self, accession, description, organism, translation_table, standard=0
    ):
        self.add_translation_table(translation_table)
        if standard:
            sql = "UPDATE reference SET standard = 0 WHERE standard != 0"
            self.cursor.execute(sql)
        sql = "INSERT INTO reference (id, accession, description, organism, translation_id, standard) VALUES(?, ?, ?, ?, ?, ?);"
        self.cursor.execute(
            sql, [None, accession, description, organism, translation_table, standard]
        )
        return self.cursor.lastrowid

    # INSERTING DATA

    def insert_property(self, sample_id, property_name, property_value):
        sql = (
            "INSERT OR REPLACE INTO sample2property (sample_id, property_id, value_"
            + self.properties[property_name]["datatype"]
            + ") VALUES(?, ?, ?);"
        )
        # print([sample_id, self.properties[property_name]['id'], property_value])
        self.cursor.execute(
            sql, [sample_id, self.properties[property_name]["id"], property_value]
        )
        return

    def insert_sequence(self, seqhash):
        sql = "INSERT OR IGNORE INTO sequence (seqhash) VALUES(?);"
        self.cursor.execute(sql, [seqhash])
        return seqhash

    def insert_sample(self, sample_name, seqhash):
        self.insert_sequence(seqhash)
        sql = "INSERT OR REPLACE INTO sample (name, seqhash, datahash) VALUES(?, ?, ?);"
        self.cursor.execute(sql, [sample_name, seqhash, ""])
        sql = "SELECT id FROM sample WHERE name = ?;"
        sid = self.cursor.execute(sql, [sample_name]).fetchone()["id"]
        for pname in self.properties:
            if not self.properties[pname]["standard"] is None:
                self.insert_property(sid, pname, self.properties[pname]["standard"])
        return

    def insert_alignment(self, seqhash, element_id):
        sql = (
            "INSERT OR IGNORE INTO alignment (id, seqhash, element_id) VALUES(?, ?, ?);"
        )
        self.cursor.execute(sql, [None, seqhash, element_id])
        sql = "SELECT id FROM alignment WHERE element_id = ? AND seqhash = ?;"
        return self.cursor.execute(sql, [element_id, seqhash]).fetchone()["id"]

    def insert_molecule(
        self,
        reference_id,
        type,
        accession,
        symbol,
        description,
        segment,
        length,
        standard=0,
    ):
        if symbol.strip() == "":
            symbol = accession
        if standard:
            sql = "UPDATE molecule SET standard = ? WHERE reference_id = ? AND standard = 1"
            self.cursor.execute(sql, [0, reference_id])
        sql = "INSERT INTO molecule (id, reference_id, type, accession, symbol, description, segment, length, standard) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?);"
        self.cursor.execute(
            sql,
            [
                None,
                reference_id,
                type,
                accession,
                symbol,
                description,
                segment,
                length,
                standard,
            ],
        )
        sql = "SELECT id FROM molecule WHERE accession = ?"
        mid = self.cursor.execute(sql, [accession]).fetchone()["id"]
        return mid

    def insert_element(
        self,
        molecule_id,
        type,
        accession,
        symbol,
        description,
        start,
        end,
        strand,
        sequence,
        standard=0,
        parent_id=" ",
        parts=None,
    ):
        if symbol.strip() == "":
            symbol = accession
        if standard:
            sql = (
                "UPDATE element SET standard = ? WHERE molecule_id = ? AND standard = 1"
            )
            self.cursor.execute(sql, [0, molecule_id])
        sql = "INSERT INTO element (id, molecule_id, type, accession, symbol, description, start, end, strand, sequence, standard, parent_id) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);"
        self.cursor.execute(
            sql,
            [
                None,
                molecule_id,
                type,
                accession,
                symbol,
                description,
                start,
                end,
                strand,
                sequence,
                standard,
                parent_id,
            ],
        )
        sql = "SELECT id FROM element WHERE accession = ?"
        mid = self.cursor.execute(sql, [accession]).fetchone()["id"]
        if parts is not None:
            for part in parts:
                sql = "INSERT OR IGNORE INTO elempart (element_id, start, end, strand, base, segment) VALUES(?, ?, ?, ?, ?, ?);"
                self.cursor.execute(sql, [mid] + part)
        return mid

    def insert_variant(
        self, alignment_id, element_id, ref, alt, start, end, label, parent_id=""
    ):
        sql = "INSERT OR IGNORE INTO variant (id, element_id, start, end, ref, alt, label, parent_id) VALUES(?, ?, ?, ?, ?, ?, ?, ?);"
        self.cursor.execute(
            sql, [None, element_id, start, end, ref, alt, label, parent_id]
        )
        vid = self.get_variant_id(element_id, start, end, ref, alt)
        sql = "INSERT OR IGNORE INTO alignment2variant (alignment_id, variant_id) VALUES(?, ?);"
        self.cursor.execute(sql, [alignment_id, vid])
        return vid

    # DELETING DATA

    def delete_samples(self, *sample_names):
        sql = (
            "DELETE FROM sample WHERE name IN ("
            + ", ".join(["?"] * len(sample_names))
            + ");"
        )
        self.cursor.execute(sql, sample_names)

    def delete_property(self, property_name):
        sql = "DELETE FROM sample2property WHERE property_id = ?;"
        self.cursor.execute(sql, [self.properties[property_name]["id"]])
        sql = "DELETE FROM property WHERE name = ?;"
        self.cursor.execute(sql, [property_name])

    # SELECTING DATA

    def sample_exists(self, sample_name):
        sql = "SELECT EXISTS(SELECT 1 FROM sample WHERE name=? LIMIT 1) as found"
        return bool(self.cursor.execute(sql, [sample_name]).fetchone()["found"])

    def get_sample_id(self, sample_name):
        sql = "SELECT id FROM sample WHERE name = ? LIMIT 1;"
        row = self.cursor.execute(sql, [sample_name]).fetchone()
        return row["id"] if row else None

    def get_sample_data(self, sample_name):
        sql = "SELECT id, seqhash FROM sample WHERE name = ? LIMIT 1;"
        row = self.cursor.execute(sql, [sample_name]).fetchone()
        return (row["id"], row["seqhash"]) if row else (None, None)

    def seq_exists(self, seqhash):
        sql = "SELECT id FROM sample WHERE seqhash=? LIMIT 1;"
        return False if self.cursor.execute(sql, [seqhash]).fetchone() is None else True

    def get_alignment_id(self, seqhash, element_id):
        sql = "SELECT id FROM alignment WHERE element_id = ? AND seqhash = ? LIMIT 1;"
        row = self.cursor.execute(sql, [element_id, seqhash]).fetchone()
        return None if row is None else row["id"]

    def get_default_reference_accession(self):
        sql = "SELECT accession FROM reference WHERE standard=1"
        return self.cursor.execute(sql).fetchone()["accession"]

    def get_molecule_ids(self, reference_accession=None):
        if reference_accession:
            condition = '"reference.accession" = ?'
            val = [reference_accession]
        else:
            condition = '"reference.standard" = ?'
            val = [1]
        sql = (
            'SELECT "molecule.accession", "molecule.id" FROM referenceView WHERE '
            + condition
        )
        return {
            x["molecule.accession"]: x["molecule.id"]
            for x in self.cursor.execute(sql, val).fetchall()
            if x is not None
        }

    def get_molecule_id(self, molecule_accession):
        sql = "SELECT id FROM molecule WHERE accession = ?;"
        row = self.cursor.execute(sql, [molecule_accession]).fetchone()
        if row:
            row = row["id"]
        return row

    def get_molecule_data(self, *fields, reference_accession=None):
        if not fields:
            fields = "*"
        elif '"molecule.accession"' not in fields:
            fields = list(fields) + ['"molecule.accession"']
        if reference_accession:
            condition = '"reference.accession" = ?'
            vals = [reference_accession]
        else:
            condition = '"reference.standard" = ?'
            vals = [1]
        sql = (
            "SELECT "
            + ", ".join(fields)
            + " FROM referenceView WHERE "
            + condition
            + ";"
        )
        row = self.cursor.execute(sql, vals).fetchall()
        if row:
            return {x["molecule.accession"]: x for x in row}
        return {}

    def get_elements(self, molecule_id, *types):
        sql = "SELECT * FROM element WHERE molecule_id = ?"
        if types:
            sql += " AND type IN (" + ", ".join(["?"] * len(types)) + ");"
        row = self.cursor.execute(sql, [molecule_id] + list(types)).fetchall()
        if not row:
            return []
        return row

    def get_element_ids(self, reference_accession=None, type=None):
        molecule_ids = list(
            self.get_molecule_ids(reference_accession=reference_accession).values()
        )
        sql = (
            'SELECT id FROM element WHERE "molecule_id" IN ('
            + ", ".join(["?"] * len(molecule_ids))
            + ")"
        )
        if type:
            sql += " AND type = ?"
            molecule_ids.append(type)
        row = self.cursor.execute(sql, molecule_ids).fetchall()
        if not row:
            return []
        return [x["id"] for x in row]

    def get_source(self, molecule_id):
        return self.get_elements(molecule_id, "source")[0]

    def get_annotation(
        self,
        reference_accession=None,
        molecule_accession=None,
        element_accession=None,
        element_type=None,
        fields=["*"],
    ):
        conditions = []
        vals = []
        if reference_accession:
            conditions.append('"reference.accession" = ?')
            vals.append(reference_accession)
        else:
            conditions.append('"reference.standard" = ?')
            vals.append(1)
        if molecule_accession:
            conditions.append('"molecule.accession" = ?')
            vals.append(molecule_accession)
        else:
            conditions.append('"molecule.standard" = ?')
            vals.append(1)
        if element_accession:
            conditions.append('"element.accession" = ?')
            vals.append(element_accession)
        elif not element_type:
            conditions.append('"element.type" = ?')
            vals.append("source")
        if element_type:
            conditions.append('"element.type" = ?')
            vals.append(element_type)
        sql = (
            "SELECT "
            + ", ".join(fields)
            + " FROM referenceView WHERE "
            + " AND ".join(conditions)
            + ' ORDER BY "reference.id" ASC, "molecule.id" ASC, "element.id" ASC, "element.segment" ASC'
        )
        return self.cursor.execute(sql, vals).fetchall()

    def get_alignment_data(self, sample_name, reference_accession=None):
        if reference_accession is None:
            sql = 'SELECT "reference.accession" as acc FROM alignmentView WHERE "sample.name" = ? LIMIT 1'
            reference_accession = self.cursor.execute(sql, [sample_name]).fetchone()
            if reference_accession is None:
                return ""
            reference_accession = reference_accession["acc"]
        sql = 'SELECT "element.sequence", "element.symbol", "element.id" FROM alignmentView WHERE "sample.name" = ? AND "reference.accession" = ? '
        return self.cursor.execute(sql, [sample_name, reference_accession])

    def get_variant_id(self, element_id, start, end, ref, alt):
        sql = "SELECT id FROM variant WHERE element_id = ? AND start = ? AND end = ? AND ref = ? AND alt = ?;"
        row = self.cursor.execute(sql, [element_id, start, end, ref, alt]).fetchone()
        return None if row is None else row["id"]

    def iter_dna_variants(self, sample_name, *element_ids):
        if len(element_ids) == 1:
            condition = " = ?"
        elif len(element_ids) > 1:
            condition = " IN (" + ", ".join(["?"] * len(element_ids)) + ")"
        else:
            condition = ""
        # sql = "SELECT \"element.id\", \"variant.start\", \"variant.end\", \"variant.ref\", \"variant.alt\" FROM variantView WHERE \"sample.name\" = ? AND \"element.id\"" + condition
        sql = (
            """
            SELECT  variant.element_id as \"element.id\",
                    variant.start as \"variant.start\",
                    variant.end as  \"variant.end\",
                    variant.ref as  \"variant.ref\",
                    variant.alt as \"variant.alt\"
                    FROM
                        ( SELECT sample.seqhash
                        FROM sample
                        WHERE sample.name = ?
                        ) AS sample_T
                    INNER JOIN alignment
                        ON sample_T.seqhash == alignment.seqhash
                    INNER JOIN alignment2variant
                        ON alignment.id == alignment2variant.alignment_id
                    INNER JOIN    variant
                        ON alignment2variant.variant_id == variant.id
                        WHERE  variant.element_id """
            + condition
        )

        for row in self.cursor.execute(sql, [sample_name] + list(element_ids)):
            if row["variant.start"] is not None:
                yield row

    def iter_profile(self, sample_name, element_id):
        sql = 'SELECT * FROM variantView WHERE "sample.name" = ? AND "element.id" = ? ORDER BY "element.id", "variant.start";'
        for row in self.cursor.execute(sql, [sample_name, element_id]):
            yield row

    def iter_protein_variants(self, sample_name, *element_ids):
        sql = (
            "SELECT element_id, element.symbol, start, end, ref, alt FROM variant WHERE name = ? AND element.id IN ("
            + ", ".join(["?"] * len(element_ids))
            + ") AND type = 'cds' LEFT JOIN element ON element.id = element_id;"
        )
        for row in self.cursor.execute(sql, [sample_name] + element_ids):
            yield row

    def get_frameshift_vars(self, sample_name, *element_ids):
        sql = (
            "SELECT element_id, element.symbol, start, end, ref, alt FROM variant WHERE name = ? AND element.id IN ("
            + ", ".join(["?"] * len(element_ids))
            + ") AND frameshift = 1 LEFT JOIN element ON element.id = element_id;"
        )
        for row in self.cursor.execute(sql, [sample_name] + element_ids):
            yield row

    def count_samples(self):
        sql = "SELECT COUNT(*) FROM sample;"
        return self.cursor.execute(sql).fetchone()["COUNT(*)"]

    def count_sequences(self):
        sql = "SELECT COUNT(DISTINCT seqhash) FROM sample;"
        return self.cursor.execute(sql).fetchone()["COUNT(DISTINCT seqhash)"]

    def count_property(self, property_name, distinct=False, ignore_standard=False):
        d = "DISTINCT " if distinct else ""
        c = "WHERE property_id = ?"
        v = [self.properties[property_name]["id"]]
        if ignore_standard and self.properties[property_name]["standard"] is not None:
            c += " AND value_" + self.properties[property_name]["datatype"] + " != ?"
            v.append(self.properties[property_name]["standard"])
        sql = (
            "SELECT COUNT("
            + d
            + " value_"
            + self.properties[property_name]["datatype"]
            + ") as count FROM sample2property "
            + c
            + ";"
        )
        return self.cursor.execute(sql, v).fetchone()["count"]

    def get_properties(self, sample_name, property_name=None):
        if property_name:
            sql = 'SELECT * FROM propertyView WHERE "sample.name" = ? AND property_name = ?;'
            vals = [sample_name, property_name]
        else:
            sql = 'SELECT * FROM propertyView WHERE "sample.name" = ?'
            vals = [sample_name]
        properties = {}
        for row in self.cursor.execute(sql, vals):
            if row is None:
                return {}
            properties[row["property.name"]] = row[
                "sample2property.value_" + row["property.datatype"]
            ]
        return properties

    def count_properties(self):
        sql = "SELECT COUNT(*) as count FROM propertyView GROUP BY property_name;"
        return self.cursor.execute(sql).fetchall()

    def get_seqhash(self, sample_name):
        sql = 'SELECT "sample.seqhash" FROM sequenceView WHERE sample.name = ?'
        return [
            x["sample.hash"]
            for x in self.cursor.execute(sql, [sample_name]).fetchall()
            if x is not None
        ]

    def get_translation_dict(self, translation_id):
        sql = "SELECT codon, aa FROM translation WHERE id = ?"
        return {
            x["codon"]: x["aa"]
            for x in self.cursor.execute(sql, [translation_id]).fetchall()
        }

    def get_earliest_import(self):
        sql = "SELECT MIN(imported) as import FROM genome WHERE import IS NOT NULL;"
        return self.cursor.execute(sql).fetchone()["import"]

    def get_latest_import(self):
        sql = "SELECT MAX(imported) as import FROM genome WHERE import IS NOT NULL;"
        return self.cursor.execute(sql).fetchone()["import"]

    def get_earliest_date(self):
        sql = "SELECT MIN(date) as date FROM genome WHERE date IS NOT NULL;"
        return self.cursor.execute(sql).fetchone()["date"]

    def get_latest_date(self):
        sql = "SELECT MAX(date) as date FROM genome WHERE date IS NOT NULL;"
        return self.cursor.execute(sql).fetchone()["date"]

    def get_element_parts(self, element_id=None):
        sql = "SELECT start, end, strand FROM elempart WHERE element_id = ? ORDER BY segment"
        return self.cursor.execute(sql, [element_id]).fetchall()

    def get_sequence(self, element_id=None):
        sql = "SELECT sequence, type FROM element WHERE id = ?"
        row = self.cursor.execute(sql, [element_id]).fetchone()
        return None if row is None else row["sequence"]

    def extract_sequence(self, element_id=None, translation_table=None):
        sql = "SELECT sequence, type, id, parent_id FROM element WHERE id = ?"
        row = self.cursor.execute(sql, [element_id]).fetchone()
        if not row:
            return None
        element_id = row["id"]
        while row and row["type"] not in {"source", "CDS"}:
            sql = "SELECT sequence, type, parent_id FROM element WHERE id = ?"
            row = self.cursor.execute(sql, [row["parent_id"]]).fetchone()
        sequence = row["sequence"]
        parts = []
        for part in self.get_element_parts(element_id):
            parts.append(
                FeatureLocation(part["start"], part["end"], strand=part["strand"])
            )
        feat = CompoundLocation(parts) if len(parts) > 1 else parts[0]
        if translation_table is None:
            return str(feat.extract(sequence))
        return str(
            Seq(feat.extract(sequence)).translate(
                table=translation_table, stop_symbol=""
            )
        )

    # Add/Update lineages into Table
    def add_update_lineage(self, _df):
        print("Prepare:", len(_df))

        # _df.to_sql(name='lineages', con = self.connection, if_exists='replace',
        # 			 index=False)
        sql = "INSERT OR REPLACE INTO lineages (lineage, sublineage) VALUES(?, ?);"
        for ind in _df.index:
            self.cursor.execute(sql, (_df["lineage"][ind], _df["sublineage"][ind]))
        # self.commit()
        # sql = "INSERT OR REPLACE INTO sample2property (sample_id, property_id, value_" + self.properties[property_name]['datatype'] + ") VALUES(?, ?, ?);"
        # print([sample_id, self.properties[property_name]['id'], property_value])
        # self.cursor.execute(sql, [sample_id, self.properties[property_name]['id'], property_value])

    # MATCHING PROFILES
    def get_operator(self, val):
        if val.startswith("^"):
            return val[1:], self.operators["negated"]
        return val, self.operators["genuine"]

    def get_conditional_expr(self, field, operator, *vals):
        condlist = []
        vallist = []

        # transforming operators
        if operator == "=" and len(vals) > 1:
            operator = "IN"
        elif operator == "!=" and len(vals) > 1:
            operator = "NOT IN"

        # creating conditions
        if operator == "IN" or operator == "NOT IN":
            condlist.append(
                field + " " + operator + "(" + ", ".join(["?"] * len(vals)) + ")"
            )
            vallist.extend(vals)

        elif operator == "BETWEEN" or operator == "NOT BETWEEN":
            for val in vals:
                condlist.append(field + " " + operator + " ? AND ?")
                vallist.extend(val)

        else:
            condlist.extend([field + " " + operator + " ?"] * len(vals))
            vallist.extend(vals)

        return condlist, vallist

    def query_numeric(self, field, *vals, link="OR"):
        link = " " + link.strip() + " "
        defaultop = "="
        op1 = re.compile(r"^(\^*)((?:>|>=|<|<=|!=|=)?)(-?[1-9]+[0-9]*)$")
        op2 = re.compile(r"^(\^*)(-?[1-9]+[0-9]*):(-?[1-9]+[0-9]*)$")
        errmsg = (
            "query error: numeric value or range expected for field " + field + "(got: "
        )
        data = defaultdict(set)
        conditions = []
        vallist = []

        for val in vals:
            val = str(val).strip()

            # single value
            if ":" not in val:
                match = op1.match(val)
                if not match:
                    sys.exit(errmsg + val + ")")
                op = match.group(2) if match.group(2) else defaultop
                operator = (
                    self.operators["negated"][op]
                    if match.group(1)
                    else self.operators["genuine"][op]
                )
                val = match.group(3)
                data[operator].add(val)

            # range
            else:
                match = op2.match(val)
                if not match:
                    sys.exit(errmsg + val + ")")
                operator = (
                    self.operators["negated"]["BETWEEN"]
                    if match.group(1)
                    else self.operators["genuine"]["BETWEEN"]
                )
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
        op1 = re.compile(r"^(\^*)((?:>|>=|<|<=|!=|=)?)(-?[1-9]+[0-9]*(?:.[0-9]+)*)$")
        op2 = re.compile(
            r"^(\^*)(-?[1-9]+[0-9]*(?:.[0-9]+)*):(-?[1-9]+[0-9]*(?:.[0-9]+)*)$"
        )
        errmsg = (
            "query error: decimal value or range expected for field " + field + "(got: "
        )
        data = defaultdict(set)
        conditions = []
        vallist = []

        for val in vals:
            val = str(val).strip()
            # single value
            if ":" not in val:
                match = op1.match(val)
                if not match:
                    sys.exit(errmsg + val + ")")
                op = match.group(2) if match.group(2) else defaultop
                operator = (
                    self.operators["negated"][op]
                    if match.group(1)
                    else self.operators["genuine"][op]
                )
                val = match.group(3)
                data[operator].add(val)

            # range
            else:
                match = op2.match(val)
                if not match:
                    sys.exit(errmsg + val + ")")
                operator = (
                    self.operators["negated"]["BETWEEN"]
                    if match.group(1)
                    else self.operators["genuine"]["BETWEEN"]
                )
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
        op1 = re.compile(r"^(\^*)((?:>|>=|<|<=|!=|=)?)([0-9]{4}-[0-9]{2}-[0-9]{2})$")
        op2 = re.compile(
            r"^(\^*)([0-9]{4}-[0-9]{2}-[0-9]{2}):([0-9]{4}-[0-9]{2}-[0-9]{2})$"
        )
        errmsg = (
            "query error: date or date range expected for field " + field + " (got: "
        )
        data = defaultdict(set)
        conditions = []
        vallist = []

        for val in vals:
            val = str(val).strip()
            # single date
            if ":" not in val:
                match = op1.match(val)
                if not match:
                    sys.exit(errmsg + val + ")")
                op = match.group(2) if match.group(2) else defaultop
                operator = (
                    self.operators["negated"][op]
                    if match.group(1)
                    else self.operators["genuine"][op]
                )
                val = match.group(3)
                data[operator].add(val)

            # date range
            else:
                match = op2.match(val)
                if not match:
                    sys.exit(errmsg + val + ")")
                operator = (
                    self.operators["negated"]["BETWEEN"]
                    if match.group(1)
                    else self.operators["genuine"]["BETWEEN"]
                )
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
        conditions = []
        vallist = []

        for val in vals:
            if val.startswith("^"):
                val = val[1:]
                opkey = "negated"
            else:
                opkey = "genuine"
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
        conditions = []
        vallist = []

        for val in vals:
            if val.startswith("^"):
                val = val[1:]
                opkey = "negated"
            else:
                opkey = "genuine"
            val = val.strip("%") + "%"
            operator = self.operators[opkey]["LIKE"]

            data[operator].add(val)

        for operator, valset in data.items():
            condition, vals = self.get_conditional_expr(field, operator, *valset)
            conditions.extend(condition)
            vallist.extend(vals)

        return link.join(conditions), vallist

    def query_metadata(self, name, *vals):
        conditions = ['"property_id" = ?']
        valueList = [self.properties[name]["id"]]
        targetfield = "value_" + self.properties[name]["datatype"]

        # query dates
        if self.properties[name]["querytype"] == "date":
            a, b = self.query_dates(targetfield, *vals)
            conditions.append(a)
            valueList.extend(b)

        # query numeric
        elif self.properties[name]["querytype"] == "numeric":
            a, b = self.query_numeric(targetfield, *vals)
            conditions.append(a)
            valueList.extend(b)

        # query text
        elif self.properties[name]["querytype"] == "text":
            a, b = self.query_string(targetfield, *vals)
            conditions.append(a)
            valueList.extend(b)

        # query zip
        elif self.properties[name]["querytype"] == "zip":
            a, b = self.query_zip(targetfield, *vals)
            conditions.append(a)
            valueList.extend(b)

        else:
            sys.exit(
                "error: unknown query type '"
                + self.properties[name]["querytype"]
                + "' for property '"
                + name
                + "'."
            )

        return (
            'SELECT "sample_id" AS id FROM sample2property WHERE '
            + " AND ".join(conditions),
            valueList,
        )

    def query_profile(self, *vars, reference_accession=None):  # noqa: C901
        iupac_nt_code = {
            "A": set("A"),
            "C": set("C"),
            "G": set("G"),
            "T": set("T"),
            "R": set("AGR"),
            "Y": set("CTY"),
            "S": set("GCS"),
            "W": set("ATW"),
            "K": set("GTK"),
            "M": set("ACM"),
            "B": set("CGTB"),
            "D": set("AGTD"),
            "H": set("ACTH"),
            "V": set("ACGV"),
            "N": set("ACGTRYSWKMBDHVN"),
        }
        iupac_aa_code = {
            "A": set("A"),
            "R": set("R"),
            "N": set("N"),
            "D": set("D"),
            "C": set("C"),
            "Q": set("Q"),
            "E": set("E"),
            "G": set("G"),
            "H": set("H"),
            "I": set("I"),
            "L": set("L"),
            "K": set("K"),
            "M": set("M"),
            "F": set("F"),
            "P": set("P"),
            "S": set("S"),
            "T": set("T"),
            "W": set("W"),
            "Y": set("Y"),
            "V": set("V"),
            "U": set("U"),
            "O": set("O"),
            "B": set("DNB"),
            "Z": set("EQZ"),
            "J": set("ILJ"),
            "Φ": set("VILFWYMΦ"),
            "Ω": set("FWYHΩ"),
            "Ψ": set("VILMΨ"),
            "π": set("PGASπ"),
            "ζ": set("STHNQEDKRζ"),
            "+": set("KRH+"),
            "-": set("DE-"),
            "X": set("ARNDCQEGHILKMFPSTWYVUOBZJΦΩΨπζ+-X"),
        }
        del_regex = re.compile(r"^(|[^:]+:)?([^:]+:)?del:(=?[0-9]+)(|-=?[0-9]+)?$")
        snv_regex = re.compile(r"^(|[^:]+:)?([^:]+:)?([A-Z]+)([0-9]+)(=?[A-Z]+)$")

        # set variants and generate sql
        base_sql = 'SELECT "sample.id" AS id FROM variantView WHERE '
        intersect_sqls = []
        intersect_vals = []
        except_sqls = []
        except_vals = []
        for var in vars:
            c = []  # condition
            v = []  # variable

            if var.startswith("^"):
                var = var[1:]
                negate = True
            else:
                negate = False

            # variant typing
            if match := snv_regex.match(var):
                snv = True
            elif match := del_regex.match(var):
                snv = False
            else:
                sys.exit("input error: " + var + " is not a valid variant definition.")

            # set molecule
            if match.group(1):
                c.append('"molecule.symbol" = ?')
                v.append(match.group(1)[:-1])
            else:
                c.append('"molecule.standard" = ?')
                v.append(1)

            # set element
            if match.group(2):
                c.append('"element.type" = ?')
                v.append("cds")
                c.append('"element.symbol" = ?')
                v.append(match.group(2)[:-1])
                code = iupac_aa_code
            else:
                c.append('"element.standard" = ?')
                v.append(1)
                code = iupac_nt_code

            # snp and insert handling
            if snv:
                c.append('"variant.start" = ?')
                v.append(int(match.group(4)) - 1)
                c.append('"variant.end" = ?')
                v.append(int(match.group(4)))
                c.append('"variant.ref" = ?')
                v.append(match.group(3))

                # explicit alternate allele
                if match.group(5).startswith("="):
                    c.append('"variant.alt" = ?')
                    v.append(match.group(5)[1:])

                # potentially ambiguous alternate snp
                elif len(match.group(5)) == 1:
                    l = len(code[match.group(5)])
                    if l == 1:
                        c.append('"variant.alt" = ?')
                        v.append(match.group(5))
                    else:
                        c.append('"variant.alt" IN (' + ", ".join(["?"] * l) + ")")
                        v.extend(code[match.group(5)])

                # potentially ambiguous alternate insert
                else:
                    a = [
                        "".join(x)
                        for x in itertools.product(*[code[x] for x in match.group(5)])
                    ]
                    l = len(a)
                    if l == 1:
                        c.append('"variant.alt" = ?')
                        v.extend(a)
                    else:
                        c.append('"variant.alt" IN (' + ", ".join(["?"] * l) + ")")
                        v.extend(a)

            # deletion handling
            else:
                s = match.group(3)
                e = match.group(4)[1:]

                if s.startswith("="):
                    s = s[1:]
                    c.append('"variant.start" = ?')
                else:
                    c.append('"variant.start" <= ?')

                if e.startswith("="):
                    e = e[1:]
                    c.append('"variant.end" = ?')
                else:
                    c.append('"variant.end" >= ?')
                v.append(int(s) - 1)
                v.append(int(e))

                c.append('"variant.alt" = ?')
                v.append(" ")

            # assemble sub-sql
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

    def count_molecules(self, reference_accession=None):
        if reference_accession:
            condition = '"reference.accession" = ?'
            vals = [reference_accession]
        else:
            condition = '"reference.standard" = ?'
            vals = [1]
        sql = (
            'SELECT count(DISTINCT "molecule.id") AS count FROM referenceView WHERE '
            + condition
        )
        return self.cursor.execute(sql, vals).fetchone()["count"]

    def match(  # noqa: 901
        self,
        *profiles,
        reserved_props=None,
        properties=None,
        reference_accession=None,
        showN=False,
        format="csv",
    ):
        # collecting sqls for metadata-based filtering
        property_sqls = []
        property_vals = []
        # IF sublineage search is enable
        # support: include and exclude
        if "with_sublineage" in reserved_props:
            _tmp_include_lin = []  # used to keep all lienages after search.
            lineage_col = reserved_props.get("with_sublineage")
            include_lin = properties.get(lineage_col)  # get list of given lineages
            negate = False
            print(include_lin)
            while include_lin:
                in_lin = include_lin.pop(0)

                if in_lin.startswith("^"):
                    in_lin = in_lin[1:]
                    negate = True

                # have wildcard in string which mean we have to find all lineage from wildcard query
                # then we used the wildcard query result to find all sublineages agian.
                if "%" in in_lin:
                    _tobeadded_lin = self.get_list_of_lineages(in_lin)
                    for i in _tobeadded_lin:

                        # if i != in_lin: # we dont need to add same lineage agian,so we skip for the duplicate lineage.
                        if negate:  # all lineage should add not ^
                            i = "^" + i
                        include_lin.append(i)  # add more lineage to find in next round.

                value = self.lineage_sublineage_dict.get(
                    in_lin, "none"
                )  # provide a default value if the key is missing:
                # print(value)
                if value != "none":
                    if negate:
                        in_lin = "^" + in_lin
                    _tmp_include_lin.append(in_lin)

                    _list = value.split(",")
                    for i in _list:
                        if negate:  # all sublineage should add not^
                            i = "^" + i
                        include_lin.append(i)  # add more lineage to find in next round.
                        # _tmp_include_lin.append(i)
                        # if we don't find this wildcard so we discard it
                else:  # None (no child)
                    if negate:
                        in_lin = "^" + in_lin
                    _tmp_include_lin.append(in_lin)
                negate = False

            include_lin = _tmp_include_lin
            properties[lineage_col] = include_lin
        # print(properties)

        if properties:
            for pname, vals in properties.items():
                sql, val = self.query_metadata(pname, *vals)
                property_sqls.append(sql)
                property_vals.extend(val)

        property_sqls = " INTERSECT ".join(property_sqls)

        # collecting sqls for genomic profile based filtering
        profile_sqls = []
        profile_vals = []
        for profile in profiles:
            sql, val = self.query_profile(
                *profile, reference_accession=reference_accession
            )
            profile_sqls.append(sql)
            profile_vals.extend(val)

        if len(profiles) == 1:
            profile_sqls = profile_sqls[0]
        elif len(profiles) > 1:
            profile_sqls = " UNION ".join(
                ["SELECT * FROM (" + x + ")" for x in profile_sqls]
            )
        else:
            profile_sqls = ""

        if property_sqls and profile_sqls:
            if len(profiles) > 1:
                sample_selection_sql = (
                    property_sqls + " INTERSECT SELECT * FROM (" + profile_sqls + ")"
                )
            else:
                sample_selection_sql = property_sqls + " INTERSECT " + profile_sqls
        elif property_sqls or profile_sqls:
            sample_selection_sql = property_sqls + profile_sqls
        else:
            if "sample" in reserved_props:
                # if 'sample' is presented we just use only samples
                samples_condition = []
                for pname, vals in reserved_props.items():
                    if pname == "sample":
                        for x in vals:
                            samples_condition.append('"' + x + '"')
                sample_selection_sql = (
                    "SELECT id FROM sample WHERE name IN ("
                    + " , ".join(samples_condition)
                    + ")"
                )

                property_sqls = []
                property_vals = []
            else:
                sample_selection_sql = "SELECT id FROM sample"

        genome_element_condition = [
            str(x) for x in self.get_element_ids(reference_accession, "source")
        ]
        if len(genome_element_condition) == 1:
            genome_element_condition = '"element.id" = ' + genome_element_condition[0]
            m = ""
        else:
            genome_element_condition = (
                '"element.id" IN (' + ", ".join(genome_element_condition) + ")"
            )
            m = ' "molecule.symbol" || "@" || '

        if not showN:
            nn = ' AND "variant.alt" != "N" '
            nx = ' AND "variant.alt" != "X" '
        else:
            nn = ""
            nx = ""

        cds_element_condition = [
            str(x) for x in self.get_element_ids(reference_accession, "cds")
        ]
        if len(cds_element_condition) == 1:
            cds_element_condition = '"element.id" = ' + cds_element_condition[0]
        else:
            cds_element_condition = (
                '"element.id" IN (' + ", ".join(cds_element_condition) + ")"
            )

        # standard query
        if format == "csv" or format == "tsv":

            # select samples
            sql = sample_selection_sql
            sample_ids = self.cursor.execute(
                sql, property_vals + profile_vals
            ).fetchall()

            if not sample_ids:
                return []
            selected_sample_ids = ", ".join([str(x["id"]) for x in sample_ids])
            # rows = {x['id']: {"id": x['id']} for x in sample_ids}
            # print(len(sample_ids))

            #
            # Current solution:
            # After we got the selected IDs
            # We use two-stage query and then combine both results together to produce final result
            # 1. Query properties
            # 2. Query  AA/NT profile
            fields = ['"sample.name"'] + ['"' + x + '"' for x in self.properties]
            sql = "SELECT name as " + ", ".join(fields) + "FROM sample "

            joins = [
                "LEFT JOIN (SELECT sample_id, value_"
                + y["datatype"]
                + " as "
                + x
                + " FROM sample2property WHERE property_id = "
                + str(y["id"])
                + ") as t"
                + str(y["id"])
                + " ON sample.id = t"
                + str(y["id"])
                + ".sample_id"
                for x, y in self.properties.items()
            ]
            _1_final_sql = (
                sql
                + " ".join(joins)
                + " WHERE sample.id IN ("
                + selected_sample_ids
                + ")"
            )
            # print(_1_final_sql)
            _1_rows = self.cursor.execute(_1_final_sql).fetchall()
            _2_final_sql = (
                " SELECT name AS 'sample.name'  , nt_profile._profile AS NUC_PROFILE, aa_profile._profile AS AA_PROFILE \
                    FROM \
                            ( \
                              SELECT  \"sample.id\", group_concat("
                + m
                + '"variant.label") AS _profile \
                              FROM variantView WHERE "sample.id" IN ('
                + selected_sample_ids
                + ") AND "
                + genome_element_condition
                + nn
                + ' GROUP BY "sample.id" \
                            ) nt_profile, \
                            ( \
                              SELECT  "sample.id", group_concat('
                + m
                + '"element.symbol" || ":" || "variant.label") AS _profile \
                              FROM variantView WHERE "sample.id" IN ('
                + selected_sample_ids
                + ") AND "
                + cds_element_condition
                + nx
                + ' GROUP BY "sample.id" \
                            ) aa_profile, \
                            sample \
                    WHERE nt_profile."sample.id" = aa_profile."sample.id" AND  nt_profile."sample.id" = sample.id  \
                          AND sample.id  IN ('
                + selected_sample_ids
                + ")"
            )
            _2_rows = self.cursor.execute(_2_final_sql).fetchall()
            if len(_1_rows) != len(_2_rows):
                print("WARNING: There are something error between query")
            # print(set([ x['sample.name'] for x in _1_rows ]) ^ set([ x['name'] for x in _2_rows ]))
            # To combine:
            # We update list of dict (update on result from query #2)
            # merge all results
            _1_rows.extend(
                list(
                    map(
                        lambda x, y: x.update(
                            {
                                key: value
                                for key, value in y.items()
                                if (key == "NUC_PROFILE") or (key == "AA_PROFILE")
                            }
                        )
                        if x.get("sample.name") == y.get("sample.name")
                        else None,
                        _1_rows,
                        _2_rows,
                    )
                )
            )

            _1_rows = list(filter(None, _1_rows))

            # since we use "update" function (i.e. extends the dict. to include all key:value from properties base on sample name)
            # at _1_rows so we can return _1_rows only a
            return _1_rows  # list(rows.values())
            """
            sql = "WITH selected_samples AS (" + sample_selection_sql + ") \
               SELECT  *, \
                        ( \
                          SELECT group_concat(" + m + "\"variant.label\") AS nuc_profile \
                          FROM variantView WHERE \"sample.id\" IN (SELECT id FROM selected_samples) AND " + genome_element_condition + nn + " GROUP BY \"sample.name\" ORDER BY \"element.id\", \"variant.start\" \
                        ) nt_profile, \
                        ( \
                          SELECT group_concat(" + m + "\"element.symbol\" || \":\" || \"variant.label\") AS nuc_profile \
                          FROM variantView WHERE \"sample.id\" IN (SELECT id FROM selected_samples) AND " + cds_element_condition + np + " GROUP BY \"sample.name\" ORDER BY \"element.id\", \"variant.start\" \
                        ) aa_profile \
                FROM sample \
                WHERE id IN ( SELECT id FROM selected_samples )"


            sql = "SELECT * \
                   FROM variantView \
                   WHERE \"sample.id\" IN (" + sample_selection_sql + ")" # ORDER BY \"sample.id\", \"element.id\", \"variant.start\""
            """
        elif format == "count":
            sql = (
                "SELECT count(DISTINCT id) as count FROM (" + sample_selection_sql + ")"
            )
        elif format == "vcf":
            sql = (
                'SELECT "element.id",  "element.type", "molecule.accession", "variant.start", "variant.ref", "variant.alt", "variant.label", group_concat("sample.name") as samples FROM variantView WHERE "sample.id" IN ('
                + sample_selection_sql
                + ") AND "
                + genome_element_condition
                + nn
                + 'GROUP BY "molecule.accession", "variant.start", "variant.ref", "variant.alt" ORDER BY "molecule.accession", "variant.start"'
            )
        else:
            sys.exit("error: '" + format + "' is not a valid output format")
        return self.cursor.execute(sql, property_vals + profile_vals)  # .fetchall()

    def get_list_of_lineages(self, lineage):
        sql = (
            "SELECT DISTINCT lineage FROM lineages WHERE lineage LIKE '"
            + lineage
            + "';"
        )
        rows = self.cursor.execute(sql).fetchall()
        result = [i["lineage"] for i in rows]
        return result

    # MISC

    @staticmethod
    def optimize(dbfile):
        with sqlite3.connect(dbfile) as con:
            con.executescript("PRAGMA analysis_limit=400;")
            con.executescript("PRAGMA optimize;")
            con.executescript("VACUUM;")

    @staticmethod
    def dict_factory(cursor, row):
        d = {}  # OrderedDict()
        for idx, col in enumerate(cursor.description):
            d[col[0]] = row[idx]
        return d

    # Utils.
    def get_db_size(self, decimal_places=3):
        size = os.path.getsize(self.dbfile)
        for unit in ["B", "KiB", "MiB", "GiB", "TiB"]:
            if size < 1024.0:
                break
            size /= 1024.0
        return f"{size:.{decimal_places}f}{unit}"

    @staticmethod
    def upgrade_db(dbfile):
        try:
            with sqlite3.connect(dbfile) as con:
                cur = con.cursor()
                current_version = cur.execute("pragma user_version").fetchone()[0]

            print(
                "Current version:",
                current_version,
                " Upgrade to:",
                SUPPORTED_DB_VERSION,
            )
            uri = "file:" + urlquote(dbfile)
            print("Perform the Upgrade:", uri)
            while current_version < SUPPORTED_DB_VERSION:

                next_version = current_version + 1
                file_path = os.path.join(
                    os.path.dirname(os.path.realpath(__file__)),
                    "migrate/" + str(next_version) + ".sql",
                )
                if not os.path.isfile(file_path):
                    raise ValueError("Sorry, we cannot find", file_path)

                with open(file_path, "r") as handle:
                    sql = handle.read()
                with sqlite3.connect(uri + "?mode=rwc", uri=True) as con:
                    con.executescript(sql)

                current_version = next_version

        except sqlite3.Error as er:
            con.executescript("ROLLBACK")
            raise er
        except ValueError as er:
            print(er)
        finally:
            print("Database now version:", current_version)
            if current_version == SUPPORTED_DB_VERSION:
                print("Success: Database upgrade was successfully completed")
            else:
                print("Error: Upgrade was not completed")
