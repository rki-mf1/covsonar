from collections import defaultdict
from datetime import datetime
import itertools
import os
import pkgutil
import re
import sqlite3
import sys
from typing import Any
from typing import Dict
from typing import Iterator
from typing import List
from typing import Optional
from typing import Set
from typing import Tuple
from typing import Union
from urllib.parse import quote as urlquote

from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation
from Bio.SeqFeature import FeatureLocation
import pandas as pd
import sqlparse

from covsonar.logging import LoggingConfigurator


# Constants
SUPPORTED_DB_VERSION = 6  # supprted sonar database scheme version

# Initialize logger
LOGGER = LoggingConfigurator.get_logger()


class sonarDbManager:
    """
    A class to handle genomic data stored in SQLite database.

    Public attributes
        dbfile (str): SQLite database file path.
        conn (sqlite3.Connection): SQLite database connection object.
        cursor (sqlite3.Cursor): SQLite database cursor for executing SQL queries.
    """

    # CONSTANTS

    OPERATORS = {
        "standard": {
            "=": "=",
            ">": ">",
            "<": "<",
            ">=": ">=",
            "<=": "<=",
            "IN": "IN",
            "LIKE": "LIKE",
            "BETWEEN": "BETWEEN",
        },
        "inverse": {
            "=": "!=",
            ">": "<=",
            "<": ">=",
            ">=": "<",
            "<=": ">",
            "IN": "NOT IN",
            "LIKE": "NOT LIKE",
            "BETWEEN": "NOT BETWEEN",
        },
        "default": "=",
    }

    IUPAC_CODES = {
        "nt": {
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
            "n": set("N"),
        },
        "aa": {
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
            "x": set("X"),
        },
    }

    def __init__(self, dbfile: str, readonly: bool = True) -> None:
        """
        Initialize the DBOperations class.

        Args:
            dbfile (str): The path of the SQLite database file.
            readonly: If true, database connection is read-only.
        """
        # check dbfile
        if not os.path.isfile(dbfile):
            LOGGER.error("The database does not exist.")
            sys.exit(1)

        # public attributes
        self.dbfile = os.path.abspath(dbfile)
        self.con = None
        self.cursor = None

        # private attributes
        self.__timeout = -1
        self.__mode = "ro" if readonly else "rwc"
        self.__uri = self.get_uri(dbfile)
        self.__illegal_properties = {"SAMPLE"}
        self.__lineage_sublineage_dict = None

        self.__properties = None

    def __enter__(self) -> "sonarDbManager":
        """
        Enter the runtime context related to the database object.

        Returns:
            DBOperations: The current instance.
        """
        self.con = self.connect()
        self.cursor = self.con.cursor()
        self.check_db_compatibility()
        self.start_transaction()
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        """
        Exit the runtime context and close the database connection.
        In case of raised errors, the database is rolled back.
        """
        if [exc_type, exc_value, exc_traceback].count(None) != 3:
            if self.__mode == "rwc":
                LOGGER.info("rollback database")
                self.rollback()
        elif self.__mode == "rwc":
            self.commit()
        self.close()

    # PROPERTIES

    @property
    def count(self) -> int:
        """Provides the count of the rows in the 'sample' table.

        Returns:
            int: Count of rows in the 'sample' table.
        """
        return self.cursor.execute("SELECT COUNT(*) as count FROM sample").fetchone()[
            "count"
        ]

    @property
    def seq_count(self) -> int:
        """Provides the count of the sequences in the 'sequence' table.

        Returns:
            int: Count of sequences in the 'sequence' table.
        """
        return self.cursor.execute(
            "SELECT COUNT(DISTINCT seqhash) as count FROM sequence"
        ).fetchone()["count"]

    @property
    def variant_count(self) -> int:
        """Provides the count of the variants in the 'variant' table.

        Returns:
            int: Count of variants in the 'variant' table.
        """
        return self.cursor.execute("SELECT COUNT(*) as count FROM variant").fetchone()[
            "count"
        ]

    @property
    def reference_count(self) -> int:
        """Provides the count of the references in the 'reference' table.

        Returns:
            int: Count of references in the 'reference' table.
        """
        return self.cursor.execute(
            "SELECT COUNT(*) as count FROM reference"
        ).fetchone()["count"]

    @property
    def properties(self) -> Dict[str, Dict[str, Any]]:
        """
        Returns property data as a dict of dict where key is property name.
        If data is not in the cache, it fetches data from the SQLite database.

        Returns:
            Dict[str, Dict[str, Any]]: A dictionary with property names as keys
            and corresponding property data as values.
        """
        if not self.__properties:
            sql = "SELECT * FROM property_metadata;"
            rows = self.cursor.execute(sql).fetchall()
            self.__properties = {} if not rows else {row["name"]: row for row in rows}
        return self.__properties

    @property
    def lineage_sublineage_dict(self) -> Dict[str, str]:
        """
        Property that returns a dictionary mapping lineage to sublineage.
        The dictionary is created based on data read from a SQL query.
        The dictionary is cached for future use, i.e., the SQL query is executed only the first time this property is accessed.

        Returns:
            dict: A dictionary where the keys are lineage and the values are sublineage.
        """
        if not self.__lineage_sublineage_dict:
            sql = "SELECT lineage, sublineage FROM lineages"
            rows = self.cursor.execute(sql).fetchall()
            self.__lineage_sublineage_dict = (
                {} if not rows else {row["lineage"]: row["sublineage"] for row in rows}
            )
        return self.__lineage_sublineage_dict

    # BASIC OPERATIONS

    @staticmethod
    def setup(filename: str) -> None:
        """
        Set up the sonar database.

        Args:
            filename (str): Path to the SQLite database file.

        Returns:
            None

        Raises:
            FileNotFoundError: Raised if the .sqlite file does not exist.
            sqlite3.Error: Raised if there is an error connecting to the SQLite database.

        Example:
            >>> dbfile = getfixture('tmpfile_name')
            >>> sonarDbManager.setup(dbfile)

        """
        try:
            # Get SQL schema from .sqlite file in the package's data directory.
            sql_schema = pkgutil.get_data(__name__, "data/db.sqlite").decode()
        except FileNotFoundError as e:
            LOGGER.error(f"No .sqlite file found in data directory: {str(e)}")
            raise e

        # Get URI for SQLite database.
        uri = sonarDbManager.get_uri(filename)

        # Connect to the SQLite database and set up the schema.
        try:
            with sqlite3.connect(f"{uri}?mode=rwc", uri=True) as con:
                con.set_trace_callback(LOGGER.debug)
                con.executescript(sql_schema)
        except sqlite3.Error as e:
            LOGGER.error(f"Error setting up SQLite database: {str(e)}")
            raise e

    @staticmethod
    def get_uri(fname) -> str:
        """returns url-safe file uri as string

        Parameters:
            fname (str): The file name.

        Returns:
            str: The URI of the database file.

        Examples:
            >>> sonarDbManager.get_uri("test.db")
            'file:test.db'
        """
        return "file:" + urlquote(fname)

    def connect(self) -> "sqlite3.Connection":
        """
        Create a connection to the SQLite database and set the row factory.

        Returns:
            sqlite3.Connection: The connection object for the SQLite database.
        Raises
            SystemExit
        """
        con = sqlite3.connect(
            self.__uri + "?mode=" + self.__mode,
            self.__timeout,
            isolation_level=None,
            uri=True,
        )
        con.set_trace_callback(LOGGER.debug)
        con.row_factory = self.dict_factory
        return con

    def start_transaction(self) -> None:
        """Start a new transaction."""
        self.cursor.execute("BEGIN DEFERRED")

    def commit(self) -> None:
        """
        Commit the current transaction.

        This method saves any changes made since the start of the transaction
        to the SQLite database.
        """
        self.con.commit()

    def rollback(self) -> None:
        """
        Rollback the current transaction.

        This method undoes any changes made since the start of the transaction.
        """
        self.con.rollback()

    def close(self) -> None:
        """
        Close the connection to the SQLite database.

        After this method is called, the database connection will be closed
        and no further SQL statements can be executed.
        """
        self.con.close()

    def get_db_version(self) -> int:
        """
        Get the version number of the SQLite database.

        Returns:
            int: The version number of the SQLite database.
        """
        return self.cursor.execute("pragma user_version").fetchone()["user_version"]

    @staticmethod
    def optimize(dbfile: str, tmpdir: Optional[str] = None) -> None:
        """
        Optimize the SQLite database file for better performance.

        Args
            dbfile (str): SQLite database file path.
            tmpdir (str): directory to use as temp directory, optional.
        """
        # perform cleaning
        LOGGER.info("cleaning ...")
        with sonarDbManager(dbfile, readonly=False) as dbm:
            dbm.clean()

        # set tmp dir to folder where the database is stored
        if tmpdir:
            os.environ["TMPDIR"] = tmpdir

        # perform vacuum
        LOGGER.info("optimizing ...")
        with sqlite3.connect(dbfile) as con:
            con.executescript("PRAGMA analysis_limit=400;")
            con.executescript("PRAGMA optimize;")
            con.executescript("VACUUM;")
            con.commit()

    def clean(self) -> None:
        """
        Clean the SQLite database by removing data from tables where certain conditions are not met.
        """
        # SQL queries for cleanup
        cleanup_queries = [
            "DELETE FROM sequence WHERE NOT EXISTS(SELECT NULL FROM sample WHERE sample.seqhash = seqhash)",
            "DELETE FROM translation WHERE NOT EXISTS(SELECT NULL FROM reference WHERE reference.translation_id = id)",
            "DELETE FROM molecule WHERE NOT EXISTS(SELECT NULL FROM reference WHERE reference.id = reference_id)",
            "DELETE FROM element WHERE NOT EXISTS(SELECT NULL FROM molecule WHERE molecule.id = molecule_id)",
            "DELETE FROM elempart WHERE NOT EXISTS(SELECT NULL FROM element WHERE element.id = element_id)",
            "DELETE FROM variant WHERE NOT EXISTS(SELECT NULL FROM alignment2variant WHERE alignment2variant.variant_id = variant_id)",
            "DELETE FROM alignment WHERE NOT EXISTS(SELECT NULL FROM sequence WHERE sequence.seqhash = seqhash) OR NOT EXISTS(SELECT NULL FROM element WHERE element.id = element_id)",
            "DELETE FROM alignment2variant WHERE NOT EXISTS(SELECT NULL FROM alignment WHERE alignment.id = alignment_id)",
        ]

        # Execute each cleanup query
        for query in cleanup_queries:
            self.cursor.execute(query)

    @staticmethod
    def format_sql(sql: str) -> str:
        """
        Formats the given SQL string.

        Args:
            sql (str): The raw SQL string.

        Returns:
            str: The formatted SQL string.
        """
        return sqlparse.format(sql, reindent=True, keyword_case="upper")

    @staticmethod
    def is_select_query(query):
        """
        Check if the input query string consists of only SELECT statements.

        Args:
            query (str): Input SQL query string.

        Returns:
            bool: True if the input string is a SELECT query (or multiple SELECT queries), False otherwise.

        Examples:
            >>> sonarDbManager.is_select_query("SELECT * FROM table1; SELECT * FROM table2")
            True
            >>> sonarDbManager.is_select_query("SELECT * FROM table1; UPDATE table2 SET field1='value'")
            False
        """
        parsed_queries = sqlparse.split(
            query
        )  # splits multi-queries into list of individual queries)

        for query in parsed_queries:
            parsed = sqlparse.parse(query)  # parse the SQL query

            if not parsed:  # if the parsed result is empty
                return False

            for statement in parsed:
                if (
                    not statement.get_type() == "SELECT"
                ):  # check if the statement type is SELECT
                    return False

        return True

    @staticmethod
    def dict_factory(cursor: sqlite3.Cursor, row: Tuple) -> Dict:
        """
        Convert cursor rows into dictionary format.

        Args:
            cursor (sqlite3.Cursor): SQLite cursor object.
            row (Tuple): A row of data from the SQLite cursor.

        Returns
            Dict: Row data in dictionary format.
        """
        return {col[0]: row[idx] for idx, col in enumerate(cursor.description)}

    # VERSION & UPGRADES

    def get_db_size(self, decimal_places: int = 3) -> str:
        """
        Get the size of the SQLite database file in a human-readable format.

        Args:
            decimal_places (int, optional): Number of decimal places to use in the size. Defaults to 3.

        Returns:
            str: Human-readable size of the SQLite database file.
        """
        size = os.path.getsize(self.dbfile)
        for unit in ["B", "KiB", "MiB", "GiB", "TiB"]:
            if size < 1024.0:
                break
            size /= 1024.0
        return f"{size:.{decimal_places}f}{unit}"

    def check_db_compatibility(self) -> None:
        """
        Check the compatibility of the SQLite database.

        This method checks whether the version of the database is compatible with
        the software. It compares the current version of the database with the
        supported version defined in the SUPPORTED_DB_VERSION variable.

        Raises:
            SystemExit: If the database version is not identical to the supported version,
                        indicating that the software might be outdated or too new.
        """
        current_version = self.get_db_version()
        if not current_version == SUPPORTED_DB_VERSION:
            LOGGER.error(
                "The given database is not compatible with this version of sonar (database version: "
                + str(current_version)
                + "; supported database version: "
                + str(SUPPORTED_DB_VERSION)
                + ")"
            )
            sys.exit(1)

    @staticmethod
    def upgrade_db(dbfile: str) -> None:
        """
        Upgrades the database to the supported version. It executes the migration scripts
        located in the `migrate` directory, each script corresponds to a version upgrade.
        The scripts must be named as the version number and must be .sql files.

        Args:
            dbfile (str): Path to the database file.

        Raises:
            sqlite3.Error: If there's an error in SQL execution.
            ValueError: If the required migration script file is not found.
            RuntimeError: If the database version does not match the supported version after upgrade.
        """
        with sqlite3.connect(dbfile) as con:
            cur = con.cursor()
            current_version = cur.execute("pragma user_version").fetchone()[0]

        if current_version < 5:
            LOGGER.error(
                f"Sorry, but automated migration does not support databases of version {current_version}."
            )
            sys.exit(1)

        try:
            LOGGER.info(
                f"Current version: {current_version}, Upgrade to: {SUPPORTED_DB_VERSION}"
            )
            uri = "file:" + urlquote(dbfile)
            LOGGER.info(f"Perform the Upgrade: {uri}")

            while current_version < SUPPORTED_DB_VERSION:
                next_version = current_version + 1
                file_path = os.path.join(
                    os.path.dirname(os.path.realpath(__file__)),
                    "migrate",
                    f"{next_version}.sql",
                )

                if not os.path.isfile(file_path):
                    raise ValueError(
                        f"error: cannot find {file_path}, please contact us or reinstall software."
                    )

                with open(file_path, "r") as handle:
                    sql = handle.read()

                with sqlite3.connect(uri + "?mode=rwc", uri=True) as con:
                    con.executescript(sql)

                current_version = next_version

        except sqlite3.Error as er:
            con.executescript("ROLLBACK")
            raise er
        except ValueError as er:
            LOGGER.error(er)
            raise er
        finally:
            LOGGER.info(f"Database now version: {current_version}")

            if current_version == SUPPORTED_DB_VERSION:
                LOGGER.info("Success: Database upgrade was successfully completed.")
            else:
                raise RuntimeError("Error: Upgrade was not completed.")

    # QUERIES

    def direct_query(self, sql: str) -> List[Dict]:
        """
        Perform a direct SQL query on the database and return the results.

        Args:
            sql (str): SQL query to be executed.

        Returns:
            List[Dict]: A list of dictionaries representing each row of data from the executed SQL query.
        """
        # remove quotes around the query if present
        if not sonarDbManager.is_select_query(sql):
            LOGGER.error("Only SELECT statements are allowed as direct query.")
            sys.exit(1)
        return self.cursor.execute(sql).fetchall()

    # IMPORT NON-SAMPLE DATA

    def add_codon(self, translation_table: int, codon: str, amino_acid: str) -> None:
        """Adds a codon to the database.

        This method will add a codon to the database, given the name, amino acid,
        and bases that make up the codon.

        Args:
            translation_table (int): The tarbnslation table associated.
            codon (str): The codon nucleotides.
            amino_acid (str): The corresponding amino acid.

        Returns:
            None
        """
        sql = "INSERT OR IGNORE INTO translation (id, codon, aa) VALUES(?, ?, ?);"
        self.cursor.execute(sql, [translation_table, codon, amino_acid])

    def add_property(
        self,
        name: str,
        datatype: str,
        querytype: str,
        description: str,
        subject: str,
        standard: Optional[str] = None,
        check_name: bool = True,
    ) -> int:
        """
        Adds a new property and returns the property id.

        Args:
            name (str): The name of the property.
            datatype (str): The data type of the property.
            querytype (str): The query type of the property.
            description (str): The description of the property.
            subject (str): The subject of the property.
            standard (Optional[str], optional): The standard of the property. Defaults to None.
            check_name (bool, optional): Whether to check the property name. Defaults to True.

        Returns:
            int: The id of the added property.

        Raises:
            SystemExit: If an illegal or duplicate property name is used, or if there is a failure to insert data into sqlite table.

        Examples:
            >>> dbm = getfixture('init_writeable_dbm')
            >>> id = dbm.add_property("NEW_PROP", "text", "text", "my new prop stores text information", "sample")
        """
        name = name.upper()
        if name in self.__illegal_properties:
            LOGGER.error(
                f"Name '{name}' is reserved and cannot be used as property name"
            )
            sys.exit(1)

        if check_name and not re.match("^[A-Z][A-Z0-9_]+$", name):
            LOGGER.error(
                "Invalid property name (property names have to start with an letter and can contain only letters, numbers and underscores)"
            )
            sys.exit(1)

        if name in self.properties:
            LOGGER.error(
                f"A property named '{name}' already exists in the given database."
            )
            sys.exit(1)

        try:
            # TODO: set default if it is defined
            sql = f"ALTER TABLE property ADD {name} {datatype};"
            self.cursor.execute(sql)
            sql = f"CREATE INDEX IF NOT EXISTS idx_property_{name} ON property({name});"
            self.cursor.execute(sql)
            sql = "INSERT INTO property_metadata (name, datatype, querytype, description, target, standard) VALUES(?, ?, ?, ?, ?, ?);"
            self.cursor.execute(
                sql, [name, datatype, querytype, description, subject, standard]
            )
            # FIXME: Not sure if we need this id anymore
            pid = 1
            self.__properties = False

        except sqlite3.Error as error:
            LOGGER.error(f"Failed to insert data into sqlite table ({str(error)}).")
            sys.exit(1)
        return pid

    def add_translation_table(self, translation_table: int) -> None:
        """
        Adds codon amino acid relationship for a given translation table to database.
        None-sense codons including gaps are assigned to a 0-length string.

        Args:
            translation_table (int): The translation table to be added to the database.

        Raises:
            sqlite3.Error: If an SQLite error occurs.

        Example usage:
            >>> dbm = getfixture('init_writeable_dbm')
            >>> dbm.add_translation_table(1)
        """
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
        self,
        accession: str,
        description: str,
        organism: str,
        translation_table: int,
        standard: int = 0,
    ) -> int:
        """
        Adds a reference to a database and returns the assigned row id.
        None-sense codons including gaps are assigned to a 0-length string.
        The reference is set to the default reference if standard is 1.

        Args:
            accession (str): The accession of the reference.
            description (str): The description of the reference.
            organism (str): The organism of the reference.
            translation_table (int): The translation table of the reference.
            standard (int, optional): Whether the reference is standard. Defaults to 0.

        Returns:
            int: The row ID of the newly added reference.

        Raises:
            sqlite3.Error: If an error occurs while interacting with the database.

        Example usage:
            >>> dbm = getfixture('init_writeable_dbm')
            >>> rowid = dbm.add_reference("REF1", "my new reference", "virus X", 1)
        """
        self.add_translation_table(translation_table)
        if standard:
            sql = "UPDATE reference SET standard = 0 WHERE standard = 1"
            self.cursor.execute(sql)
        sql = "INSERT INTO reference (id, accession, description, organism, translation_id, standard) VALUES(?, ?, ?, ?, ?, ?);"
        self.cursor.execute(
            sql, [None, accession, description, organism, translation_table, standard]
        )
        return self.cursor.lastrowid

    def add_update_lineage(self, lineage_df: pd.DataFrame) -> None:
        """
        Updates the lineages table in the database with given DataFrame.

        Args:
            lineage_df (pd.DataFrame): DataFrame with columns 'lineage' and 'sublineage'.
        """
        LOGGER.info(f"Prepare: {len(lineage_df)}")
        sql = "INSERT OR REPLACE INTO lineages (lineage, sublineage) VALUES (?, ?);"
        data = list(zip(lineage_df["lineage"], lineage_df["sublineage"]))
        self.cursor.executemany(sql, data)

    # IMPORT SAMPLE DATA

    def insert_property(
        self, sample_id: int, property_name: str, property_value: Union[str, int, float]
    ) -> None:
        """
        Inserts/Updates a property value of a given sample in the database.

        Args:
            sample_id (int): The ID of the sample for which the property is being updated.
            property_name (str): The name of the property being updated.
            property_value (Union[str, int, float]): The value of the property being updated.

        Raises:
            sqlite3.Error: If an error occurs while interacting with the database.

        Example usage:
            >>> dbm = getfixture('init_writeable_dbm')
            >>> dbm.insert_property(1, "LINEAGE", "BA.5")
        """
        illegal = {
            "GENOMIC_PROFILE",
            "SAMPLE_NAME",
            "PROTEOMIC_PROFILE",
            "FRAMESHIFT_MUTATION",
        }
        if property_name in illegal:
            LOGGER.error("This proprty name is reserved and cannot be used.")
            sys.exit(1)
        # Make sure a row with this sample_id exists in the property table
        sql = f"INSERT INTO property(sample_id,{property_name}) VALUES ({sample_id},'{property_value}') ON CONFLICT(sample_id) DO UPDATE SET {property_name}='{property_value}';"
        self.cursor.execute(sql)
        # self.cursor.execute(
        #     sql, [property_name, property_value, sample_id]
        # )

    def insert_sequence(self, seqhash: str) -> None:
        """
        Inserts a sequence represented by its hash to the database. If the hash is already known, it is ignored.

        Args:
            seqhash (str): The hash of the sequence to be inserted into the database.

        Raises:
            sqlite3.Error: If an error occurs while interacting with the database.

        Example usage:
            >>> dbm = getfixture('init_writeable_dbm')
            >>> dbm.insert_sequence("1a1f34ef4318911c2f98a7a1d6b7e9217c4ae1d1")
        """
        sql = "INSERT OR IGNORE INTO sequence (seqhash) VALUES(?);"
        self.cursor.execute(sql, [seqhash])

    def insert_sample(self, sample_name: str, seqhash: str) -> int:
        """
        Inserts or updates a sample/genome in the database and returns the sample id.

        Args:
            sample_name (str): The name of the sample to be inserted or updated in the database.
            seqhash (str): The hash of the sequence to be associated with the sample.

        Returns:
            int: The sample id of the inserted or updated sample.

        Raises:
            sqlite3.Error: If an error occurs while interacting with the database.

        Example usage:
            >>> dbm = getfixture('init_writeable_dbm')
            >>> rowid = dbm.insert_sample("my_new_sample", "1a1f34ef4318911c2f98a7a1d6b7e9217c4ae1d1")
        """
        sql = "INSERT OR REPLACE INTO sample (name, seqhash, datahash) VALUES(?, ?, ?);"
        self.cursor.execute(sql, [sample_name, seqhash, ""])
        sql = "SELECT id FROM sample WHERE name = ?;"
        sid = self.cursor.execute(sql, [sample_name]).fetchone()["id"]
        self.insert_sequence(seqhash)
        return sid

    def insert_alignment(self, seqhash: str, element_id: int) -> int:
        """
        Inserts a sequence-alignment relation into the database if not existing and returns the row id.

        Args:
            seqhash (str): The hash of the sequence to be associated with the alignment.
            element_id (int): The element id to be associated with the alignment.

        Returns:
            int: The row id of the inserted alignment.

        Raises:
            sqlite3.Error: If an error occurs while interacting with the database.

        Example usage:
            >>> dbm = getfixture('init_writeable_dbm')
            >>> rowid = dbm.insert_alignment("1a1f34ef4318911c2f98a7a1d6b7e9217c4ae1d1", 1)
        """
        # Construct the SQL query
        sql_insert = (
            "INSERT OR IGNORE INTO alignment (id, seqhash, element_id) VALUES(?, ?, ?);"
        )

        # Execute the SQL query
        self.cursor.execute(sql_insert, [None, seqhash, element_id])

        # Retrieve the alignment ID
        sql_retrieve = "SELECT id FROM alignment WHERE element_id = ? AND seqhash = ?;"
        alignment_id = self.cursor.execute(
            sql_retrieve, [element_id, seqhash]
        ).fetchone()["id"]

        return alignment_id

    def insert_molecule(
        self,
        reference_id: int,
        type: str,
        accession: str,
        symbol: str,
        description: str,
        segment: int,
        length: int,
        standard: int = 0,
    ) -> int:
        """
        Inserts a molecule into the database. If standard is 1, the molecule is set as default molecule of the respective reference. If the molecule already exists in
        the database, its properties will be updated. Returns the rowid of the inserted/updated molecule.

        Args:
            reference_id (int): The reference ID.
            type (str): The type of molecule.
            accession (str): The accession number of the molecule.
            symbol (str): The symbol representing the molecule.
            description (str): A description of the molecule.
            segment (int): The segment number.
            length (int): The length of the molecule.
            standard (int, optional): If set to 1, the molecule is set as the default. Defaults to 0.

        Returns:
            int: The rowid of the inserted/updated molecule.

        Raises:
            sqlite3.Error: If an error occurs while interacting with the database.

        Example usage:
            >>> dbm = getfixture('init_writeable_dbm')
            >>> rowid = dbm.insert_molecule(1, "plasmid", "CP028427.1", "pARLON1", "Gulosibacter molinativorax strain ON4 plasmid pARLON1, complete sequence", 1, 37013)
        """
        if not symbol.strip():
            symbol = accession
        if standard:
            self.cursor.execute(
                "UPDATE molecule SET standard = ? WHERE reference_id = ? AND standard = 1",
                [0, reference_id],
            )
        else:
            result = self.cursor.execute(
                'SELECT COUNT(*) AS count FROM referenceView WHERE "reference.id" = ? AND "molecule.standard" = 1',
                [reference_id],
            ).fetchone()
            if result["count"] == 0:
                standard = 1

        self.cursor.execute(
            "INSERT OR REPLACE INTO molecule (id, reference_id, type, accession, symbol, description, segment, length, standard) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?);",
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
        return self.cursor.execute(
            "SELECT id FROM molecule WHERE accession = ?", [accession]
        ).fetchone()["id"]

    def insert_element(
        self,
        molecule_id: int,
        type: str,
        accession: str,
        symbol: str,
        description: str,
        start: int,
        end: int,
        strand: int,
        sequence: str,
        standard: int = 0,
        parent_id: int = 0,
        parts: Optional[List[Tuple[int]]] = None,
    ) -> int:
        """
        Inserts an element (source, CDS, or protein) into the database. If the element already exists in the database,
        its properties will be updated. Returns the rowid of the inserted/updated element.

        Args:
            molecule_id (int): The molecule ID.
            type (str): The type of element.
            accession (str): The accession number of the element.
            symbol (str): The symbol representing the element.
            description (str): A description of the element.
            start (int): The starting position of the element.
            end (int): The ending position of the element.
            strand (int): The strand of the element.
            sequence (str): The sequence of the element.
            standard (int, optional): If set to 1, the element is set as the default. Defaults to 0.
            parent_id (int, optional): The ID of the parent element. Defaults to 0.
            parts (Optional[List[Tuple[int]]], optional): List of tuples containing start and end coordinates if the element
            is not linearly encoded on its molecule. Defaults to None.

        Returns:
            int: The rowid of the inserted/updated element.

        Raises:
            sqlite3.Error: If an error occurs while interacting with the database.

        Example usage:
            >>> dbm = getfixture('init_writeable_dbm')
            >>> rowid = dbm.insert_element(1, "protein", "GMOLON4_3257", "NlpD", "M23/M37 family peptidase", 5579, 6199, 1, "MKGLRSSNPKGEASD")
        """
        if not symbol.strip():
            symbol = accession
        if standard:
            self.cursor.execute(
                "UPDATE element SET standard = ? WHERE molecule_id = ? AND standard = 1",
                [0, molecule_id],
            )
        self.cursor.execute(
            "INSERT OR REPLACE INTO element (id, molecule_id, type, accession, symbol, description, start, end, strand, sequence, standard, parent_id) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);",
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
                str(sequence),
                standard,
                parent_id,
            ],
        )
        eid = self.cursor.execute(
            "SELECT id FROM element WHERE accession = ?", [accession]
        ).fetchone()["id"]
        if parts is not None:
            for part in parts:
                self.cursor.execute(
                    "INSERT OR IGNORE INTO elempart (element_id, start, end, strand, base, segment) VALUES(?, ?, ?, ?, ?, ?);",
                    [eid] + list(part),
                )
        return eid

    def insert_variant(
        self,
        alignment_id: int,
        element_id: int,
        ref: str,
        alt: str,
        start: int,
        end: int,
        label: str,
        parent_id: str = "",
        frameshift: int = 0,
    ) -> int:
        """
        Inserts a variant into the database if it does not already exist. Returns the rowid of the inserted variant.

        Args:
            alignment_id (int): The alignment ID.
            element_id (int): The element ID.
            ref (str): The reference sequence.
            alt (str): The alternate sequence.
            start (int): The starting position of the variant.
            end (int): The ending position of the variant.
            label (str): The label of the variant.
            parent_id (str, optional): The ID of the parent element. Defaults to "".
            frameshift (int, optional): If set to 1, the variant is a frameshift variant. Defaults to 0.

        Returns:
            int: The rowid of the inserted variant.

        Raises:
            sqlite3.Error: If an error occurs while interacting with the database.

        Example usage:
            >>> dbm = getfixture('init_writeable_dbm')
            >>> rowid = dbm.insert_variant(1, 1, "A", "T", 0, 1, "A1T", "", 0)
        """
        self.cursor.execute(
            "INSERT OR IGNORE INTO variant (id, element_id, start, end, ref, alt, label, parent_id, frameshift) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?);",
            [None, element_id, start, end, ref, alt, label, parent_id, frameshift],
        )
        vid = self.get_variant_id(element_id, start, end, ref, alt)
        self.cursor.execute(
            "INSERT OR IGNORE INTO alignment2variant (alignment_id, variant_id) VALUES(?, ?);",
            [alignment_id, vid],
        )
        return vid

    # DELETE DATA

    def delete_samples(self, *sample_names: str) -> None:
        """
        Deletes one or more given samples based on their names if they exist in the database.

        Args:
            sample_names (str): The names of the samples to be deleted.

        Raises:
            sqlite3.Error: If an error occurs while interacting with the database.

        Example usage:
            >>> dbm = getfixture('init_writeable_dbm')
            >>> dbm.delete_samples("NC_045512")
        """
        sample_names = list(set(sample_names))
        sql = (
            "DELETE FROM sample WHERE name IN ("
            + ", ".join(["?"] * len(sample_names))
            + ");"
        )
        self.cursor.execute(sql, sample_names)
        self.clean()

    def delete_property(self, property_name: str) -> None:
        """
        Deletes a property and all related data linked to samples based on the property name
        from the database if the property exists.

        Args:
            property_name (str): The name of the property to be deleted.

        Raises:
            sqlite3.Error: If an error occurs while interacting with the database.
        """
        if property_name in self.properties:
            sql = f"DROP INDEX idx_property_{property_name};"
            self.cursor.execute(sql)
            sql = f"ALTER TABLE property DROP COLUMN {property_name};"
            self.cursor.execute(sql)
            sql = f"DELETE FROM property_metadata WHERE name='{property_name}';"
            self.cursor.execute(sql)
            self.__properties = False
            self.clean()

    # GET DATA

    def get_sample_id(self, sample_name: str) -> Optional[int]:
        """
        Returns the rowid of a sample based on its name if it exists, else None is returned.

        Args:
            sample_name (str): Name of the sample.

        Returns:
            Optional[int]: The id of the sample if exists, else None.

        Example usage:
            >>> dbm = getfixture('init_readonly_dbm')
            >>> id = dbm.get_sample_id("seq01")
        """
        sql = "SELECT id FROM sample WHERE name = ? LIMIT 1;"
        row = self.cursor.execute(sql, [sample_name]).fetchone()
        return None if row is None else row["id"]

    def get_sample_data(self, sample_name: str) -> Tuple[Optional[int], Optional[str]]:
        """
        Returns a tuple of rowid and seqhash of a sample based on its name if it exists, else a tuple of Nones is returned.

        Args:
            sample_name (str): Name of the sample.

        Returns:
            Optional[Tuple[int,str]]: Tuple of id and seqhash of the sample, if exists, else a Tuple of None, None.

        """
        sql = "SELECT id, seqhash FROM sample WHERE name = ? LIMIT 1;"
        row = self.cursor.execute(sql, [sample_name]).fetchone()
        return (row["id"], row["seqhash"]) if row else (None, None)

    def iter_sample_names(self) -> Iterator[str]:
        """
        Iterates over all sample names stored in the database.

        Returns:
            Iterator[str]: An iterator yielding sample names from the database.

        """
        sql = "SELECT name FROM sample WHERE 1;"
        for row in self.cursor.execute(sql):
            yield row["name"]

    def get_alignment_id(self, seqhash: str, element_id: int) -> Optional[int]:
        """
        Returns the rowid of a sample based on the respective seqhash and element. If no
        alignment of the given sequence to the given element has been stored, None is returned.

        Args:
            seqhash (str): The seqhash of the sample.
            element_id (int): The element id.

        Returns:
            Optional[int]: The id of the alignment if exists, else None.

        """
        sql = "SELECT id FROM alignment WHERE element_id = ? AND seqhash = ? LIMIT 1;"
        row = self.cursor.execute(sql, [element_id, seqhash]).fetchone()
        return None if row is None else row["id"]

    def get_default_reference_accession(self) -> str:
        """
        Returns accession of the reference defined as default in the database.

        Returns:
            str: The default reference accession from the database.
        """
        sql = "SELECT accession FROM reference WHERE standard=1"
        return self.cursor.execute(sql).fetchone()["accession"]

    def get_molecule_ids(
        self, reference_accession: Optional[str] = None
    ) -> Dict[str, int]:
        """
        Returns a dictionary with accessions as keys and respective rowids as values for
        all molecules related to a given (or the default) reference.

        Args:
            reference_accession (str, optional): The reference accession. Defaults to None, which means the default reference.

        Returns:
            Dict[str, int]: Dictionary with molecule accessions as keys and their ids as values.
        """
        if reference_accession is not None:
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

    def get_molecule_data(
        self, *fields, reference_accession: Optional[str] = None
    ) -> Dict[str, Dict[str, Any]]:
        """
        Returns a dictionary with molecule accessions as keys and sub-dicts as values for all molecules
        of a given (or the default) reference. The sub-dicts store all table field names
        (or, alternatively, the given table field names only) as keys and the stored data as values.

        Args:
            *fields: Fields to be included in the sub-dicts.
            reference_accession (str, optional): The reference accession. Defaults to None, which means the default reference.

        Returns:
            Dict[str, Dict[str, Any]]: Dictionary with molecule accessions as keys and their data as values.
        """
        if not fields:
            fields = "*"
        elif "molecule.accession" not in fields:
            fields = list(set(fields)) + ['"molecule.accession"']

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
        return {x["molecule.accession"]: x for x in row} if row else {}

    def get_elements(self, molecule_id: int, *types: str) -> List[Dict[str, Any]]:
        """
        Returns a list of elements for a given molecule_id and optional types.

        Args:
            molecule_id (int): The id of the molecule.
            *types: Types of the elements.

        Returns:
            List[Dict[str, Any]]: List of elements as dictionaries.
        """
        sql = "SELECT * FROM element WHERE molecule_id = ?"
        if types:
            sql += " AND type IN (" + ", ".join(["?"] * len(types)) + ");"
        row = self.cursor.execute(sql, [molecule_id] + list(types)).fetchall()
        return row if row else []

    def get_element_ids(
        self,
        reference_accession: Optional[str] = None,
        element_type: Optional[str] = None,
    ) -> List[str]:
        """
        Returns ids of elements given a reference accession and a type.

        Args:
            reference_accession (str, optional): The accession of the reference. Defaults to None (standard reference ist used).
            element_type (str, optional): The type of the element. Defaults to None.

        Returns:
            List[str]: A list of element ids.
        """
        # Get the molecule ids given a reference accession
        molecule_ids = list(
            self.get_molecule_ids(reference_accession=reference_accession).values()
        )

        # Construct SQL query
        query = (
            'SELECT id FROM element WHERE "molecule_id" IN ('
            + ", ".join(["?"] * len(molecule_ids))
            + ")"
        )

        # Add type to the query if provided
        if element_type:
            query += " AND type = ?"
            molecule_ids.append(element_type)

        # Execute query
        rows = self.cursor.execute(query, molecule_ids).fetchall()

        # Return an empty list if no result found, otherwise return the list of ids
        return [] if not rows else [x["id"] for x in rows]

    def get_source(self, molecule_id: int) -> Optional[str]:
        """
        Returns the source data given a molecule id.

        Args:
            molecule_id (int): The id of the molecule.

        Returns:
            Optional[str]: The source if it exists, None otherwise.

        Example usage:
        >>> dbm = getfixture('init_readonly_dbm')
        >>> dbm.get_source(molecule_id=5)
        >>> dbm.get_source(molecule_id=1)['accession']
        'MN908947.3'
        """
        # Get the source elements given a molecule id
        source_elements = self.get_elements(molecule_id, "source")

        # Return None if no source elements found, otherwise return the first source element
        return None if not source_elements else source_elements[0]

    def get_annotation(
        self,
        reference_accession: Optional[str] = None,
        molecule_accession: Optional[str] = None,
        element_accession: Optional[str] = None,
        element_type: Optional[str] = None,
        fields: List[str] = ["*"],
    ) -> List[Dict[str, Any]]:
        """
        Retrieves the annotation based on the given accessions and a type.

        Args:
            reference_accession (str, optional): Accession of the reference. Defaults to None.
            molecule_accession (str, optional): Accession of the molecule. Defaults to None.
            element_accession (str, optional): Accession of the element. Defaults to None.
            element_type (str, optional): Type of the element. Defaults to None.
            fields (List[str], optional): Fields to be selected in the query. Defaults to ["*"].

        Returns:
            List[Dict[str, Any]]: A list of results, each result is a dictionary containing field names as keys and the corresponding data as values.
        """
        # Prepare the conditions and values for the query
        conditions = []
        values = []

        # Assemble conditions defining accessions or standard elements
        if reference_accession:
            conditions.append('"reference.accession" = ?')
            values.append(reference_accession)
        else:
            conditions.append('"reference.standard" = ?')
            values.append(1)

        if molecule_accession:
            conditions.append('"molecule.accession" = ?')
            values.append(molecule_accession)
        else:
            conditions.append(
                '("molecule.standard" IS NULL OR "molecule.standard" = 1)'
            )

        if element_accession:
            conditions.append('"element.accession" = ?')
            values.append(element_accession)
        elif not element_type:
            conditions.append('("element.type" IS NULL OR "element.type" = \'source\')')

        if element_type:
            conditions.append('"element.type" = ?')
            values.append(element_type)

        # Construct SQL query to fetch the annotation
        sql = (
            "SELECT "
            + ", ".join(fields)
            + " FROM referenceView WHERE "
            + " AND ".join(conditions)
            + ' ORDER BY "reference.id" ASC, "molecule.id" ASC, "element.id" ASC, "element.segment" ASC'
        )

        # Execute the query and return the results
        return self.cursor.execute(sql, values).fetchall()

    def get_alignment_data(
        self, sample_name: str, reference_accession: Optional[str] = None
    ) -> Union[sqlite3.Cursor, str]:
        """
        Retrieves the alignment data for the given sample and reference accession.

        Args:
            sample_name (str): Name of the sample.
            reference_accession (str, optional): Accession of the reference. Defaults to None.

        Returns:
            sqlite3.Cursor or str: Cursor object with query results or an empty string if no reference accession is found.
        """
        # If reference_accession is not provided, fetch it from the database
        if reference_accession is None:
            sql = 'SELECT "reference.accession" as acc FROM alignmentView WHERE "sample.name" = ? LIMIT 1'
            reference_accession = self.cursor.execute(sql, [sample_name]).fetchone()
            if reference_accession is None:
                return ""
            reference_accession = reference_accession["acc"]

        # Construct SQL query to fetch the alignment data
        sql = 'SELECT "element.sequence", "element.symbol", "element.id" FROM alignmentView WHERE "sample.name" = ? AND "reference.accession" = ?'

        # Execute the query and return the results
        return self.cursor.execute(sql, [sample_name, reference_accession])

    def get_variant_id(
        self, element_id: int, start: int, end: int, ref: str, alt: str
    ) -> Optional[int]:
        """
        Retrieves the ID of the variant based on the given parameters.

        Args:
            element_id (int): ID of the element.
            start (int): Start position of the variant.
            end (int): End position of the variant.
            ref (str): Reference base(s).
            alt (str): Alternative base(s).

        Returns:
            int or None: ID of the variant if found, else None.
        """
        # Construct SQL query to fetch the variant ID
        sql = "SELECT id FROM variant WHERE element_id = ? AND start = ? AND end = ? AND ref = ? AND alt = ?;"

        # Execute the query and fetch the result
        row = self.cursor.execute(sql, [element_id, start, end, ref, alt]).fetchone()

        # Return the variant ID if found, else None
        return None if row is None else row["id"]

    def iter_dna_variants(
        self, sample_name: str, *element_ids: int
    ) -> Iterator[sqlite3.Row]:
        """
        Iterator over DNA variants for a given sample and a list of element IDs.

        Args:
            sample_name (str): Name of the sample to be analyzed.
            *element_ids (int): Variable length argument list of element IDs.

        Returns:
            Iterator[sqlite3.Row]: An iterator that yields rows from the executed SQL query.

        Yields:
            sqlite3.Row: The next row from the executed SQL query.
        """

        if element_ids:
            condition = f" IN ({', '.join(['?']*len(element_ids))})"
        else:
            condition = ""

        sql = f"""SELECT  variant.element_id as \"element.id\",
                    variant.start as \"variant.start\",
                    variant.end as \"variant.end\",
                    variant.ref as \"variant.ref\",
                    variant.alt as \"variant.alt\"
                FROM
                    (SELECT sample.seqhash
                    FROM sample
                    WHERE sample.name = ?) AS sample_T
                INNER JOIN alignment
                    ON sample_T.seqhash == alignment.seqhash
                INNER JOIN alignment2variant
                    ON alignment.id == alignment2variant.alignment_id
                INNER JOIN variant
                    ON alignment2variant.variant_id == variant.id
                WHERE  variant.element_id{condition}"""

        for row in self.cursor.execute(sql, [sample_name] + list(element_ids)):
            if row["variant.start"] is not None:
                yield row

    def get_translation_dict(self, translation_id: int) -> Dict[str, str]:
        """
        Returns a dictionary of codon to amino acid mappings for a given translation ID.

        Args:
            translation_id (int): The ID of the translation table.

        Returns:
            Dict[str, str]: A dictionary where each key-value pair represents a codon and its corresponding amino acid.
        """
        sql = "SELECT codon, aa FROM translation WHERE id = ?"
        return {
            row["codon"]: row["aa"]
            for row in self.cursor.execute(sql, [translation_id]).fetchall()
        }

    def get_element_parts(self, element_id: int) -> List[Dict[str, Any]]:
        """
        Returns 'start', 'end', and 'strand' information of segments of a given element.

        Args:
            element_id (int): The ID of the element.

        Returns:
            List[Dict[str, Any]]: A list of dictionaries where each dictionary contains the 'start', 'end', and 'strand' of a part of the element.
        """
        sql = "SELECT start, end, strand FROM elempart WHERE element_id = ? ORDER BY segment"
        return self.cursor.execute(sql, [element_id]).fetchall()

    def get_sequence(self, element_id: int) -> Optional[str]:
        """
        Returns the sequence of a given element.

        Args:
            element_id (int): The ID of the element.

        Returns:
            Optional[str]: The sequence of the element if it exists, else None.
        """
        sql = "SELECT sequence FROM element WHERE id = ?"
        row = self.cursor.execute(sql, [element_id]).fetchone()
        return None if row is None else row["sequence"]

    def extract_sequence(
        self,
        element_id: int = None,
        translation_table: Optional[int] = None,
        molecule_id: int = None,
    ) -> Optional[str]:
        """
        Extracts the sequence of an element except CDS or SOURCE. If a translation table is provided, the sequence is translated into protein sequence.

        Args:
            element_id (int, optional): The ID of the element. Defaults to None.
            translation_table (int, optional): The ID of the translation table. Defaults to None.
            molecule_id (int, optional): The ID of the molecule. Defaults to None.

        Returns:
            Optional[str]: The extracted (and possibly translated) sequence if the element exists, else None.
        """
        # fetch relevant data from database
        sql = "SELECT sequence, type, id, parent_id FROM element WHERE id = ? AND molecule_id = ?;"
        row = self.cursor.execute(sql, [element_id, molecule_id]).fetchone()
        if not row:
            return None
        element_id = row["id"]

        # select most parent element (but not (source)
        while row and row["type"] not in {"source", "CDS"}:
            sql = "SELECT sequence, type, parent_id FROM element WHERE id = ?"
            row = self.cursor.execute(sql, [row["parent_id"]]).fetchone()
        sequence = row["sequence"]

        # assemble element segments (e.g. introns)
        parts = [
            FeatureLocation(part["start"], part["end"], strand=part["strand"])
            for part in self.get_element_parts(element_id)
        ]
        feat = CompoundLocation(parts) if len(parts) > 1 else parts[0]
        if translation_table is None:
            return str(feat.extract(sequence))
        return str(
            Seq(feat.extract(sequence)).translate(
                table=translation_table, stop_symbol=""
            )
        )

    # COUNTING DATA

    def count_samples(self) -> int:
        """
        Count the number of samples in the database.

        Returns:
            int: The total number of samples.
        """
        sql = "SELECT COUNT(*) FROM sample;"
        return self.cursor.execute(sql).fetchone()["COUNT(*)"]

    def count_sequences(self) -> int:
        """
        Count the number of distinct sequences in the samples.

        Returns:
            int: The total number of distinct sequences.
        """
        sql = "SELECT COUNT(DISTINCT seqhash) FROM sample;"
        return self.cursor.execute(sql).fetchone()["COUNT(DISTINCT seqhash)"]

    def count_property(
        self, property_name: str, distinct: bool = False, ignore_standard: bool = False
    ) -> int:
        """
        Count the number of properties.

        Args:
            property_name (str): The name of the property.
            distinct (bool, optional): Whether to count distinct values. Defaults to False.
            ignore_standard (bool, optional): Whether to ignore standard values. Defaults to False.

        Returns:
            int: The total count of the specified property.
        """
        # Form the SQL query based on the given arguments
        distinct_flag = "DISTINCT " if distinct else ""
        conditions = f"WHERE {property_name} != ?" if ignore_standard else ""
        vals = [self.properties[property_name]["standard"]] if ignore_standard else []

        sql = (
            f"SELECT COUNT({distinct_flag} {property_name}) "
            f"as count FROM property {conditions};"
        )

        return self.cursor.execute(sql, vals).fetchone()["count"]

    def count_variants(self, protein_level: bool = False) -> int:
        """
        Count the number of nucleotide-level mutations stored in the database.

        Args:
            protein_level (bool): If true, coutn protein-level mutations intead of nucleotide level mutations.

        Returns:
            int: The total number of nucleotide-level mutations.
        """
        element_type = "source" if not protein_level else "cds"
        sql = f"""
                SELECT COUNT(DISTINCT "element.id" || '_' || "variant.start" || '_' || "variant.end") AS count
                    FROM variantView
                WHERE "element.type" = '{element_type}';
            """
        return self.cursor.execute(sql).fetchone()["count"]

    # CONDITIONING SAMPLE PROPERTIES
    @staticmethod
    def get_operator(op: Optional[str] = None, inverse: bool = False) -> str:
        """Returns the appropriate operator for a given operation type and inversion flag.

        Args:
            op (Optional[str]): The operation to be performed. If not specified or empty, the default operation '=' will be used.
            inverse (bool, optional): A flag indicating whether to use the inverse of the operation. Defaults to False.

        Returns:
            str: The operator corresponding to the operation type and inversion flag.
        """
        op = op if op else sonarDbManager.OPERATORS["default"]
        return (
            sonarDbManager.OPERATORS["inverse"][op]
            if inverse
            else sonarDbManager.OPERATORS["standard"][op]
        )

    def get_conditional_expr(
        self, field: str, operator: str, *vals: Union[str, Tuple[str, str]]
    ) -> Tuple[List[str], List[str]]:
        """
        Generates SQL conditional expressions and their corresponding values for given field, operator, and values.

        Args:
            field (str): Name of the SQL field.
            operator (str): SQL operator to be used in the conditional expressions.
            vals (Union[str, Tuple[str, str]]): Values to be used in the conditional expressions. Can be either a single
            value or a tuple for range-based operations.

        Returns:
            Tuple[List[str], List[str]]: Tuple of lists of conditional expressions and corresponding values.
        """
        conditions = []
        values = []

        # transforming operators
        if len(vals) > 1:
            if operator in ["=", "!="]:
                operator = "IN" if operator == "=" else "NOT IN"
        elif operator in ["IN", "NOT IN"]:
            operator = "=" if operator == "IN" else "!="

        # creating conditions
        if operator in ["IN", "NOT IN"]:
            conditions.append(f"{field} {operator} ({', '.join(['?' for _ in vals ])})")
            values.extend(vals)

        elif operator in ["BETWEEN", "NOT BETWEEN"]:
            conditions.extend(
                [f"{field} {operator} ? AND ?" for _ in range(0, len(vals), 2)]
            )
            values.extend(vals)

        else:
            conditions.extend([f"{field} {operator} ?" for _ in vals])
            values.extend(vals)

        return conditions, values

    def build_numeric_condition(
        self, field: str, *vals: str, logic_link: str = "AND"
    ) -> Tuple[str, List[Union[str, Tuple[str, str]]]]:
        """
        Construct conditions for numeric fields.

        Args:
            field (str): The name of the field to query.
            vals (str): The values to search for in the field.
            logic_link (str, optional): The logica link for subqueries. Defaults to "AND".

        Returns:
            Tuple[str, List[Union[str, Tuple[str, str]]]]: Tuple containing the formatted query and the list of values.

        Raises:
            SystemExit: Raises error if the provided values don't match the expected format.
        """

        logic_link = f" {logic_link.strip()} "
        pattern_single = re.compile(r"^(\^*)((?:>|>=|<|<=|!=|=)?)(-?[1-9]+[0-9]*)$")
        pattern_range = re.compile(r"^(\^*)(-?[1-9]+[0-9]*):(-?[1-9]+[0-9]*)$")
        error_msg = (
            f"query error: numeric value or range expected for field {field}(got: "
        )
        data = defaultdict(list)

        for val in vals:
            val = str(val).strip()

            # Processing single value
            if ":" not in val:
                match = pattern_single.match(val)
                if not match:
                    LOGGER.error(f"{error_msg}{val})")
                    sys.exit(1)
                operator = sonarDbManager.get_operator(match.group(2), match.group(1))
                num = int(match.group(3))

                data[operator].append(num)

            # Processing value range
            else:
                match = pattern_range.match(val)
                if not match:
                    LOGGER.error(f"{error_msg}{val})")
                    sys.exit(1)
                operator = sonarDbManager.get_operator(
                    sonarDbManager.OPERATORS["standard"]["BETWEEN"], match.group(1)
                )

                num1 = int(match.group(2))
                num2 = int(match.group(3))

                # Plausibility check
                if num1 >= num2:
                    LOGGER.error(f"Invalid range ({match.group(0)}).")
                    sys.exit(1)

                data[operator] += [num1, num2]

        # Assemble query conditions and values
        conditions = []
        values = []
        for operator, vals in data.items():
            c, v = self.get_conditional_expr(field, operator, *vals)
            conditions.extend(c)
            values.extend(v)

        # Return query conditions and values as strings
        return logic_link.join(conditions), values

    def build_float_condition(
        self, field: str, *vals: Union[float, str], logic_link: str = "AND"
    ) -> Tuple[str, List[float]]:
        """
        Builds a float condition for a given field and values.

        Args:
            field (str): Field for the condition.
            vals (Union[float, str]): Values to be used in the condition. Each value can be a float or a string representing a single value or a range.
            logic_link (str, optional): Logical operator to link the conditions. Defaults to 'AND'.

        Returns:
            Tuple[str, List[float]]: Returns the constructed condition and a list of corresponding values.
        """

        logic_link = f" {logic_link.strip()} "
        single_value_pattern = re.compile(
            r"^(\^*)((?:>|>=|<|<=|!=|=)?)(-?[1-9]+[0-9]*(?:.[0-9]+)*)$"
        )
        range_value_pattern = re.compile(
            r"^(\^*)(-?[1-9]+[0-9]*(?:.[0-9]+)*):(-?[1-9]+[0-9]*(?:.[0-9]+)*)$"
        )

        data = defaultdict(list)
        for val in vals:
            val = str(val).strip()
            if ":" not in val:
                # Processing single value
                match = single_value_pattern.match(val)
                if not match:
                    LOGGER.error(
                        f"Decimal value or range expected for field {field} (got: {val})"
                    )
                    sys.exit(1)

                operator = sonarDbManager.get_operator(match.group(2), match.group(1))
                decinum = float(match.group(3))

                data[operator].append(decinum)
            else:
                # Processing range
                match = range_value_pattern.match(val)
                if not match:
                    LOGGER.error(
                        f"Decimal value or range expected for field {field} (got: {val})"
                    )
                    sys.exit(1)

                operator = sonarDbManager.get_operator(
                    sonarDbManager.OPERATORS["standard"]["BETWEEN"], match.group(1)
                )
                decinum1 = float(match.group(2))
                decinum2 = float(match.group(3))

                # Plausibility check
                if decinum1 >= decinum2:
                    LOGGER.error(f"Invalid range ({match.group(0)}).")
                    sys.exit(1)

                data[operator] += [decinum1, decinum2]

        # Assemble query conditions and values
        conditions = []
        values = []
        for operator, vals in data.items():
            c, v = self.get_conditional_expr(field, operator, *vals)
            conditions.extend(c)
            values.extend(v)

        # Return query conditions and values as strings
        return logic_link.join(conditions), values

    def build_date_condition(
        self, field: str, *vals, logic_link: str = "AND"
    ) -> Tuple[str, List[str]]:
        """
        Construct conditions for date fields.

        Args:
            field (str): The field to query.
            *vals: The values to query.
            link (str): The link operator for the query conditions. Default is "AND".

        Returns:
            Tuple[str, List[str]]: The conditions for the query and the values to be used in the query.
        """
        logic_link = f" {logic_link.strip()} "
        pattern_single = re.compile(
            r"^(\^*)((?:>|>=|<|<=|!=|=)?)([0-9]{4}-[0-9]{2}-[0-9]{2})$"
        )
        pattern_range = re.compile(
            r"^(\^*)([0-9]{4}-[0-9]{2}-[0-9]{2}):([0-9]{4}-[0-9]{2}-[0-9]{2})$"
        )
        error_msg = (
            "query error: date or date range expected for field " + field + " (got: "
        )
        data = defaultdict(list)

        for val in vals:
            # processing single value
            if ":" not in val:
                match = pattern_single.match(val)
                if not match:
                    LOGGER.error(f"{error_msg}{val})")
                    sys.exit(1)
                operator = sonarDbManager.get_operator(match.group(2), match.group(1))
                data[operator].append(match.group(3))
            # processing value range
            else:
                match = pattern_range.match(val)
                if not match:
                    LOGGER.error(f"{error_msg}{val})")
                    sys.exit(1)
                operator = sonarDbManager.get_operator(
                    sonarDbManager.OPERATORS["standard"]["BETWEEN"], match.group(1)
                )

                date1 = datetime.strptime(match.group(2), "%Y-%m-%d")
                date2 = datetime.strptime(match.group(3), "%Y-%m-%d")

                # Plausibility check
                if date1 >= date2:
                    LOGGER.error("Invalid range (" + match.group(0) + ").")
                    sys.exit(1)

                data[operator] += [match.group(2), match.group(3)]

        # Assemble query conditions and values
        conditions = []
        values = []
        for operator, vals in data.items():
            c, v = self.get_conditional_expr(field, operator, *vals)
            conditions.extend(c)
            values.extend(v)

        # Return query conditions and values as strings
        return logic_link.join(conditions), values

    @staticmethod
    def custom_strip(s, char_to_remove):
        if s.startswith(char_to_remove):
            s = s[1:]
        if s.endswith(char_to_remove):
            s = s[:-1]
        return s

    def build_string_condition(
        self, field: str, *vals: str, logic_link: str = "AND"
    ) -> Tuple[str, List[Any]]:
        """
        Construct conditions for string fields.

        Args:
            field (str): The field to query.
            vals (str): The values to query for.
            link (str, optional): The logic link to use between conditions. Defaults to "AND".

        Returns:
            Tuple[str, List[Any]]: The assembled query conditions as a string and the corresponding values as a list.
        """
        logic_link = f" {logic_link.strip()} "
        data = defaultdict(list)

        for val in vals:
            # Determine the operation key and strip the inverse symbol if necessary
            negate = True if val.startswith("^") else False
            val = val[1:] if negate else val

            # Determine the operator
            operator = self.get_operator(
                "LIKE" if val.startswith("%") or val.endswith("%") else "=", negate
            )

            # Add the value to the corresponding operator set
            data[operator].append(val)

        # Assemble query conditions and values
        conditions = []
        values = []
        for operator, vals in data.items():
            c, v = self.get_conditional_expr(field, operator, *vals)
            conditions.extend(c)
            values.extend(v)

        # Return query conditions and values as strings
        return logic_link.join(conditions), values

    def build_zip_condition(
        self, field: str, *vals: str, logic_link: str = "AND"
    ) -> Tuple[str, List[Any]]:
        """
        Construct conditions for zip code fields.

        Args:
            field (str): The field to query.
            vals (str): The values to query for. '%' can be used for wildcard matches.
            logic_link (str, optional): The logic link to use between conditions. Defaults to "AND".

        Returns:
            Tuple[str, List[Any]]: The assembled query conditions as a string and the corresponding values as a list.
        """
        logic_link = f" {logic_link.strip()} "
        data = defaultdict(list)

        for val in vals:
            # Determine the operation key and strip the inverse symbol if necessary
            negate = True if val.startswith("^") else False
            val = val[1:] if negate else val

            # Add the '%' symbol at the end of the zip code for LIKE queries
            orig_val = val
            val = sonarDbManager.custom_strip(val, "%")

            # Plausibility check
            try:
                _ = int(val)
            except ValueError:
                LOGGER.error(f"Invalid range ({orig_val}).")
                sys.exit(1)

            # Get the operator for LIKE queries
            operator = self.get_operator("LIKE", negate)

            # Add the value to the corresponding operator set
            data[operator].append(val)

        # Assemble query conditions and values
        conditions = []
        values = []
        for operator, vals in data.items():
            c, v = self.get_conditional_expr(field, operator, *vals)
            conditions.extend(c)
            values.extend(v)

        # Return query conditions and values as strings
        return logic_link.join(conditions), values

    def resolve_pango_sublineages(self, *lineages: str) -> Set[str]:
        """
        Resolve PANGO sublineages from given lineages.

        Args:
            lineages (str): Lineages to which sublineages will be added.

        Returns:
            Set[str]: A set of all lineages including added sublineages.
        """
        lineage_pool = set()
        lineages = list(set(lineages))

        while lineages:
            lineage = lineages.pop()

            negate = True if lineage.startswith("^") else False
            op = "^" if negate else ""
            if negate:
                lineage = lineage[1:]

            # If the lineage ends with "*", retrieve sublineages
            if lineage.endswith("*"):
                lineage = lineage[:-1]
                sublineages = self.lineage_sublineage_dict.get(lineage, "").split(",")
                lineages.extend((op + sublineage for sublineage in sublineages))

            lineage_pool.add(op + lineage)

        return lineage_pool

    def resolve_pango_wildcards(self, *lineages: str) -> Set[str]:
        """
        Resolve wildcard characters in the given lineages.

        Args:
            lineages (str): Lineages containing wildcard characters.

        Returns:
            Set[str]: A set of all lineages with wildcard characters resolved.
        """
        lineage_pool = set()

        for lineage in lineages:
            negate = True if lineage.startswith("^") else False
            op = "^" if negate else ""
            if negate:
                lineage = lineage[1:]

            if "%" in lineage:
                lineage = lineage.replace("~", "~~").replace("_", "~_")

                if lineage.endswith("*"):
                    suff = "*"
                    lineage = lineage[:-1]
                else:
                    suff = ""

                # Prepare the SQL query
                sql = "SELECT DISTINCT lineage FROM lineages WHERE lineage LIKE ? ESCAPE '~';"

                # Execute the query and update the lineage pool
                lineage_pool.update(
                    (
                        op + res["lineage"] + suff
                        for res in self.cursor.execute(sql, [lineage]).fetchall()
                    )
                )
            else:
                lineage_pool.add(lineage)

        return lineage_pool

    def build_pango_condition(
        self, field: str, *vals: str, logic_link: str = "AND"
    ) -> Tuple[str, List[str]]:
        """
        Construct conditions of PANGO fields.

        Args:
            field (str): The field to query.
            vals (str): The values to query for.
            logic_link (str, optional): The logic link to use between conditions. Defaults to "AND".

        Returns:
            Tuple[str, List[str]]: The assembled query conditions as a string and the corresponding values as a list.
        """
        logic_link = f" {logic_link.strip()} "
        data = defaultdict(list)

        # Resolve PANGO wildcards and add sublineages
        vals = self.resolve_pango_wildcards(*vals)
        vals = self.resolve_pango_sublineages(*vals)

        for val in vals:
            # If the lineage starts with "^", add it to the "!=" set, else add it to the "=" set
            data["!=" if val.startswith("^") else "="].append(
                val[1:] if val.startswith("^") else val
            )

        # Assemble query conditions and values
        conditions = []
        values = []
        for operator, vals in data.items():
            c, v = self.get_conditional_expr(field, operator, *vals)
            conditions.extend(c)
            values.extend(v)

        # Return query conditions and values as strings
        return logic_link.join(conditions), values

    # CONDITIONING VARIANTS

    def build_generic_variant_condition(
        self, match: re.Match
    ) -> Tuple[List[str], List[str], Dict[str, Set[str]]]:
        """
        Extracts variant information based on the given regex match object and generates conditions for SQL query.

        Args:
            match (re.Match): A regex match object containing parsed variant details.

        Returns:
            Tuple[List[str], List[str], Dict[str, Set[str]]]: A tuple containing three elements:
                - A list of conditions for the SQL query.
                - A list of corresponding values to be used in the SQL query.
                - A dictionary representing the IUPAC code set that corresponds to the variant type (either nucleotide or amino acid).
        """
        conditions = []
        values = []

        # Set molecule symbol or standard
        conditions.append(
            "molecule.symbol = ?" if match.group(1) else "molecule.standard = ?"
        )
        values.append(match.group(1)[:-1] if match.group(1) else 1)

        # Set element type and symbol or standard
        if match.group(2):
            conditions.extend(["element.type = ?", "element.symbol = ?"])
            values.extend(["cds", match.group(2)[:-1]])
            iupac_code = self.IUPAC_CODES["aa"]
        else:
            conditions.append("element.standard = ?")
            values.append(1)
            iupac_code = self.IUPAC_CODES["nt"]

        return conditions, values, iupac_code

    def build_snp_and_insert_condition(
        self,
        match: re.Match,
    ) -> Tuple[List[str], List[Any]]:
        """
        Helper method to process SNP and insertions.

        Args:
            match (re.Match): Regex match object.

        Returns:
            Tuple[List[str], List[Any]]: Tuple of conditions and corresponding values.
        """
        conditions, values, iupac_code = self.build_generic_variant_condition(match)

        # process general variant information
        conditions.extend(["variant.start = ?", "variant.end = ?", "variant.ref = ?"])
        values.extend([int(match.group(4)) - 1, int(match.group(4)), match.group(3)])

        # handling different forms of alternate allele
        alt = match.group(5)
        if alt.startswith("="):
            conditions.append("variant.alt = ?")
            values.append(alt[1:])
        else:
            try:
                if len(alt) == 1:
                    resolved_alt = iupac_code[alt]
                else:
                    resolved_alt = [
                        "".join(x)
                        for x in itertools.product(*[iupac_code[x] for x in alt])
                    ]

                if len(resolved_alt) == 1:
                    conditions.append("variant.alt = ?")
                else:
                    conditions.append(
                        f"variant.alt IN ({', '.join('?' for _ in range(len(resolved_alt)))})"
                    )

                values.extend(resolved_alt)
            except KeyError:
                LOGGER.error(f"Invalid alternate allele notation '{alt}'.")
                sys.exit(1)

        return conditions, values

    def build_deletion_condition(
        self, match: re.Match
    ) -> Tuple[List[str], List[Union[str, int]]]:
        """
        Helper method to process deletions.

        Args:
            match (re.Match): Regex match object.
        Returns:
            Tuple[List[str], List[Union[str, int]]]: Updated lists of conditions and values.
        """
        conditions, values, _ = self.build_generic_variant_condition(match)

        start, end = match.group(3), match.group(4)[1:]

        # set deletion start
        if start.startswith("="):
            conditions.append("variant.start = ?")
            values.append(int(start[1:]) - 1)
        else:
            conditions.append("variant.start <= ?")
            values.append(int(start) - 1)

        # set deletion end
        if end.startswith("="):
            conditions.append("variant.end = ?")
            end = int(end[1:])
        else:
            conditions.append("variant.end >= ?")
        values.append(int(end))

        conditions.append("variant.alt = ?")
        values.append(" ")

        return conditions, values

    # MATCHING QUERIES

    def build_sample_property_condition(
        self, name: str, *vals: Union[str, Tuple[str, str]]
    ) -> Tuple[List[str], List[Union[str, int, float]]]:
        """
        Creates SQL WHERE clause and corresponding list of values to define sample properties based on property name and values.

        Args:
            name (str): The property name to query.
            *vals (Union[str, Tuple[str, str]]): The values to apply the query.

        Returns:
            Tuple[List[str], List[Union[str, int, float]]]: A tuple containing a list of conditional clauses and a list of corresponding values.
        """
        if name not in self.properties:
            LOGGER.error(f"Property '{name}' is unkown.")
            sys.exit(1)

        conditions = []
        values = []
        data_field = name
        query_type = self.properties[name]["querytype"]

        # map between the query type and the corresponding method
        query_functions = {
            "date": self.build_date_condition,
            "numeric": self.build_numeric_condition,
            "text": self.build_string_condition,
            "zip": self.build_zip_condition,
            "float": self.build_float_condition,
            "pango": self.build_pango_condition,
        }

        query_function = query_functions.get(query_type)
        if query_function is None:
            LOGGER.error(f"Unknown query type '{query_type}' for property '{name}'.")
            sys.exit(1)

        condition, vals = query_function(data_field, *vals)
        conditions.append(condition)
        values.extend(vals)

        return conditions, values

    def create_sample_property_case(
        self,
        properties: Optional[
            Dict[str, Tuple[str, List[Union[str, int, float]]]]
        ] = None,
    ) -> Tuple[List[str], List[str], List[Union[str, int, float]]]:
        """
        Create SQL WHERE clause for sample metadata-based filtering.

        Args:
            properties: A dictionary mapping property names to a list of their values.

        Returns:
            A tuple where:
                - first element is a list of CASE statements,
                - second element is a list of WHERE conditions,
                - third element is a list of values associated with the queries.
        """
        property_cases = []
        property_conditions = []
        property_vals = []

        if not properties:
            return property_conditions, property_cases, property_vals

        pid = 0
        for pname, vals in properties.items():
            if not vals:
                continue
            pid += 1
            case, val = self.build_sample_property_condition(pname.lstrip("."), *vals)
            property_vals.extend(val)
            property_cases.append(
                f"SUM(CASE WHEN {' AND '.join(case)} THEN 1 ELSE 0 END) AS property_{pid}"
            )
            property_conditions.append(f"property_{pid} >= 1")

        return property_cases, property_conditions, property_vals

    def create_profile_cases(
        self, *profiles: Tuple[str, ...]
    ) -> Tuple[List[str], List[str], List[Union[str, int, float]]]:
        """
        Create SQL CASE and WHERE clauses for genomic profile-based filtering.

        Args:
            profiles: A list of tuples, where each tuple contains variant notations representing a genomic profile.

        Returns:
            A tuple where:
                - first element is a list of CASE statements,
                - second element is a list of WHERE conditions,
                - third element is a list of values associated with the query.
        """
        regexes = {
            "snv": re.compile(r"^(|[^:]+:)?([^:]+:)?([A-Z]+)([0-9]+)(=?[A-Zxn]+)$"),
            "del": re.compile(r"^(|[^:]+:)?([^:]+:)?del:(=?[0-9]+)(|-=?[0-9]+)?$"),
        }

        processing_funcs = {
            "snv": self.build_snp_and_insert_condition,
            "del": self.build_deletion_condition,
        }

        ids = {}
        cases = []
        wheres = []
        vals = []

        for profile in profiles:
            where_conditions = []

            for mutation in profile:
                count = 1 if not mutation.startswith("^") else 0
                mutation = mutation.lstrip("^")

                # create case if novel mutation
                if mutation not in ids:
                    ids[mutation] = len(ids) + 1
                    for mutation_type, regex in regexes.items():
                        match = regex.match(mutation)
                        if match:
                            case, val = processing_funcs[mutation_type](match)
                            cases.append(
                                f"SUM(CASE WHEN {' AND '.join(case)} THEN 1 ELSE 0 END) AS mutation_{ids[mutation]}"
                            )
                            vals.extend(val)
                            break
                    if not match:
                        LOGGER.error(f"Invalid mutation notation '{mutation}'.")
                        sys.exit(1)

                if count == 0:
                    where_conditions.append(f"mutation_{ids[mutation]} = {count}")
                else:
                    where_conditions.append(f"mutation_{ids[mutation]} >= {count}")

                if len(where_conditions) == 1:
                    wheres.extend(where_conditions)
                elif len(where_conditions) > 1:
                    wheres.append("(" + " AND ".join(where_conditions) + ")")

        return cases, wheres, vals

    def create_sample_selection_sql(
        self,
        samples: Optional[List[str]] = None,
        properties: Optional[Dict[str, List[str]]] = None,
        profiles: Optional[Dict[str, List[str]]] = None,
        frameshifts_only: bool = None,
    ) -> Tuple[str, List[str]]:
        """
        Create a SQL query and corrspond value list to rertieve sample IDs based on the given sample names, properties, and genomic profiles.

        Args:
            samples (Optional[List[str]]): A list of samples to consider for query creation. Default is None.
            properties (Optional[Dict[str, List[str]]]): A dictionary of properties for query creation. Default is None.
            profiles (Optional[Dict[str, List[str]]]): A dictionary of profiles for query creation. Default is None.
            frameshifts_only (bool): If true, consider samples with frameshift mutations only.

        Returns:
            Tuple[str, List[str]]: A SQL query to retrieve matching sample IDs and a list of property and profile values.
        """

        # breakpoint()
        conditions = []
        vals = []

        # set framsehift-related table data
        if frameshifts_only:
            table = """(
                        SELECT DISTINCT sample.id AS id, sample.name AS name, sample.seqhash AS seqhash
                        FROM sample
                        JOIN alignment ON sample.seqhash = alignment.seqhash
                        JOIN alignment2variant ON alignment.id = alignment2variant.alignment_id
                        JOIN variant ON alignment2variant.variant_id = variant.id
                        WHERE variant.frameshift = 1
                    )"""
        else:
            table = "sample"

        # add sample-related condition
        if len(samples) == 1:
            conditions.append("s.name = ?")
            vals.extend(samples)
        elif len(samples) > 1:
            sample_wildcards = ", ".join(["?"] * len(samples))
            conditions.append(f"s.name IN ({sample_wildcards})")
            vals.extend(samples)

        # add properties- and profile-related conditions
        cases = []

        property_cases, property_conditions, property_vals = (
            self.create_sample_property_case(properties) if properties else ([], [], [])
        )

        profile_cases, profile_conditions, profile_vals = (
            self.create_profile_cases(*profiles) if profiles else ([], [], [])
        )

        sql = "SELECT sub.sample_id, sub.name, sub.seqhash FROM (SELECT s.id AS 'sample_id', s.name , s.seqhash"
        joins = ""

        # add joins and cases for sample propetries
        if property_cases:
            joins += """JOIN property ON s.id = property.sample_id\n
                     """
            cases.extend(property_cases)
            conditions.extend(property_conditions)
            vals.extend(property_vals)

        # add joins and cases for genome profiles
        if profile_cases:
            joins += """JOIN alignment ON s.seqhash = alignment.seqhash
                      JOIN alignment2variant ON alignment.id = alignment2variant.alignment_id
                      JOIN variant ON alignment2variant.variant_id = variant.id
                      JOIN element ON variant.element_id = element.id
                      JOIN molecule ON element.molecule_id = molecule.id\n"""

            cases.extend(profile_cases)
            if len(profile_conditions) == 1:
                conditions.extend(profile_conditions)
            else:
                conditions.append("(" + " OR ".join(profile_conditions) + ")")
            vals.extend(profile_vals)

        if cases:
            sql += ", " + ", ".join(cases)
            sql += f" FROM {table} s"
            sql += " " + joins + " "

        else:
            sql += f" FROM {table} s"

        if conditions:
            conditions = " AND ".join(conditions)
            sql += f" GROUP BY s.id HAVING {conditions}"

        sql += ") AS sub"

        # sql += " GROUP BY sub.sample_id"

        sql = self.format_sql(sql)

        return sql, vals

    def create_genomic_element_conditions(
        self,
        reference_accession: str,
        element_alias: Optional[str] = "element",
        molecule_alias: Optional[str] = "molecule",
    ) -> Tuple[str, str]:
        """
        Generate the genomic element conditions for the SQL query.

        Args:
            reference_accession (str): The reference accession to create genomic element conditions.
            element_alias: Alias used for element table
            molecule_alias: Alias used for molecule table

        Returns:
            Tuple[str, str]: A tuple containing the genomic element condition and the molecule prefix.
        """
        element_ids = self.get_element_ids(reference_accession, "source")
        if len(element_ids) == 1:
            genome_element_condition = f"{element_alias}.id = {element_ids[0]}"
            molecule_prefix = ""
        else:
            formatted_ids = ", ".join(map(str, element_ids))
            genome_element_condition = f"{element_alias}.id IN ({formatted_ids})"
            molecule_prefix = f'{molecule_alias}.symbol || "@" || '

        return genome_element_condition, molecule_prefix

    def create_special_variant_filter_conditions(
        self, filter_n: bool, filter_x: bool, ignore_terminal_gaps: bool
    ) -> Tuple[str, str]:
        """
        Create query conditions for special variants such as ambiguities and terminal gaps.

        Args:
            filter_n (bool): Flag for handling 'N' variants.
            filter_x (bool): Flag for handling 'X' variants.
            ignore_terminal_gaps (bool): Flag for handling terminal gap variants.

        Returns:
            Tuple[str, str]: A tuple containing the nucleotide filter and amino acid filter conditions.
        """

        # process mutation display filters
        aa_filter = " AND alt != 'X'" if filter_x else ""
        nuc_filter = []
        if filter_n:
            nuc_filter.append("alt != 'N'")
        if ignore_terminal_gaps:
            nuc_filter.append("alt != '.'")
        nuc_filter = " AND " + " AND ".join(nuc_filter) if nuc_filter else ""
        return nuc_filter, aa_filter

    def handle_csv_tsv_format(
        self,
        sample_selection_query: str,
        sample_selection_values: List[str],
        reference_accession: Optional[str],
        filter_n: bool,
        filter_x: bool,
        ignore_terminal_gaps: bool,
        output_columns: Optional[List[str]] = None,
    ) -> List[Dict[str, str]]:
        """
        Fetches data in csv/tsv format based on the given parameters.

        Args:
            sample_selection_query (str): SQL query to select the samples.
            sample_selection_values (List[str]): Values for property-based filtering.
            reference_accession (str): Reference accession to construct the conditions for the SQL query.
            filter_n (bool): Flag for handling 'N' variants.
            filter_x (bool): Flag for handling 'X' variants.
            ignore_terminal_gaps (bool): Flag for handling terminal gap variants.
            output_columns (Optional[List[str]]): List of output columns to include in the results.
                                                If not provided or empty, all columns are included.

        Returns:
            List[Dict[str, str]]: List of dictionaries where each dictionary represents a row in the resulting csv/tsv file.
        """
        # process property output
        property_cols = []
        for prop_name, prop_data in self.properties.items():
            prop_name = prop_name.lstrip(".")
            property_cols.append("p." + prop_name)
        property_cols_str = "\n, ".join(property_cols)

        # process mutation display filters
        nuc_filter, aa_filter = self.create_special_variant_filter_conditions(
            filter_n, filter_x, ignore_terminal_gaps
        )

        # set genomic and cds elements to consider
        genome_condition, molecule_prefix = self.create_genomic_element_conditions(
            reference_accession, element_alias="e", molecule_alias="m"
        )

        genome_condition = genome_condition.replace("e.id", "element_id")

        molecule_col = "m.symbol," if molecule_prefix else ""
        molecule_join = (
            "LEFT JOIN molecule m ON e.molecule_id = m.id" if molecule_prefix else ""
        )
        molecule_order = "m.symbol, " if molecule_prefix else ""

        cds_ids = ", ".join(
            [str(x) for x in self.get_element_ids(reference_accession, "cds")]
        )

        # assemble sql
        props_str = "rows." + ", rows.".join(
            sorted([s.lstrip(".") for s in self.properties.keys()])
        )

        sql = f"""SELECT
                name AS SAMPLE_NAME,
                {props_str},
                GROUP_CONCAT(CASE WHEN {genome_condition}{nuc_filter} THEN label END, ' ' ) AS GENOMIC_PROFILE,
                GROUP_CONCAT(CASE WHEN element_id IN ({cds_ids}){aa_filter} THEN {molecule_prefix}symbol || ':' || label END, ' ' ) AS PROTEOMIC_PROFILE,
                GROUP_CONCAT(CASE WHEN {genome_condition} AND frameshift = 1 THEN label END, ' ' ) AS FRAMESHIFT_MUTATIONS
                FROM
                (
                    SELECT
                        fs.sample_id,
                        fs.name,
                        fs.seqhash,
                        v.start,
                        v.label,
                        v.end,
                        v.ref,
                        v.alt,
                        v.element_id,
                        v.frameshift,
                        e.symbol,
                        {molecule_col}
                        {property_cols_str}
                    FROM (
                        {sample_selection_query}
                    ) AS fs
                    LEFT JOIN property p ON fs.sample_id = p.sample_id
                    LEFT JOIN alignment a ON fs.seqhash = a.seqhash
                    LEFT JOIN alignment2variant a2v ON a.id = a2v.alignment_id
                    LEFT JOIN variant v ON a2v.variant_id = v.id
                    LEFT JOIN element e ON v.element_id = e.id
                    {molecule_join}
                    ORDER BY {molecule_order}e.symbol, v.start
                ) as rows
                GROUP BY SAMPLE_NAME
                ORDER BY SAMPLE_NAME
                """
        sql = sqlparse.format(sql, reindent=True, keyword_case="upper")
        return self.cursor.execute(sql, sample_selection_values).fetchall()

    def handle_count_format(
        self, sample_selection_sql: str, sample_selection_values: List[str]
    ) -> int:
        """
        Count the distinct sample IDs in a selected sample set.

        Args:
            sample_selection_sql (str): The SQL command to select the samples.
            sample_selection_values (List[str]): Values for property-based filtering.

        Returns:
            int: The count of distinct sample IDs in the selected sample set.

        This function counts the number of distinct sample IDs in a selected sample set.
        It combines the sample selection SQL command with count operation and executes the resulting SQL command.
        """
        sql = f"SELECT COUNT(DISTINCT sample_id) AS count FROM ({sample_selection_sql})"
        result = self.cursor.execute(sql, sample_selection_values).fetchone()

        return result["count"] if result else 0

    def handle_vcf_format(
        self,
        sample_selection_sql: str,
        sample_selection_values: List[str],
        reference_accession: str,
        filter_n: str,
        ignore_terminal_gaps: str,
    ) -> List[Dict[str, str]]:
        """
        Fetches data in VCF format based on the given parameters.

        Args:
            sample_selection_sql (str): SQL query to select the samples.
            sample_selection_values (List[str]): Values for property-based filtering.
            reference_accession (str): Reference accession to construct the conditions for the SQL query.
            filter_n (str): Condition for handling 'N' variants.
            filter_terminal_gaps (str): Condition for handling terminal gap variants.
        Returns:
            List[Dict[str, str]]: List of dictionaries where each dictionary represents a row in the resulting VCF file.
        """
        genome_condition, _ = self.create_genomic_element_conditions(
            reference_accession, element_alias="e", molecule_alias="m"
        )
        nuc_filter, _ = self.create_special_variant_filter_conditions(
            filter_n, False, ignore_terminal_gaps
        )

        sql = f"""SELECT
                    rows.element_id AS "element.id",
                    rows.type AS "element.type",
                    rows.accession AS "molecule.accession",
                    rows.start AS "variant.start",
                    rows.end AS "variant.end",
                    rows.ref AS "variant.ref",
                    rows.alt AS "variant.alt",
                    rows.label AS "variant.label",
                    GROUP_CONCAT(rows.name) as samples
                FROM
                (
                    SELECT
                        fs.sample_id,
                        fs.name,
                        fs.seqhash,
                        v.start,
                        v.label,
                        v.end,
                        v.ref,
                        v.alt,
                        v.element_id,
                        e.type,
                        m.accession
                    FROM (
                        {sample_selection_sql}
                    ) AS fs
                    LEFT JOIN alignment a ON fs.seqhash = a.seqhash
                    LEFT JOIN alignment2variant a2v ON a.id = a2v.alignment_id
                    LEFT JOIN variant v ON a2v.variant_id = v.id
                    LEFT JOIN element e ON v.element_id = e.id
                    LEFT JOIN molecule m ON e.molecule_id = m.id
                    WHERE {genome_condition}{nuc_filter}
                ) as rows
                GROUP BY
                    "molecule.accession", "variant.start", "variant.ref", "variant.alt"
                ORDER BY
                    "molecule.accession", "variant.start"
        """
        return self.cursor.execute(sql, sample_selection_values)

    def match(
        self,
        profiles: Optional[Tuple[str, ...]] = None,
        samples: Optional[List[str]] = None,
        reference_accession: Optional[str] = None,
        properties: Optional[Dict[str, List[str]]] = None,
        frameshifts_only: Optional[bool] = False,
        filter_n: Optional[bool] = True,
        filter_x: Optional[bool] = True,
        ignore_terminal_gaps: Optional[bool] = True,
        output_columns: Optional[List[str]] = None,
        format: Optional[str] = "csv",
    ) -> Union[List[Dict[str, str]], str]:
        """
        This function matches samples based on their metadata and genomic profiles.

        Args:
            profiles (Optional[Tuple[str, ...]], optional): List of profiles to match. Defaults to None.
            samples (Optional[List[str]], optional): List of samples to consider. Defaults to None.
            reference_accession (str): Reference accession to construct the conditions for the SQL query.
            properties (Optional[Dict[str, List[str]]], optional): Dictionary of properties to filter samples. Defaults to None.
            frameshifts_only (Optional[bool], optional): Whether to consider only frameshifts. Defaults to False.
            filter_n (Optional[bool], optional): Whether to filter 'N'. Defaults to True.
            filter_x (Optional[bool], optional): Whether to filter 'X'. Defaults to True.
            ignore_terminal_gaps (Optional[bool], optional): Whether to ignore terminal gaps. Defaults to True.
            output_columns (Optional[List[str]], optional): List of output columns. Defaults to None.
            format (Optional[str], optional): Output format. It can be "csv", "tsv", "count", "vcf". Defaults to "csv".

        Returns:
            Union[List[Dict[str, str]], str]: Returns the matched samples in the requested format.

        Raises:
            SystemExit: If the provided output format is not supported.
        """

        (
            sample_selection_query,
            sample_selection_values,
        ) = self.create_sample_selection_sql(
            samples, properties, profiles, frameshifts_only
        )

        if format == "csv" or format == "tsv":
            return self.handle_csv_tsv_format(
                sample_selection_query,
                sample_selection_values,
                reference_accession,
                filter_n,
                filter_x,
                ignore_terminal_gaps,
                output_columns,
            )
        elif format == "count":
            return self.handle_count_format(
                sample_selection_query, sample_selection_values
            )
        elif format == "vcf":
            return self.handle_vcf_format(
                sample_selection_query,
                sample_selection_values,
                reference_accession,
                filter_n,
                ignore_terminal_gaps,
            )
        else:
            LOGGER.error(f"'{format}' is not a valid output format.")
            sys.exit(1)
