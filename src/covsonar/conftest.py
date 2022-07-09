from sonar import sonarBasics, sonarDBManager, sonarCache, sonarLinmgr, sonarAligner
import pytest
import tempfile
import os

# PYTEST FIXTURES
@pytest.fixture(autouse=True)
def mock_workerpool_imap_unordered(monkeypatch):
    """Mock mpire's WorkerPool.imap_unordered function
    This is necessary to work around crashes caused by trying to calculate
    coverage with multiprocess subprocesses, and also to make the tests
    reproducible (ordered).
    """
    monkeypatch.setattr(
        "mpire.WorkerPool.imap_unordered",
        lambda self, func, args=(), kwds={}, callback=None, error_callback=None: (
            func(arg) for arg in args
        ),
    )

@pytest.fixture
def setup_db(tmpfile_name, scope="session"):
	""" Fixture to set up a temporay session-lasting test database with test data """
	dbfile = tmpfile_name
	sonarBasics.setup_db(dbfile, quiet=True)
	with sonarDBManager(dbfile, readonly=False) as dbm:
		dbm.add_property("LINEAGE", "text", "text", " ")
		dbm.add_property("CITY", "zip", "zip", " ")
		dbm.add_property("Ct", "integer", "integer", " ")
		dbm.add_property("CONC", "float", "float", " ")
		dbm.add_property("SAMPLING", "date", "date", " ")
	return dbfile

@pytest.fixture
def tmpfile_name(tmpdir_factory, scope="session"):
	yield str(tmpdir_factory.mktemp('dbm_test').join(next(tempfile._get_candidate_names())))

@pytest.fixture
def testdb(setup_db, scope="session"):
	db = setup_db
	script_dir = os.path.dirname(os.path.abspath(__file__))
	fasta = os.path.join(script_dir, "..", "test", "test.fasta")
	sonarBasics.import_data(db, fasta=[fasta], tsv=[], cols={}, cachedir = None, autodetect=True, progress=False, update=True, debug=False, quiet=True)
	with sonarDBManager(db, readonly=True) as dbm:
			pass
			#exit("seq count " + str(dbm.count_sequences()))
	return db

@pytest.fixture
def init_readonly_dbm(testdb, scope="session"):
	""" Fixture to set up a read-only dbm object """
	with sonarDBManager(testdb, readonly=True) as dbm:
		yield dbm

@pytest.fixture
def init_writeable_dbm(testdb, scope="session"):
	""" Fixture to set up a wirte-able dbm object """
	with sonarDBManager(testdb, readonly=False) as dbm:
		yield dbm