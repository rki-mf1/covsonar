from sonar import sonarBasics, sonarDBManager, sonarCache, sonarLinmgr, sonarAligner
import pytest
import tempfile

# PYTEST FIXTURES
@pytest.fixture
def setup_testdb(tmpdir_factory, scope="session"):
	""" Fixture to set up a temporay session-lasting test database with test data """
	dbfile = str(tmpdir_factory.mktemp('dbm_test').join("testdb"))
	sonarBasics.setup_db(dbfile, quiet=True)
	with sonarDBManager(dbfile, readonly=False) as dbm:
		dbm.add_property("LINEAGE", "text", "text", " ")
	return dbfile

@pytest.fixture
def setup_db(get_tmpfile_name, scope="session"):
	""" Fixture to set up a temporay test database with test data """
	sonarBasics.setup_db(get_tmpfile_name, quiet=True)
	return get_tmpfile_name

@pytest.fixture
def testdb(setup_testdb, scope="session"):
	""" Fixture to set up a temporay test database with test data """
	return setup_testdb

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

@pytest.fixture
def get_tmpfile_name(tmpdir_factory):
	return str(tmpdir_factory.mktemp('dbm_test').join(next(tempfile._get_candidate_names())))
