import tempfile

import pytest

from covsonar.basics import sonarBasics
from covsonar.dbm import sonarDBManager


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


@pytest.fixture(scope="session")
def setup_testdb(tmpdir_factory):
    """Fixture to set up a temporay session-lasting test database with test data"""
    dbfile = str(tmpdir_factory.mktemp("dbm_test").join("testdb"))
    # sonarBasics.setup_db(dbfile, quiet=True)
    sonarBasics.setup_db(dbfile)
    with sonarDBManager(dbfile, readonly=False) as dbm:
        dbm.add_property("LINEAGE", "text", "text", " ")
    return dbfile


@pytest.fixture
def setup_db(get_tmpfile_name):
    """Fixture to set up a temporay test database with test data"""
    # sonarBasics.setup_db(get_tmpfile_name, quiet=True)
    sonarBasics.setup_db(get_tmpfile_name)
    return get_tmpfile_name


@pytest.fixture
def testdb(setup_testdb):
    """Fixture to set up a temporay test database with test data"""
    return setup_testdb


@pytest.fixture
def init_readonly_dbm(testdb):
    """Fixture to set up a read-only dbm object"""
    with sonarDBManager(testdb, readonly=True) as dbm:
        yield dbm


@pytest.fixture
def init_writeable_dbm(testdb):
    """Fixture to set up a wirte-able dbm object"""
    with sonarDBManager(testdb, readonly=False) as dbm:
        yield dbm


@pytest.fixture
def get_tmpfile_name(tmpdir_factory):
    return str(
        tmpdir_factory.mktemp("dbm_test").join(next(tempfile._get_candidate_names()))
    )
