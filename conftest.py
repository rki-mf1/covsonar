import os
import tempfile

import pytest
from src.covsonar.utils import sonarUtils
from src.covsonar.dbm import sonarDBManager


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
def setup_db(tmp_path_factory):
    """Fixture to set up a temporay session-lasting test database with test data"""
    dbfile = str(tmp_path_factory.mktemp("data") / "test.db")
    sonarUtils.setup_db(dbfile, quiet=True)
    with sonarDBManager(dbfile, readonly=False) as dbm:
        dbm.add_property("LINEAGE", "text", "text", " ", "sample")
        dbm.add_property("CITY", "zip", "zip", " ", "sample")
        dbm.add_property("Ct", "integer", "integer", " ", "sample")
        dbm.add_property("CONC", "float", "float", " ", "sample")
        dbm.add_property("SAMPLING", "date", "date", " ", "sample")
        dbm.add_property("SEQ_TYPE", "text", "text", " ", "sample")
    return dbfile


@pytest.fixture
def tmpfile_name(tmpdir_factory):
    yield str(
        tmpdir_factory.mktemp("dbm_test").join(next(tempfile._get_candidate_names()))
    )


@pytest.fixture(scope="session")
def testdb(setup_db):
    db = setup_db
    script_dir = os.path.dirname(os.path.abspath(__file__))
    fasta = os.path.join(script_dir, "tests", "data", "test.fasta")
    sonarUtils.import_data(
        db,
        fasta=[fasta],
        tsv_files=[],
        prop_links={},
        cachedir=None,
        autolink=True,
        progress=False,
        update=True,
        debug=False,
        quiet=True,
    )
    return db


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
