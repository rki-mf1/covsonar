import logging
import os
from pathlib import Path
import shutil

import pytest
from src.covsonar.cache import sonarCache
from src.covsonar.dbm import sonarDbManager


def test_get_refseq_id_failcase(tmpfile_name, testdb):
    # We use query condition where a record is nonexistent.
    # with pytest.raises(Exception) as e:
    with sonarCache(db=testdb, logfile=tmpfile_name) as sc:
        result = sc.get_refseq_id(refmol_acc=type("Test"))
    assert result is None


def test_write_checkref_log(testdb, tmpfile_name, tmp_path, logger, caplog):
    data = {"header": "HEader", "name": "NAme"}

    # ignore_errors
    with sonarCache(db=testdb, ignore_errors=True, outdir=tmp_path) as sc:
        with caplog.at_level(logging.WARNING, logger=logger.name):
            sc.write_checkref_log(data=data, refseq_id=None)
            assert (
                "Skipping 'NAme' referring to an unknown reference ('HEader')."
                == caplog.records[-1].message
            )

    # not ignore_errors
    with sonarCache(
        db=testdb, ignore_errors=False, logfile=tmpfile_name, outdir=tmp_path
    ) as sc:
        with caplog.at_level(logging.ERROR, logger=logger.name), pytest.raises(
            SystemExit
        ) as pytest_wrapped_e:
            sc.write_checkref_log(data=data, refseq_id=None)
            assert (
                "ERROR: Fasta header refers to an unknown reference ('HEader')."
                == caplog.records[-1].message
            )
            assert pytest_wrapped_e.type == SystemExit
            assert pytest_wrapped_e.value.code == 1


def test_get_refhash_failcase(tmpfile_name, tmp_path, testdb):
    # We use query condition where a record is nonexistent.

    with sonarCache(db=testdb, outdir=tmp_path) as sc:
        result = sc.get_refhash(refmol_acc=tmpfile_name)
    assert result is None


def test_get_add_fasta(tmp_path, monkeypatch):
    """
    This add_fasta function is quite complex and hard to perform the unit test;
    we might consider rewriting this function.
    """
    monkeypatch.chdir(Path(__file__).parent)
    # tsv = "data/meta.tsv"
    testdb = "data/test-with-seqs.db"
    fasta = "data/test.fasta"
    # bad_fasta = "data/bad.fasta"
    db_path = os.path.join(tmp_path, "import-test.db")

    shutil.copy(testdb, db_path)
    # no update
    result = sonarCache(db=db_path, allow_updates=False).add_fasta(fasta)
    assert result is None
    # Update
    result = sonarCache(db=db_path, allow_updates=True).add_fasta(fasta)
    assert result is None

    # fail
    # result = sonarCache(db=db_path, allow_updates=True).add_fasta(bad_fasta)
    # assert result is None


def test_add_data_files(monkeypatch):
    monkeypatch.chdir(Path(__file__).parent)
    testdb = "data/test-with-seqs.db"
    data = {
        "seqhash": "wVHEydRZjNKWn3ubf/+OMXxQ1WQ",
        "algnid": 1,
        "refseq_id": 0,
        "sequence": "wVHEydRZjNKWn3ubf/+OMXxQ1WQ",
    }
    seqhash = "Test201"
    with sonarDbManager(testdb) as dbm:
        return_data = sonarCache(db=testdb).add_data_files(data, seqhash, dbm)

    logging.debug(return_data)
    data1_truth = {
        "reffile": None,
        "ttfile": None,
        "liftfile": None,
        "algnfile": None,
        "varfile": None,
    }

    assert set(data1_truth.items()).issubset(set(return_data.items()))

    # seqhash !=
    data = {
        "seqhash": "Test201",
        "algnid": 1,
        "refseq_id": 0,
        "sequence": "wVHEydRZjNKWn3ubf/+OMXxQ1WQ",
    }
    with sonarDbManager(testdb) as dbm:
        return_data = sonarCache(db=testdb).add_data_files(data, seqhash, dbm)
    data2_truth = {
        "seqfile": None,
        "seqhash": None,
        "reffile": None,
        "ttfile": None,
        "liftfile": None,
        "algnfile": None,
        "varfile": None,
    }
    assert set(data2_truth.items()).issubset(set(return_data.items()))


def test_paranoid_test(monkeypatch, testdb, logger, caplog):
    """
    # passed case needed
    with open("data/cache-test/ref-seq-passed") as f:
        data = f.read()
    refseqs= eval(data)
    with open("data/cache-test/sample-data-passed") as f:
        data = f.read()
    sample_data = eval(data)

    with sonarDbManager(testdb) as dbm:
        _return = sonarCache(db=testdb).paranoid_test(refseqs=refseqs,sample_data=sample_data,dbm=dbm)
    """
    monkeypatch.chdir(Path(__file__).parent)

    # fail case
    with open("data/cache-test/ref-seq") as f:
        data = f.read()
    refseqs = eval(data)
    with open("data/cache-test/sample-data") as f:
        data = f.read()
    sample_data = eval(data)

    with sonarDbManager(testdb) as dbm, caplog.at_level(
        logging.ERROR, logger=logger.name
    ), pytest.raises(SystemExit) as pytest_wrapped_e:
        sonarCache(db=testdb).paranoid_test(
            refseqs=refseqs, sample_data=sample_data, dbm=dbm
        )
        assert "ERROR: Paranoid test failed." == caplog.records[-1].message
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 1

    os.remove("paranoid.alignment.fna")


def test_cache_sequence(testdb, monkeypatch, logger, caplog):
    # fail case
    monkeypatch.chdir(Path(__file__).parent)

    seqhash = "84nCd4KD+0OVpmoZCSDTnc/mn08"
    sequence = "test"
    with caplog.at_level(logging.ERROR, logger=logger.name), pytest.raises(
        SystemExit
    ) as pytest_wrapped_e:
        sonarCache(db=testdb, outdir="data/cache-test").cache_sequence(
            seqhash, sequence
        )
        assert "ERROR: Sequences differ for seqhash" == caplog.records[-1].message
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 1
