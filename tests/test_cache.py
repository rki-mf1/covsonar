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


def test_write_checkref_log(testdb, tmpfile_name, tmp_path, capsys):
    # tmpfile_name = "./tesst.txt"
    data = {"header": "HEader", "name": "NAme"}
    # ignore_errors
    with sonarCache(db=testdb, ignore_errors=True, outdir=tmp_path) as sc:
        sc.write_checkref_log(data=data, refseq_id=None)
        captured = capsys.readouterr()
    assert captured.err == "skipping NAme referring to an unknown reference (HEader)"

    # not ignore_errors
    with sonarCache(
        db=testdb, ignore_errors=False, logfile=tmpfile_name, outdir=tmp_path
    ) as sc:
        sc.write_checkref_log(data=data, refseq_id=None)


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


def test_paranoid_test(monkeypatch, testdb):
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

    with pytest.raises(ValueError) as pytest_wrapped_e:
        with sonarDbManager(testdb) as dbm:
            sonarCache(db=testdb).paranoid_test(
                refseqs=refseqs, sample_data=sample_data, dbm=dbm
            )
    assert pytest_wrapped_e.type == ValueError
    assert "cannot be restored" in str(pytest_wrapped_e.value)
    os.remove("paranoid.alignment.fna")


def test_cache_sequence(testdb, monkeypatch):
    # fail case
    monkeypatch.chdir(Path(__file__).parent)

    seqhash = "84nCd4KD+0OVpmoZCSDTnc/mn08"
    sequence = "test"
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        sonarCache(db=testdb, outdir="data/cache-test").cache_sequence(
            seqhash, sequence
        )
    assert pytest_wrapped_e.type == SystemExit
    assert "sequences differ for seqhash" in pytest_wrapped_e.value.code
