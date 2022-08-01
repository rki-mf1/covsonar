import os
from pathlib import Path
import shutil

from covsonar.cache import sonarCache


def test_get_refseq_id_failcase(tmpfile_name, testdb):
    # We use query condition where a record is nonexistent.
    result = sonarCache(db=testdb).get_refseq_id(refmol_acc=tmpfile_name)

    assert result is None


def test_get_refhash_failcase(tmpfile_name, testdb):
    # We use query condition where a record is nonexistent.
    result = sonarCache(db=testdb).get_refhash(refmol_acc=tmpfile_name)

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
    db_path = os.path.join(tmp_path, "import-test.db")

    shutil.copy(testdb, db_path)
    # no update
    result = sonarCache(db=db_path, allow_updates=False).add_fasta(fasta)
    assert result is None
    # Update
    result = sonarCache(db=db_path, allow_updates=True).add_fasta(fasta)
    assert result is None
