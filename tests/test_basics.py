import lzma
from pathlib import Path
import re

import pytest
from src.covsonar.basics import sonarBasics
from src.covsonar.utils import sonarUtils


def test_setup_and_file_exists(tmpfile_name):
    fname = tmpfile_name
    open(fname, "w").close()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        sonarUtils.setup_db(fname)
    assert pytest_wrapped_e.type == SystemExit
    assert (
        re.match(r"^setup error: .* does already exist.$", pytest_wrapped_e.value.code)
        is not None
    )


def test_autocreate_with_givenRef(tmpfile_name, monkeypatch):
    monkeypatch.chdir(Path(__file__).parent)
    ref = Path("data/ref.gb")
    assert (
        sonarUtils.setup_db(tmpfile_name, default_props=True, reference_gb=ref) is None
    )


def test_notautocreate_and_notsetup(tmpfile_name):
    fname = tmpfile_name
    assert sonarUtils.setup_db(fname, reference_gb=None, default_props=False) is None


def test_basicObject():
    obj = sonarUtils()
    assert isinstance(obj, sonarUtils) is True


def test_match(testdb):
    # "Wrong outputformat."
    with pytest.raises(SystemExit) as e:
        sonarUtils().match(
            testdb,
            format="HDFS",
            debug="False",
        )
    assert e.type == SystemExit
    assert "is not a valid output format" in e.value.code


def test_import_data(testdb, monkeypatch):
    """ """
    monkeypatch.chdir(Path(__file__).parent)
    tsv = "data/meta.tsv"
    fasta = "data/test.fasta"
    with pytest.raises(SystemExit) as e:
        sonarUtils().import_data(
            testdb,
            fasta=[],
            tsv_files=[],
            prop_links={},
            cachedir=None,
            autolink=False,
            progress=True,
            update=True,
            threads=1,
            debug=True,
            quiet=True,
        )
    # "Nothing to import."
    assert e.type == SystemExit
    assert e.value.code == 0

    # " not a valid property assignment."
    with pytest.raises(SystemExit) as e:
        sonarUtils().import_data(
            testdb,
            fasta=[fasta],
            tsv_files=[],
            prop_links=["STUDY", "STE"],
            autolink=False,
            progress=True,
            update=True,
            debug=True,
            quiet=True,
        )
    assert e.type == SystemExit
    assert (
        e.value.code
        == "input error: STUDY is not a valid column-to-property assignment."
    )

    # unknown column to the selected DB
    with pytest.raises(SystemExit) as e:
        sonarUtils().import_data(
            testdb,
            fasta=[fasta],
            tsv_files=[],
            prop_links=["STE=STUDY", "sample=IMS_ID"],
            autolink=False,
            progress=True,
            update=True,
            debug=True,
            quiet=True,
        )
    assert e.type == SystemExit
    assert "is unknown to the selected database" in e.value.code

    # sample column has to be assigned.
    with pytest.raises(SystemExit) as e:
        sonarUtils().import_data(
            testdb,
            fasta=[],
            tsv_files=[tsv],
            prop_links=[],
            cachedir=None,
            autolink=False,
            progress=True,
            update=True,
            debug=True,
            quiet=True,
        )
    assert e.type == SystemExit
    assert "input error: missing 'sample' column assignment." in e.value.code

    # tsv does not provide any informative column
    with pytest.raises(SystemExit) as e:
        sonarUtils().import_data(
            testdb,
            fasta=[],
            tsv_files=[tsv],
            prop_links=["sample=IMS_ID"],
            autolink=False,
        )
    assert e.type == SystemExit
    assert (
        "input error: the file does not provide any informative column." == e.value.code
    )

    # is not an unique column
    bad_tsv = "data/meta.bad.tsv"
    with pytest.raises(SystemExit) as e:
        sonarUtils().import_data(
            testdb,
            fasta=[],
            tsv_files=[bad_tsv],
            prop_links=["sample=IMS_ID", "SEQ_TYPE=SEQ_TYPE"],
            autolink=False,
        )
    assert e.type == SystemExit
    assert "error: SEQ_TYPE is not a unique column." == e.value.code

    # in 'tsv file does not contain required sample column.'
    with pytest.raises(SystemExit) as e:
        sonarUtils().import_data(
            testdb,
            fasta=[],
            tsv_files=[tsv],
            prop_links=["sample=INVISIBLE_COl"],
            autolink=True,
        )
    assert e.type == SystemExit
    assert "input error: missing 'sample' column assignment." == e.value.code


def test_open_file(tmpfile_name):
    # create dump file. xz
    xz_file = tmpfile_name + ".xz"
    msg = "The quick brown fox jumps over the lazy dog"
    with lzma.open(xz_file, "wb") as iFile:
        iFile.write(msg.encode())

    # decompress data
    with sonarBasics.open_file_autodetect(xz_file) as handle:
        assert handle.read() == msg
