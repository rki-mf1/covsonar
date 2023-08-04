import logging
import lzma
from pathlib import Path
import re

import pytest
from src.covsonar.basics import sonarBasics
from src.covsonar.utils import sonarUtils


def test_setup_and_file_exists(logger, caplog, tmpfile_name):
    fname = tmpfile_name
    open(fname, "w").close()
    with caplog.at_level(logging.ERROR, logger=logger.name), pytest.raises(
        SystemExit
    ) as pytest_wrapped_e:
        sonarUtils.setup_db(fname)
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 1
        assert (
            re.match(
                r".* does already exist.$",
                sonarBasics.get_last_log_entry(caplog),
            )
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


def test_match(caplog, testdb, logger):
    # "Wrong outputformat."
    with caplog.at_level(logging.ERROR, logger=logger.name), pytest.raises(
        SystemExit
    ) as pytest_wrapped_e:
        sonarUtils().match(
            testdb,
            format="HDFS",
        )
        assert "'HDFS' is not a valid output format." == caplog.records[-1].message
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 1


def test_import_data(testdb, monkeypatch, logger, caplog):
    """ """
    monkeypatch.chdir(Path(__file__).parent)
    tsv = "data/meta.tsv"
    fasta = "data/test.fasta"
    # "Nothing to import."
    with caplog.at_level(logging.ERROR, logger=logger.name), pytest.raises(
        SystemExit
    ) as pytest_wrapped_e:
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
            quiet=True,
        )
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 0

    # " not a valid property assignment."
    with caplog.at_level(logging.ERROR, logger=logger.name), pytest.raises(
        SystemExit
    ) as pytest_wrapped_e:
        sonarUtils().import_data(
            testdb,
            fasta=[fasta],
            tsv_files=[],
            prop_links=["STUDY", "STE"],
            autolink=False,
            progress=True,
            update=True,
            quiet=True,
        )
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 1
        assert (
            "'STUDY' is not a valid column-to-property assignment."
            == caplog.records[-1].message
        )

    # unknown column to the selected DB
    with caplog.at_level(logging.ERROR, logger=logger.name), pytest.raises(
        SystemExit
    ) as pytest_wrapped_e:
        sonarUtils().import_data(
            testdb,
            fasta=[fasta],
            tsv_files=[],
            prop_links=["STE=STUDY", "sample=IMS_ID"],
            autolink=False,
            progress=True,
            update=True,
            quiet=True,
        )
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 1
        assert (
            "Sample property 'STUDY' is unknown to the selected database. Use list-props to see all valid properties"
            == caplog.records[-1].message
        )

    # sample column has to be assigned.
    with caplog.at_level(logging.ERROR, logger=logger.name), pytest.raises(
        SystemExit
    ) as pytest_wrapped_e:
        sonarUtils().import_data(
            testdb,
            fasta=[],
            tsv_files=[tsv],
            prop_links=[],
            cachedir=None,
            autolink=False,
            progress=True,
            update=True,
            quiet=True,
        )
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 1
        assert "Missing 'sample' column assignment." == caplog.records[-1].message

    # tsv does not provide any informative column
    with caplog.at_level(logging.ERROR, logger=logger.name), pytest.raises(
        SystemExit
    ) as pytest_wrapped_e:
        sonarUtils().import_data(
            testdb,
            fasta=[],
            tsv_files=[tsv],
            prop_links=["sample=IMS_ID"],
            autolink=False,
        )
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 1
        assert (
            "The file does not provide any informative column."
            == caplog.records[-1].message
        )

    # is not an unique column
    bad_tsv = "data/meta.bad.tsv"
    with caplog.at_level(logging.ERROR, logger=logger.name), pytest.raises(
        SystemExit
    ) as pytest_wrapped_e:
        sonarUtils().import_data(
            testdb,
            fasta=[],
            tsv_files=[bad_tsv],
            prop_links=["sample=IMS_ID", "SEQ_TYPE=SEQ_TYPE"],
            autolink=False,
        )
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 1
        assert "'SEQ_TYPE' is not a unique column." == caplog.records[-1].message

    # in 'tsv file does not contain required sample column.'
    with caplog.at_level(logging.ERROR, logger=logger.name), pytest.raises(
        SystemExit
    ) as pytest_wrapped_e:
        sonarUtils().import_data(
            testdb,
            fasta=[],
            tsv_files=[tsv],
            prop_links=["sample=INVISIBLE_COl"],
            autolink=True,
        )

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 1
        assert "Missing 'sample' column assignment." == caplog.records[-1].message


def test_open_file(tmpfile_name):
    # create dump file. xz
    xz_file = tmpfile_name + ".xz"
    msg = "The quick brown fox jumps over the lazy dog"
    with lzma.open(xz_file, "wb") as iFile:
        iFile.write(msg.encode())

    # decompress data
    with sonarBasics.open_file_autodetect(xz_file) as handle:
        assert handle.read() == msg
