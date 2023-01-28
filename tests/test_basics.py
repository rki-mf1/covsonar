import lzma
from pathlib import Path
import re

import pytest
from src.covsonar.basics import sonarBasics


def test_setup_and_file_exists(tmpfile_name):
    fname = tmpfile_name
    open(fname, "w").close()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        sonarBasics.setup_db(fname)
    assert pytest_wrapped_e.type == SystemExit
    assert (
        re.match(r"^setup error: .* does already exist.$", pytest_wrapped_e.value.code)
        is not None
    )


# def test_autocreate_with_givenRef(tmpfile_name):
#    fname = tmpfile_name
#    # ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
#    ref = pkg_resources.resource_filename("covsonar",'data/ref.gb')
#    print(ref)
#    assert sonarBasics.setup_db(fname,auto_create=True,reference_gb=ref) is None


def test_notautocreate_and_notsetup(tmpfile_name):
    fname = tmpfile_name
    assert sonarBasics.setup_db(fname, auto_create=False, default_setup=False) is None


def test_basicObject():
    obj = sonarBasics()
    assert isinstance(obj, sonarBasics) is True


def test_match(testdb):
    # "Wrong outputformat."
    with pytest.raises(SystemExit) as e:
        sonarBasics().match(
            testdb,
            profiles=[],
            reserved_props_dict={},
            propdict={},
            reference=None,
            outfile=None,
            format="HDFS",
            debug="False",
        )
    assert e.type == SystemExit
    assert "is not a valid output format" in e.value.code


def test_import_data(testdb, monkeypatch):
    """ """
    monkeypatch.chdir(Path(__file__).parent)
    tsv = "data/meta.tsv"
    # fasta = "data/test.fasta"
    with pytest.raises(SystemExit) as e:
        sonarBasics().import_data(
            testdb,
            fasta=[],
            tsv=[],
            cols={},
            cachedir=None,
            autodetect=False,
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
        sonarBasics().import_data(
            testdb,
            fasta=[],
            tsv=[""],
            cols=["STUDY", "STE"],
            autodetect=False,
            progress=True,
            update=True,
            debug=True,
            quiet=True,
        )
    assert e.type == SystemExit
    assert "STUDY is not a valid sample property assignment" in e.value.code

    # unknown column to the selected DB
    with pytest.raises(SystemExit) as e:
        sonarBasics().import_data(
            testdb,
            fasta=[],
            tsv=[""],
            cols=["STE=STUDY", "sample=IMS_ID"],
            autodetect=False,
            progress=True,
            update=True,
            debug=True,
            quiet=True,
        )
    assert e.type == SystemExit
    assert "is unknown to the selected database" in e.value.code

    # sample column has to be assigned.
    with pytest.raises(SystemExit) as e:
        sonarBasics().import_data(
            testdb,
            fasta=[],
            tsv=[tsv],
            cols=[],
            cachedir=None,
            autodetect=False,
            progress=True,
            update=True,
            debug=True,
            quiet=True,
        )
    assert e.type == SystemExit
    assert "sample column has to be assigned" in e.value.code

    # tsv does not provide any informative column
    with pytest.raises(SystemExit) as e:
        sonarBasics().import_data(
            testdb, fasta=[], tsv=[tsv], cols=["sample=IMS_ID"], autodetect=False
        )
    assert e.type == SystemExit
    assert "tsv does not provide any informative column" in e.value.code

    # is not an unique column
    bad_tsv = "data/meta.bad.tsv"
    with pytest.raises(SystemExit) as e:
        sonarBasics().import_data(
            testdb,
            fasta=[],
            tsv=[bad_tsv],
            cols=["sample=IMS_ID", "SEQ_TYPE=SEQ_TYPE"],
            autodetect=False,
        )
    assert e.type == SystemExit
    assert "is not an unique column" in e.value.code

    # in 'tsv file does not contain required sample column.'
    with pytest.raises(SystemExit) as e:
        sonarBasics().import_data(
            testdb, fasta=[], tsv=[tsv], cols=["sample=INVISIBLE_COl"], autodetect=True
        )
    assert e.type == SystemExit
    assert "tsv file does not contain required sample column" in e.value.code


def test_open_file(tmpfile_name):

    # input error: " + fname + " does not exist.
    with pytest.raises(SystemExit) as e:
        sonarBasics().open_file(tmpfile_name, compressed="xz")
    assert e.type == SystemExit
    assert "does not exist" in e.value.code

    # create dump file. xz
    xz_file = tmpfile_name + ".xz"
    msg = "The quick brown fox jumps over the lazy dog"
    with lzma.open(xz_file, "wb") as iFile:
        iFile.write(msg.encode())

    # decompress data
    answer = sonarBasics().open_file(xz_file, compressed="xz")
    assert answer.read() == msg

    # "input error: " + fname + " cannot be opened."
    with pytest.raises(SystemExit) as e:
        answer = sonarBasics().open_file(xz_file, compressed="zip", encoding="utf-2022")
    assert e.type == SystemExit
    assert "cannot be opened" in e.value.code
