import os
from pathlib import Path

import pytest

from covsonar import sonar


def test_check_file_not_exist(tmpfile_name):
    fname = "no/real/file/existence"
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        sonar.check_file(fname)
    assert pytest_wrapped_e.type == SystemExit
    assert (
        pytest_wrapped_e.value.code == "iput error: " + fname + " is not a valid file."
    )


def test_check_file_exist(monkeypatch):
    monkeypatch.chdir(Path(__file__).parent)
    fname = "data/test.fasta"
    assert sonar.check_file(fname) is None


def test_check_db_not_exist(monkeypatch):
    monkeypatch.chdir(Path(__file__).parent)
    parsed_args = sonar.parse_args(["db-upgrade", "--db", "/no/db/is/here.db"])
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        sonar.main(parsed_args)
    assert pytest_wrapped_e.type == SystemExit


def test_unknown_properties(monkeypatch):
    monkeypatch.chdir(Path(__file__).parent)
    db_path = "data/test-with-seqs.db"
    parsed_args = sonar.parse_args(
        ["delete-prop", "--db", db_path, "--name", "WHAT_IS_THIS"]
    )
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        sonar.main(parsed_args)
    assert pytest_wrapped_e.type == SystemExit


def test_match_with_sample(monkeypatch):
    monkeypatch.chdir(Path(__file__).parent)
    db_path = "data/test-with-seqs.db"
    sample_file = "data/sample_list.txt"
    parsed_args = sonar.parse_args(
        [
            "match",
            "--db",
            db_path,
            "--sample",
            "IMS-10025-CVDP-00960",
            "--sample-file",
            sample_file,
            "--count",
        ]
    )
    assert sonar.main(parsed_args) == 0


def test_delete_nothing(tmp_path, monkeypatch, caplog):
    monkeypatch.chdir(Path(__file__).parent)
    db_path = "data/test-with-seqs.db"
    _tmp_sample_file = os.path.join(tmp_path, "_tmp_sample_file.txt")
    with open(_tmp_sample_file, "w") as file:
        file.write("NoOneKnowWhatYouAreLookingFor")

    parsed_args = sonar.parse_args(
        [
            "delete",
            "--db",
            db_path,
            "--sample",
            "nothing_to_delete",
            "--sample-file",
            _tmp_sample_file,
        ]
    )
    # with caplog.at_level(logging.INFO):
    assert sonar.main(parsed_args) == 0
    # assert '0 of 2 samples found and deleted. 5 samples remain in the database.' in caplog.text

    parsed_args = sonar.parse_args(
        [
            "delete",
            "--db",
            db_path,
        ]
    )
    # with caplog.at_level(logging.INFO):
    assert sonar.main(parsed_args) == 0
    # assert 'Nothing to delete.' in caplog.text
