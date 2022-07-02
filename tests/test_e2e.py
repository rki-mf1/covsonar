from pathlib import Path
import re

import pytest

from covsonar import sonar


def split_cli(s):
    """Split a string into a list of individual arguments, respecting quotes"""
    return re.findall(r'(?:[^\s,"]|"(?:\\.|[^"])*")+', s)


def run_cli(s):
    """Helper function to simulate running the command line ./sonar <args>"""
    return sonar.main(sonar.parse_args(split_cli(s)))


def test_help():
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        parsed_args = sonar.parse_args(["--help"])
        sonar.main(parsed_args)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0


def test_setup_db(tmp_path):
    parsed_args = sonar.parse_args(["setup", "--db", str(tmp_path / "test.db")])
    retval = sonar.main(parsed_args)
    assert retval == 0


# The following two tests run the commands from the example test script, but
# are split in two to skip the "import" command which is currently causing a
# crash when run via pytest
def test_valid_beginning(tmp_path, monkeypatch):
    """The test example provided by other devs, up to the import command"""
    monkeypatch.chdir(Path(__file__).parent)

    db_path = str(tmp_path / "test.db")
    sonar.main(sonar.parse_args(split_cli(f"setup --db {db_path}")))

    run_cli(f"add-prop --db {db_path} --name SENDING_LAB --dtype integer --descr descr")
    run_cli(f"add-prop --db {db_path} --name DATE_DRAW --dtype date --descr descr")
    run_cli(f"add-prop --db {db_path} --name SEQ_TYPE --dtype text --descr descr")
    run_cli(f"add-prop --db {db_path} --name SEQ_REASON --dtype text --descr descr")
    run_cli(f"add-prop --db {db_path} --name SAMPLE_TYPE --dtype text --descr descr")
    run_cli(f"add-prop --db {db_path} --name OWN_FASTA_ID --dtype text --descr descr")
    run_cli(f"add-prop --db {db_path} --name DOWNLOAD_ID --dtype text --descr descr")
    run_cli(f"add-prop --db {db_path} --name DEMIS_ID --dtype integer --descr descr")
    run_cli(f"add-prop --db {db_path} --name RECEIVE_DATE --dtype date --descr descr")
    run_cli(
        f"add-prop --db {db_path} --name PROCESSING_DATE --dtype date --descr descr"
    )
    run_cli(
        f"add-prop --db {db_path} --name PUBLICATION_STATUS --dtype text --descr descr"
    )
    run_cli(
        f"add-prop --db {db_path} --name HASHED_SEQUENCE --dtype text --descr descr"
    )
    run_cli(f"add-prop --db {db_path} --name TIMESTAMP --dtype text --descr descr")
    run_cli(f"add-prop --db {db_path} --name STUDY --dtype text --descr descr")
    run_cli(
        f"add-prop --db {db_path} --name DOWNLOADING_TIMESTAMP --dtype text --descr descr"
    )
    run_cli(f"add-prop --db {db_path} --name SENDING_LAB_PC --dtype zip --descr descr")
    run_cli(f"add-prop --db {db_path} --name DEMIS_ID_PC --dtype zip --descr descr")
    run_cli(f"add-prop --db {db_path} --name VERSION --dtype integer --descr descr")
    run_cli(f"add-prop --db {db_path} --name DESH_QC_PASSED --dtype text --descr descr")
    run_cli(
        f"add-prop --db {db_path} --name DESH_REJECTION_REASON --dtype text --descr descr"
    )
    run_cli(f"add-prop --db {db_path} --name DUPLICATE_ID --dtype text --descr descr")
    run_cli(f"add-prop --db {db_path} --name LINEAGE --dtype text --descr descr")

    run_cli(f"update-lineage-info --db {db_path}")


def test_valid_end(tmp_path, monkeypatch):
    """The test example provided by other devs, after the import command"""
    monkeypatch.chdir(Path(__file__).parent)

    db_path = "data/test.db"

    run_cli(f"match --db {db_path} --profile S:A67G --DEMIS_ID 10013")
    run_cli(
        f"match --db {db_path} --DATE_DRAW 2021-11-01:2022-12-15 -o {tmp_path}/temp.tsv"
    )
    run_cli(f"match --db {db_path} --LINEAGE B.1.1.7 --with-sublineage LINEAGE --count")
