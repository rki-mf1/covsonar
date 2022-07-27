import filecmp
from pathlib import Path
import re
import shutil

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
    parsed_args = sonar.parse_args(["setup", "--db", str(tmp_path / "test.db"), "-a"])
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


def test_import(tmp_path, monkeypatch):
    """The test example provided by other devs, after the import command"""
    monkeypatch.chdir(Path(__file__).parent)

    db_path_orig = Path("data/test.db")
    db_path = tmp_path / "import-test.db"

    shutil.copy(db_path_orig, db_path)

    run_cli(
        f"import --db {db_path} --fasta data/seqs.fasta.gz --tsv data/meta.tsv --cache {tmp_path} --cols sample=IMS_ID --threads 2"
    )


def test_valid_end(tmp_path, monkeypatch):
    """The test example provided by other devs, after the import command"""
    monkeypatch.chdir(Path(__file__).parent)

    db_path_orig = Path("data/test-with-seqs.db")
    db_path = tmp_path / "test-with-seqs.db"
    shutil.copy(db_path_orig, db_path)
    run_cli(f"match --db {db_path} --profile ^A3451T A3451TGAT ")
    run_cli(f"match --db {db_path} --profile del:28363-28371  --profile A3451N ")
    run_cli(f"match --db {db_path} --profile ^S:A67X S:E484K")
    run_cli(f"match --db {db_path} --profile S:A67G --profile S:N501Y --debug")
    run_cli(f"match --db {db_path} --profile S:A67G --DEMIS_ID 10013 --debug")
    run_cli(
        f"match --db {db_path} --DATE_DRAW 2021-03-01:2022-03-15 -o {tmp_path}/temp.tsv"
    )
    run_cli(f"match --db {db_path} --LINEAGE B.1.1.7 --with-sublineage LINEAGE --count")


# the following functions, we try to extend the test cases to make
# covsonar executes all command tools and also increase test coverage.(test reliability )
# However, the test is not for assessment validity.
def test_valid_extend(tmp_path, monkeypatch):
    monkeypatch.chdir(Path(__file__).parent)

    db_path = "data/test-with-seqs.db"
    # sonar.parse_args(["--version"])
    run_cli(
        f"match --db {db_path} --LINEAGE ^B.1.1.7 --with-sublineage LINEAGE --count"
    )
    run_cli(f"match --db {db_path} --LINEAGE ^B.1.1% AY.4% --with-sublineage LINEAGE")
    run_cli(f"match --db {db_path} --format csv -o {tmp_path}/out.csv")
    run_cli(f"match --db {db_path} --format vcf -o {tmp_path}/out.vcf")
    run_cli(
        f"restore --db {db_path} --sample IMS-10025-CVDP-00960 IMS-10087-CVDP-D484F3AD-CD8F-473C-8A5E-DB5D6A710BE5 IMS-10004-CVDP-0672526C-BAEA-4FE9-A57B-941CBCC13343 IMS-10013-CVDP-69DF29F4-D7E3-4954-94F4-65C20BE7B850 IMS-10013-CVDP-37E0BD5A-03D8-42CE-95C0-7B900B714B95"
    )

    assert filecmp.cmp(f"{tmp_path}/out.csv", "data/out.csv")
    assert filecmp.cmp(f"{tmp_path}/out.vcf", "data/out.vcf")


def test_info(tmp_path, monkeypatch):
    monkeypatch.chdir(Path(__file__).parent)
    db_path_orig = Path("data/test-with-seqs.db")
    db_path = tmp_path / "import-test.db"

    shutil.copy(db_path_orig, db_path)
    # sonar.parse_args(["--version"])
    run_cli(f" info --db {db_path}")
    run_cli(f" list-prop --db {db_path}")
    run_cli(f" dev --db {db_path}")


def test_db_management(tmp_path, monkeypatch):
    monkeypatch.chdir(Path(__file__).parent)
    db_path_orig = Path("data/test-with-seqs.db")
    db_path = tmp_path / "import-test.db"

    shutil.copy(db_path_orig, db_path)
    monkeypatch.setattr("builtins.input", lambda _: "YES")
    run_cli(f"db-upgrade --db {db_path}")
    # number_inputs = StringIO('2\n3\n')
    run_cli(f"optimize  --db {db_path}")


def test_edit_sample(tmp_path, monkeypatch):
    monkeypatch.chdir(Path(__file__).parent)
    db_path_orig = Path("data/test-with-seqs.db")
    db_path = tmp_path / "import-test.db"

    shutil.copy(db_path_orig, db_path)
    run_cli(f"add-prop --db {db_path} --name AGE --dtype float --descr descr")
    run_cli(
        f"delete --db {db_path} --sample IMS-10004-CVDP-0672526C-BAEA-4FE9-A57B-941CBCC13343 IMS-10013-CVDP-37E0BD5A-03D8-42CE-95C0-7B900B714B95"
    )
    monkeypatch.setattr("builtins.input", lambda _: "YES")
    run_cli(f"delete-prop --db {db_path}  --name AGE")
