from Bio.Seq import Seq
import pytest

from covsonar.dbm import sonarDBManager


def test_add_property_that_exists(init_writeable_dbm):
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        init_writeable_dbm.add_property(
            "newprop", datatype="string", querytype="string", description="new prop"
        )
        init_writeable_dbm.add_property(
            "NEWPROP", datatype="string", querytype="string", description="new prop"
        )
    assert pytest_wrapped_e.type == SystemExit
    assert (
        pytest_wrapped_e.value.code
        == "error: a property named NEWPROP already exists in the given database."
    )


def test_translationtable(init_writeable_dbm):
    for tt in [1, 2, 11]:
        init_writeable_dbm.add_translation_table(tt)
    for codon, aa in init_writeable_dbm.get_translation_dict(1).items():
        try:
            exp = str(Seq.translate(codon, table=tt))
        except Exception:
            exp = ""
        assert aa == exp


def test_add_reference(testdb):
    with sonarDBManager(testdb, readonly=False) as dbm:
        dbm.add_reference("REF2", "my new reference", "virus X", 1, 1)
        assert dbm.get_default_reference_accession() == "REF2"
        dbm.add_reference("REF3", "my new reference", "virus X", 1, 1)
        assert dbm.get_default_reference_accession() == "REF3"
        dbm.add_reference("REF4", "my new reference", "virus X", 1)
        assert dbm.get_default_reference_accession() == "REF3"


def test_db_writeablity(testdb):
    with sonarDBManager(testdb, readonly=True) as dbm:
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            dbm.add_property(
                "YET_ANOTHER_PROP",
                "text",
                "text",
                "my new prop stores text information",
            )
        assert pytest_wrapped_e.type == SystemExit
        assert (
            pytest_wrapped_e.value.code
            == "error: failed to insert data into sqlite table (attempt to write a readonly database)"
        )
    with sonarDBManager(testdb, readonly=False) as dbm:
        dbm.add_property(
            "YET_ANOTHER_PROP", "text", "text", "my new prop stores text information"
        )
        assert "YET_ANOTHER_PROP" in dbm.properties
        dbm.delete_property("YET_ANOTHER_PROP")
        assert "YET_ANOTHER_PROP" not in dbm.properties
