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
        except:
            exp = ""
        assert aa == exp


def test_new_default_red(setup_db):
    with sonarDBManager(setup_db, readonly=False) as dbm:
        dbm.add_reference("REF1", "my new reference", "virus X", 1, 1)
        assert dbm.get_default_reference_accession() == "REF1"
        dbm.add_reference("REF2", "my new reference", "virus X", 1, 1)
        assert dbm.get_default_reference_accession() == "REF2"
        dbm.add_reference("REF3", "my new reference", "virus X", 1)
        assert dbm.get_default_reference_accession() == "REF2"
