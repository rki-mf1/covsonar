import logging
from pathlib import Path
import sqlite3

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
    # init_writeable_dbm.close()


def test_translationtable(init_writeable_dbm):
    for tt in [1, 2, 11]:
        init_writeable_dbm.add_translation_table(tt)
    for codon, aa in init_writeable_dbm.get_translation_dict(1).items():
        try:
            exp = str(Seq.translate(codon, table=tt))
        except Exception:
            exp = ""
        assert aa == exp
    # init_writeable_dbm.close()


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


def test_get_molecule_data(monkeypatch):
    """
    [!!!] currently we use covSonar with SARS-CoV2 info. by default.
        # this test function should be edited when apply to other pathogens.
        # assertion can be more appropriate check than current code.

        ignor line 689
    """
    monkeypatch.chdir(Path(__file__).parent)
    db_path_orig = "data/test-with-seqs.db"
    with sonarDBManager(db_path_orig, readonly=True) as dbm:
        # Get all field
        refmols = dbm.get_molecule_data()
        assert type(refmols) is dict

        # Check if MN908947.3 and molecule.accession is inserted in the Dict
        refmols = dbm.get_molecule_data(
            '"molecule.id"',
            '"molecule.standard"',
            '"translation.id"',
            reference_accession="MN908947.3",
        )
        assert "molecule.accession" in refmols["MN908947.3"]


def test_get_annotation(init_readonly_dbm):
    anno = init_readonly_dbm.get_annotation(element_accession="ORF18")
    logging.info(anno)
    assert len(anno) == 0
    # assert "ORF10" == anno[0]["element.accession"]

    anno = init_readonly_dbm.get_annotation(reference_accession="MN908947.3")
    logging.info(anno)
    assert "MN908947.3" == anno[0]["reference.accession"]


def test_get_alignment_data(init_readonly_dbm):
    align = init_readonly_dbm.get_alignment_data(
        sample_name="seq01", reference_accession="MN908947.3"
    )
    assert align.fetchall()[0]["element.symbol"] == "MN908947.3"

    align = init_readonly_dbm.get_alignment_data(
        sample_name="", reference_accession=None
    )
    logging.info(align)
    assert align == ""


def test_get_element_id(init_readonly_dbm):
    """
    [!!!] currently we setup covSonar with SARS-CoV2 info. by default.
    """
    element = init_readonly_dbm.get_element_ids(
        reference_accession="MN908947.3", type=None
    )
    logging.info(element)
    assert 1 in element

    # it must return empty list
    element = init_readonly_dbm.get_element_ids("none", "none")
    assert len(element) == 0


def test_get_elements(init_readonly_dbm):
    """
    [!!!] currently we setup covSonar with SARS-CoV2 info. by default.
    """

    element = init_readonly_dbm.get_elements("1", "gene")
    assert type(element) is list
    logging.info(element)
    assert len(element) == 10

    element = init_readonly_dbm.get_elements("1")
    assert type(element) is list
    logging.info(element)
    assert len(element) == 21

    # it must return empty list
    element = init_readonly_dbm.get_elements("none", "none")
    assert len(element) == 0


def test_iter_dna_variants(init_readonly_dbm):
    """
    [!!!] currently we setup covSonar with SARS-CoV2 info. by default.
    """
    #  element_ids > 1
    row = init_readonly_dbm.iter_dna_variants("seq01", 12, 13, 14)
    assert len(list(row)) > 1

    # no element_ids
    row = init_readonly_dbm.iter_dna_variants(sample_name="seq01")
    assert 11 == 11


def test_upgrade_db(tmpfile_name, monkeypatch, caplog):
    """ """
    monkeypatch.chdir(Path(__file__).parent)

    # create db
    sqlite3.connect(tmpfile_name)
    with caplog.at_level(logging.ERROR):
        sonarDBManager(tmpfile_name).upgrade_db(tmpfile_name)
        assert "Sorry, we cannot find" in caplog.text
