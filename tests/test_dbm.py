import logging
from pathlib import Path
import sqlite3

from Bio.Seq import Seq
import pytest
from src.covsonar.dbm import sonarDbManager


def test_add_property_that_exists(init_writeable_dbm):
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        init_writeable_dbm.add_property(
            "newprop",
            datatype="string",
            querytype="string",
            description="new prop",
            subject="sample",
        )
        init_writeable_dbm.add_property(
            "NEWPROP",
            datatype="string",
            querytype="string",
            description="new prop",
            subject="sample",
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
    with sonarDbManager(testdb, readonly=False) as dbm:
        dbm.add_reference("REF2", "my new reference", "virus X", 1, 1)
        assert dbm.get_default_reference_accession() == "REF2"
        dbm.add_reference("REF3", "my new reference", "virus X", 1, 1)
        assert dbm.get_default_reference_accession() == "REF3"
        dbm.add_reference("REF4", "my new reference", "virus X", 1)
        assert dbm.get_default_reference_accession() == "REF3"


def test_db_writeablity(testdb):
    with sonarDbManager(testdb, readonly=True) as dbm:
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            dbm.add_property(
                "YET_ANOTHER_PROP",
                "text",
                "text",
                "my new prop stores text information",
                "sample",
            )
        assert pytest_wrapped_e.type == SystemExit
        assert (
            pytest_wrapped_e.value.code
            == "error: failed to insert data into sqlite table (attempt to write a readonly database)"
        )
    with sonarDbManager(testdb, readonly=False) as dbm:
        dbm.add_property(
            "YET_ANOTHER_PROP",
            "text",
            "text",
            "my new prop stores text information",
            "sample",
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
    with sonarDbManager(db_path_orig, readonly=True) as dbm:
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
    element = init_readonly_dbm.get_element_ids(reference_accession="MN908947.3")
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
    assert len(element) in (21, 22)

    # it must return empty list
    element = init_readonly_dbm.get_elements("none", "none")
    assert len(element) == 0


def test_iter_dna_variants(init_readonly_dbm):
    """
    [!!!] currently we setup covSonar with SARS-CoV2 info. by default.
    """
    #  element_ids > 1
    row = [x for x in init_readonly_dbm.iter_dna_variants("seq01", 12, 13, 14)]
    assert len(list(row)) == 1

    # no element_ids
    row = [x for x in init_readonly_dbm.iter_dna_variants(sample_name="seq01")]
    x = list(row)
    assert len(x) in {10, 11}  # check why it is varying


def test_upgrade_db(tmpfile_name, monkeypatch, caplog):
    """ """
    monkeypatch.chdir(Path(__file__).parent)

    # create db
    sqlite3.connect(tmpfile_name)
    with pytest.raises((ValueError, RuntimeError)) as pytest_wrapped_e:
        sonarDbManager(tmpfile_name).upgrade_db(tmpfile_name)
    assert str(pytest_wrapped_e.value) == "Error: Upgrade was not completed."


def test_query_profile(init_readonly_dbm):
    """unit test"""
    # snp
    _sql_val = init_readonly_dbm.build_genomic_profile_query("A3451T")
    logging.info(_sql_val)
    logging.info(_sql_val[1])
    logging.info(type(_sql_val[1]))
    assert all(i in _sql_val[1] for i in [1, 1, 3450, 3451, "A", "T"])
    _sql_val = init_readonly_dbm.build_genomic_profile_query("S:N501Y")
    logging.info(_sql_val[1])
    assert all(i in _sql_val[1] for i in [1, "cds", "S", 500, 501, "N", "Y"])
    # del
    _sql_val = init_readonly_dbm.build_genomic_profile_query("ORF1ab:del:3001-3004")
    logging.info(_sql_val[1])
    assert all(i in _sql_val[1] for i in [1, "cds", "ORF1ab", 3000, 3003, " "])
    _sql_val = init_readonly_dbm.build_genomic_profile_query("del:11288-11296")
    logging.info(_sql_val[1])
    assert all(i in _sql_val[1] for i in [11287, 11295, " "])
    # any AA
    _sql_val = init_readonly_dbm.build_genomic_profile_query("S:K517X")
    logging.info(_sql_val[1])
    assert all(
        i in _sql_val[1]
        for i in [
            "O",
            "K",
            "P",
            "Ω",
            "Ψ",
            "U",
            "L",
            "V",
            "B",
            "M",
            "Z",
            "E",
            "Φ",
            "X",
            "π",
            "Y",
            "J",
            "T",
            "ζ",
            "-",
            "A",
        ]
    )  # we reduce some match characters
    # any NT
    _sql_val = init_readonly_dbm.build_genomic_profile_query("C417N")
    assert all(
        i in _sql_val[1]
        for i in [
            "Y",
            "H",
            "G",
            "W",
            "R",
            "N",
            "T",
            "B",
            "K",
            "D",
            "M",
            "V",
            "C",
            "S",
            "A",
        ]
    )
    # exact AA
    _sql_val = init_readonly_dbm.build_genomic_profile_query("S:K417x")
    assert all(i in _sql_val[1] for i in ["S", 416, 417, "K", "X"])
    # excat NT
    _sql_val = init_readonly_dbm.build_genomic_profile_query("T418n")
    assert all(i in _sql_val[1] for i in [418, "T", "N"])


def test_query_profile_complexcase(init_readonly_dbm):
    """unit test"""
    # NT special case
    _sql_val = init_readonly_dbm.build_genomic_profile_query("T418nN")
    assert all(
        i in _sql_val[1]
        for i in [
            417,
            418,
            "T",
            "NR",
            "NA",
            "ND",
            "NK",
            "NV",
            "NB",
            "NN",
            "NT",
            "NW",
            "NM",
            "NC",
            "NY",
            "NS",
            "NH",
            "NG",
        ]
    )  # we reduce some match characters
    # AA special case
    _sql_val = init_readonly_dbm.build_genomic_profile_query("S:K418XX")
    assert all(
        i in _sql_val[1]
        for i in ["CC", "MZ", "ZQ", "EN", "NN", "EQ", "WΦ", "XW", "ζW", "πζ", "πY"]
    )
    # Combine AA and NT
    _sql_val = init_readonly_dbm.build_genomic_profile_query("S:K418I", "T418A")
    assert all(
        i in _sql_val[1] for i in ["S", 417, 418, "K", "I", 1, 1, 417, 418, "T", "A"]
    )


def test_query_profile_failcase(init_readonly_dbm):
    """unit test"""
    # invalid deletion of NT
    with pytest.raises(ValueError) as e:
        init_readonly_dbm.build_genomic_profile_query("del:-10000")
    assert e.type == ValueError
    assert "Please check the query statement" in str(e.value.args[0])

    # invalid deletion of AA
    with pytest.raises(ValueError) as e:
        init_readonly_dbm.build_genomic_profile_query("S:dl:501-1000")
    assert e.type == ValueError
    assert "Please check the query statement" in str(e.value.args[0])
    # K500IX malform of AA
    with pytest.raises(ValueError) as e:
        init_readonly_dbm.build_genomic_profile_query("K500IX")
    assert e.type == ValueError
    assert "Please check the query statement" in str(e.value.args[0])

    # A417xT' cannot combine AA and NT query, x is not in NT code
    with pytest.raises(ValueError) as e:
        init_readonly_dbm.build_genomic_profile_query("A417xT")
    assert e.type == ValueError
    assert "Please check the query statement" in str(e.value.args[0])

    # K417X~', ~ is not in AA code
    with pytest.raises(ValueError) as e:
        init_readonly_dbm.build_genomic_profile_query("K417X~")
    assert e.type == ValueError
    assert "Please check the query statement" in str(e.value.args[0])


def test_setup_database(tmpfile_name):
    # _failcase
    with pytest.raises(SystemExit) as e:
        sonarDbManager(dbfile=tmpfile_name + "APOLLO")
    assert e.type == SystemExit
    assert "database error" in e.value.code
