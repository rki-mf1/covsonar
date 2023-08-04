import logging
from pathlib import Path
import sqlite3

from Bio.Seq import Seq
import pytest
from src.covsonar.dbm import sonarDbManager


def test_add_property_that_exists(init_writeable_dbm, logger, caplog):
    with caplog.at_level(logging.ERROR, logger=logger.name), pytest.raises(
        SystemExit
    ) as pytest_wrapped_e:
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
        assert (
            "A property named 'NEWPROP' already exists in the given database.."
            == caplog.records[-1].message
        )
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 1


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


def test_db_writeablity(testdb, logger, caplog):
    with sonarDbManager(testdb, readonly=True) as dbm, caplog.at_level(
        logging.ERROR, logger=logger.name
    ), pytest.raises(SystemExit) as pytest_wrapped_e:
        dbm.add_property(
            "YET_ANOTHER_PROP",
            "text",
            "text",
            "my new prop stores text information",
            "sample",
        )
        assert (
            "Failed to insert data into sqlite table (attempt to write a readonly database)."
            == caplog.records[-1].message
        )
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 1
    assert pytest_wrapped_e.type == SystemExit
    assert (
        caplog.records[-1].message
        == "Failed to insert data into sqlite table (attempt to write a readonly database)."
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
    assert len(anno) == 0
    # assert "ORF10" == anno[0]["element.accession"]

    anno = init_readonly_dbm.get_annotation(reference_accession="MN908947.3")
    assert "MN908947.3" == anno[0]["reference.accession"]


def test_get_alignment_data(init_readonly_dbm):
    align = init_readonly_dbm.get_alignment_data(
        sample_name="seq01", reference_accession="MN908947.3"
    )
    assert align.fetchall()[0]["element.symbol"] == "MN908947.3"

    align = init_readonly_dbm.get_alignment_data(
        sample_name="", reference_accession=None
    )
    assert align == ""


def test_get_element_id(init_readonly_dbm):
    """
    [!!!] currently we setup covSonar with SARS-CoV2 info. by default.
    """
    element = init_readonly_dbm.get_element_ids(reference_accession="MN908947.3")
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
    assert len(element) == 10

    element = init_readonly_dbm.get_elements("1")
    assert type(element) is list
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


def test_upgrade_db(tmpfile_name, monkeypatch, logger, caplog):
    """ """
    monkeypatch.chdir(Path(__file__).parent)

    # create db
    sqlite3.connect(tmpfile_name)
    with pytest.raises(SystemExit) as pytest_wrapped_e, caplog.at_level(
        logging.ERROR, logger=logger.name
    ):
        sonarDbManager(tmpfile_name).upgrade_db(tmpfile_name)
        assert (
            "Sorry, but automated migration does not support databases of version 0."
            == caplog.records[-1].message
        )
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 1


def test_query_profile(init_readonly_dbm):
    """unit test"""
    # snp
    _sql_val = init_readonly_dbm.create_genomic_profile_sql("A3451T")
    assert all(i in _sql_val[1] for i in [1, 1, 3450, 3451, "A", "T"])
    _sql_val = init_readonly_dbm.create_genomic_profile_sql("S:N501Y")
    assert all(i in _sql_val[1] for i in [1, "cds", "S", 500, 501, "N", "Y"])
    # del
    _sql_val = init_readonly_dbm.create_genomic_profile_sql("ORF1ab:del:3001-3004")
    assert all(i in _sql_val[1] for i in [1, "cds", "ORF1ab", 3000, 3003, " "])
    _sql_val = init_readonly_dbm.create_genomic_profile_sql("del:11288-11296")
    assert all(i in _sql_val[1] for i in [11287, 11295, " "])
    # any AA
    _sql_val = init_readonly_dbm.create_genomic_profile_sql("S:K517X")
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
    _sql_val = init_readonly_dbm.create_genomic_profile_sql("C417N")
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
    _sql_val = init_readonly_dbm.create_genomic_profile_sql("S:K417x")
    assert all(i in _sql_val[1] for i in ["S", 416, 417, "K", "X"])
    # excat NT
    _sql_val = init_readonly_dbm.create_genomic_profile_sql("T418n")
    assert all(i in _sql_val[1] for i in [418, "T", "N"])


def test_query_profile_complexcase(init_readonly_dbm):
    """unit test"""
    # NT special case
    _sql_val = init_readonly_dbm.create_genomic_profile_sql("T418nN")
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
    _sql_val = init_readonly_dbm.create_genomic_profile_sql("S:K418XX")
    assert all(
        i in _sql_val[1]
        for i in ["CC", "MZ", "ZQ", "EN", "NN", "EQ", "WΦ", "XW", "ζW", "πζ", "πY"]
    )
    # Combine AA and NT
    _sql_val = init_readonly_dbm.create_genomic_profile_sql("S:K418I", "T418A")
    assert all(
        i in _sql_val[1] for i in ["S", 417, 418, "K", "I", 1, 1, 417, 418, "T", "A"]
    )


def test_query_profile_failcase(init_readonly_dbm, logger, caplog):
    """unit test"""
    # invalid deletion of NT
    with caplog.at_level(logging.ERROR, logger=logger.name), pytest.raises(
        SystemExit
    ) as pytest_wrapped_e:
        init_readonly_dbm.create_genomic_profile_sql("del:-10000")
        assert (
            "The alternate allele notation 'del:-10000' is invalid."
            == caplog.records[-1].message
        )
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 1

    # invalid deletion of AA
    with caplog.at_level(logging.ERROR, logger=logger.name), pytest.raises(
        SystemExit
    ) as pytest_wrapped_e:
        init_readonly_dbm.create_genomic_profile_sql("S:dl:501-1000")
        assert "Please check the query statement." == caplog.records[-1].message
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 1

    # K500IX malform of AA
    with caplog.at_level(logging.ERROR, logger=logger.name), pytest.raises(
        SystemExit
    ) as pytest_wrapped_e:
        init_readonly_dbm.create_genomic_profile_sql("K500IX")
        assert (
            "The alternate allele notation 'IX' is invalid."
            == caplog.records[-1].message
        )
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 1

    # A417xT' cannot combine AA and NT query, x is not in NT code
    with caplog.at_level(logging.ERROR, logger=logger.name), pytest.raises(
        SystemExit
    ) as pytest_wrapped_e:
        init_readonly_dbm.create_genomic_profile_sql("A417xT")
        assert (
            "The alternate allele notation 'xT' is invalid."
            == caplog.records[-1].message
        )
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 1

    # K417X~', ~ is not in AA code
    with caplog.at_level(logging.ERROR, logger=logger.name), pytest.raises(
        SystemExit
    ) as pytest_wrapped_e:
        init_readonly_dbm.create_genomic_profile_sql("K417X~")
        assert (
            "The alternate allele notation 'X~' is invalid."
            == caplog.records[-1].message
        )
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 1


def test_setup_database(tmpfile_name, logger, caplog):
    # _failcase
    with caplog.at_level(logging.ERROR, logger=logger.name), pytest.raises(
        SystemExit
    ) as pytest_wrapped_e:
        sonarDbManager(dbfile=tmpfile_name + "APOLLO")
        assert "Please check the query statement." == caplog.records[-1].message
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 1
