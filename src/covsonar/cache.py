#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

import base64
from collections import defaultdict
import hashlib
from itertools import zip_longest
import os
import pickle
import re
import shutil
import sys
from tempfile import mkdtemp
from tempfile import TemporaryDirectory
from typing import Any
from typing import DefaultDict
from typing import Dict
from typing import Iterator
from typing import List
from typing import Optional
from typing import Tuple
from typing import Union

import pandas as pd
from tqdm import tqdm

from covsonar.align import sonarAligner
from covsonar.basics import sonarBasics
from covsonar.dbm import sonarDbManager
from covsonar.logging import LoggingConfigurator


# Initialize logger
LOGGER = LoggingConfigurator.get_logger()


class sonarCache:
    """A class for managing a file-based cache for data processing."""

    def __init__(
        self,
        db: Optional[str] = None,
        outdir: Optional[str] = None,
        refacc: Optional[str] = None,
        logfile: Optional[str] = None,
        allow_updates: bool = True,
        ignore_errors: bool = False,
        temp: bool = False,
        disable_progress: bool = False,
    ):
        """
        Initialize the sonarCache object.

        Args:
            db (str): The database to import.
            outdir (str): The output directory to cache data.
            refacc (str): The reference accession.
            logfile (str): The log file to use.
            allow_updates (bool): Whether to allow updates or not.
            ignore_errors (bool): Whether to skip import errors and keep on going.
            temp (bool): Whether to set cache dir temporary or not.
            disable_progress (bool): Whether to disable progress display or not.
        """

        # arg-derived information
        self.db = db
        self.allow_updates = allow_updates
        self.refacc = refacc
        self.temp = temp
        self.ignore_errors = ignore_errors
        self.disable_progress = disable_progress

        # database-derived information
        with sonarDbManager(self.db) as dbm:
            self.refmols = dbm.get_molecule_data(
                '"molecule.accession"',
                '"molecule.id"',
                '"molecule.standard"',
                '"translation.id"',
                reference_accession=self.refacc,
            )
            self.default_refmol_acc = [
                x for x in self.refmols if self.refmols[x]["molecule.standard"] == 1
            ][0]
            self.default_refmol_id = [
                x for x in self.refmols if self.refmols[x]["molecule.id"] == 1
            ][0]
            self.sources = {
                x["molecule.accession"]: dbm.get_source(x["molecule.id"])
                for x in self.refmols.values()
            }
            self.properties = dbm.properties

        # compiled regex pattern
        self._propregex = re.compile(
            r"\[(" + "|".join(self.properties.keys()) + r")=([^\[\]=]+)\]"
        )
        self._molregex = re.compile(r"\[molecule=([^\[\]=]+)\]")

        # (file) paths
        self.basedir = (
            os.path.abspath(mkdtemp(prefix=".sonarCache_"))
            if not outdir
            else os.path.abspath(outdir)
        )
        if not os.path.exists(self.basedir):
            os.makedirs(self.basedir)

        self.logfile = (
            open(os.path.join(self.basedir, logfile), "a+") if logfile else None
        )
        self.smk_config = os.path.join(self.basedir, "config.yaml")
        self.sample_dir = os.path.join(self.basedir, "samples")
        self.seq_dir = os.path.join(self.basedir, "seq")
        self.algn_dir = os.path.join(self.basedir, "algn")
        self.var_dir = os.path.join(self.basedir, "var")
        self.ref_dir = os.path.join(self.basedir, "ref")

        os.makedirs(self.basedir, exist_ok=True)
        os.makedirs(self.seq_dir, exist_ok=True)
        os.makedirs(self.ref_dir, exist_ok=True)
        os.makedirs(self.var_dir, exist_ok=True)
        os.makedirs(self.sample_dir, exist_ok=True)

        # memories for processed files
        self._samplefiles = set()
        self._samplefiles_to_profile = set()
        self._refs = set()
        self._lifts = set()
        self._cds = set()
        self._tt = set()

        # others
        self.subdir_len = 2

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        if os.path.isdir(self.basedir) and self.temp:
            shutil.rmtree(self.basedir)
        if self.logfile:
            self.logfile.close()

    @staticmethod
    def slugify(string: str) -> str:
        """
        Slugify a string.

        Args:
            string (str): The input string.

        Returns:
            str: The slugified string.
        """
        return (
            base64.urlsafe_b64encode(string.encode("UTF-8")).decode("UTF-8").rstrip("=")
        )

    @staticmethod
    def write_pickle(fname: str, data: Any) -> None:
        """
        Write data to a pickle file.

        Args:
            fname (str): The file name.
            data (Any): The data to be pickled.
        """
        with open(fname, "wb") as handle:
            pickle.dump(data, handle)

    @staticmethod
    def read_pickle(fname: str) -> Any:
        """
        Read data from a pickle file.

        Args:
            fname (str): The file name.

        Returns:
            any: The unpickled data.
        """
        with open(fname, "rb") as handle:
            return pickle.load(handle, encoding="bytes")

    def file_collision(self, fname: str, data: Any) -> bool:
        """
        Check if the contents of a file collide with the provided data.

        Args:
            fname (str): The file name.
            data (any): The data to compare with the file contents.

        Returns:
            bool: True if there is a collision, False otherwise.
        """
        with open(fname, "r") as handle:
            if handle.read() != data:
                return True
        return False

    def get_seq_fname(self, seqhash: str) -> str:
        """
        Get the file name for a sequence file.

        Args:
            seqhash (str): The sequence hash.

        Returns:
            str: The file name.
        """
        fn = self.slugify(seqhash)
        return os.path.join(self.seq_dir, fn[: self.subdir_len], fn + ".seq")

    def get_ref_fname(self, refid: int) -> str:
        """
        Get the file name for a reference file.

        Args:
            refid (int): The reference ID.

        Returns:
            str: The file name.
        """
        return os.path.join(self.ref_dir, str(refid) + ".seq")

    def get_lift_fname(self, refid: int) -> str:
        """
        Get the file name for a lift file.

        Args:
            refid (int): The reference ID.

        Returns:
            str: The file name.
        """
        return os.path.join(self.ref_dir, str(refid) + ".lift")

    def get_cds_fname(self, refid: int) -> str:
        """
        Get the file name for a coding sequence (CDS) file.

        Args:
            refid (int): The reference ID.

        Returns:
            str: The file name.
        """
        return os.path.join(self.ref_dir, str(refid) + ".cds")

    def get_tt_fname(self, refid: int) -> str:
        """
        Get the file name for a translation table file.

        Args:
            refid (int): The reference ID.

        Returns:
            str: The file name.
        """
        return os.path.join(self.ref_dir, str(refid) + ".tt")

    def get_algn_fname(self, seqhash: str) -> str:
        """
        Get the file name for an alignment file.

        Args:
            seqhash (str): The sequence hash.

        Returns:
            str: The file name.
        """
        fn = self.slugify(seqhash)
        return os.path.join(self.algn_dir, fn[: self.subdir_len], fn + ".algn")

    def get_var_fname(self, seqhash: str) -> str:
        """
        Get the file name for a variant file.

        Args:
            seqhash (str): The sequence hash.

        Returns:
            str: The file name.
        """
        fn = self.slugify(seqhash)
        return os.path.join(self.var_dir, fn[: self.subdir_len], fn + ".var")

    def get_sample_fname(self, sample_name: str) -> str:
        """
        Get the file name for a sample file.

        Args:
            sample_name (str): The sample name.

        Returns:
            str: The file name.
        """
        fn = self.slugify(hashlib.sha1(sample_name.encode("utf-8")).hexdigest())
        return os.path.join(self.sample_dir, fn[: self.subdir_len], fn + ".sample")

    def cache_sample(
        self,
        name: str,
        seqhash: str,
        header: str,
        refmol: str,
        refmolid: int,
        refseq_id: int,
        sourceid: int,
        translation_id: int,
        sampleid: Optional[int] = None,
        algnid: Optional[int] = None,
        seqfile: Optional[str] = None,
        reffile: Optional[str] = None,
        ttfile: Optional[str] = None,
        algnfile: Optional[str] = None,
        varfile: Optional[str] = None,
        liftfile: Optional[str] = None,
        cdsfile: Optional[str] = None,
        properties: Dict[str, str] = None,
    ) -> str:
        """
        Cache a sample by saving its information to a file.

        Args:
            name (str): The sample name.
            seqhash (str): The sequence hash.
            header (str): The sample header.
            refmol (str): The reference molecule.
            refmolid (int): The reference molecule ID.
            refseq_id: The reference sequence ID.
            sourceid (int): The source ID.
            translation_id (int): The translation ID.
            sampleid (int): The sample ID, optional.
            algnid (int): The alignment ID, optional.
            seqfile (str): The sequence file path, optional
            reffile (str): The reference file path, optional
            ttfile (str): The translation table file path, optional
            algnfile (str): The alignment file path, optional
            varfile (str): The variant file path, optional
            liftfile (str): The lift file path, optional
            cdsfile (str): The CDS file path, optional
            properties (dict): The sample properties, optional.

        Returns:
            str: The file name where the sample is cached.
        """
        data = {
            "name": name,
            "sampleid": sampleid,
            "refmol": refmol,
            "refmolid": refmolid,
            "refseq_id": refseq_id,
            "sourceid": sourceid,
            "translationid": translation_id,
            "algnid": algnid,
            "header": header,
            "seqhash": seqhash,
            "seq_file": seqfile,
            "ref_file": reffile,
            "tt_file": ttfile,
            "algn_file": algnfile,
            "var_file": varfile,
            "lift_file": liftfile,
            "cds_file": cdsfile,
            "properties": properties,
        }

        fname = self.get_sample_fname(name)  # get full file path
        os.makedirs(os.path.dirname(fname), exist_ok=True)
        self.write_pickle(fname, data)

        self._samplefiles.add(fname)

        if algnid is None:
            self._samplefiles_to_profile.add(fname)  # select sample file for alignment

        return fname

    def iter_samples(self) -> Iterator[Dict[str, Any]]:
        """
        Iterate over all cached samples.

        Yields:
            dict: A dictionary containing the information for a single sample.
        """
        for fname in self._samplefiles:
            yield self.read_pickle(fname)

    def cache_sequence(self, seqhash: str, sequence: str) -> str:
        """
        Cache a sequence by saving it to a file.

        Args:
            seqhash (str): The sequence hash.
            sequence (str): The sequence.

        Returns:
            str: The file name where the sequence is cached.
        """
        fname = self.get_seq_fname(seqhash)
        if os.path.isfile(fname) and self.file_collision(fname, sequence):
            LOGGER.critical(
                "Seqhash collision: sequences differ for seqhash " + seqhash + "."
            )
            sys.exit(1)

        os.makedirs(os.path.dirname(fname), exist_ok=True)
        with open(fname, "w") as handle:
            handle.write(sequence)
        return fname

    def cache_reference(self, refid: int, sequence: str) -> str:
        """
        Cache a reference by saving it to a file.

        Args:
            refid (int): The reference ID.
            sequence (str): The reference sequence.

        Returns:
            str: The file name where the reference is cached.
        """
        fname = self.get_ref_fname(refid)
        if refid not in self._refs:
            with open(fname, "w") as handle:
                handle.write(sequence)
            self._refs.add(refid)
        return fname

    def cache_translation_table(self, translation_id: int, dbm: sonarDbManager) -> str:
        """
        Cache a translation table by saving it to a file.

        Args:
            translation_id (int): The translation table ID.
            dbm (sonarDbManager): The database manager.

        Returns:
            str: The file name where the translation table is cached.
        """
        fname = self.get_tt_fname(translation_id)
        if translation_id not in self._tt:
            self.write_pickle(fname, dbm.get_translation_dict(translation_id))
            self._tt.add(translation_id)
        return fname

    @staticmethod
    def get_cds_coord_list(cds: dict) -> List[int]:
        """
        Return a ordered list of all genomic positions for the repsctive CDS.

        Args:
            cds (dict): The cds dictionary.

        Returns:
            List[int]:
                list of integers representig the genomic coding positions in order.
        """
        coords = []
        for data in cds["ranges"]:
            if data[-1] == 1:
                coords.extend(list(range(data[0], data[1])))
            else:
                coords.extend(list(range(data[0], data[1]))[::-1])
        return coords

    def cache_cds(self, refid: int, refmol_acc: str) -> str:
        """
        Cache coding sequences (CDS) by saving them to a file.

        Args:
            refid (int): The reference ID.
            refmol_acc (str): The reference molecule accession number.

        Returns:
            str: The file name where the CDS are cached.
        """
        fname = self.get_cds_fname(refid)
        if refmol_acc not in self._cds:
            cols = ["elemid", "pos", "end"]
            rows = []
            for cds in self.iter_cds(refmol_acc):
                elemid = cds["id"]
                for pos in self.get_cds_coord_list(cds):
                    rows.append([elemid, pos, 0])
                rows[-1][2] = 1
            df = pd.DataFrame(rows, columns=cols)
            df.to_pickle(fname)
            self._cds.add(refmol_acc)
        return fname

    def cache_lift(self, refid: int, refmol_acc: str, sequence: str) -> str:
        """
        Cache lift information by saving it to a file.

        Args:
            refid (int): The reference ID.
            refmol_acc (str): The reference molecule accession number.
            sequence(str): The reference sequence.

        Returns:
            str: The file name where the lift information is cached.
        """
        fname = self.get_lift_fname(refid)
        if refmol_acc not in self._lifts:
            cols = [
                "elemid",
                "nucPos1",
                "nucPos2",
                "nucPos3",
                "ref1",
                "ref2",
                "ref3",
                "alt1",
                "alt2",
                "alt3",
                "symbol",
                "aaPos",
                "aa",
            ]
            rows = []
            for cds in self.iter_cds(refmol_acc):
                elemid = cds["id"]
                symbol = cds["symbol"]
                coords = [
                    list(group)
                    for group in zip_longest(
                        *[iter(self.get_cds_coord_list(cds))] * 3, fillvalue="-"
                    )
                ]  # provide a list of triplets for all CDS coordinates
                seq = list(reversed(cds["sequence"] + "*"))
                for aa_pos, nuc_pos_list in enumerate(coords):
                    rows.append(
                        [elemid]
                        + nuc_pos_list
                        + [
                            sequence[nuc_pos_list[0]],
                            sequence[nuc_pos_list[1]],
                            sequence[nuc_pos_list[2]],
                        ]
                        * 2
                        + [symbol, aa_pos, seq.pop()]
                    )
            df = pd.DataFrame(rows, columns=cols)
            df.to_pickle(fname)
            self._lifts.add(refmol_acc)
        return fname

    def process_fasta_entry(self, header: str, seq: str) -> Dict[str, Union[str, int]]:
        """
        Process a single entry from a FASTA file and extract the relevant information.

        Args:
            header (str): The header of the FASTA entry.
            seq (str): The sequence of the FASTA entry.

        Returns:
            dict: A dictionary containing the sample information.
        """
        sample_id = header.replace("\t", " ").replace("|", " ").split(" ")[0]
        refmol = self.get_refmol(header)
        if not refmol:
            LOGGER.error(
                sample_id
                + " refers to an unknown reference molecule ("
                + self._molregex.search(header)
                + ")."
            )
            sys.exit(1)
        seq = sonarBasics.harmonize_seq(seq)
        seqhash = sonarBasics.hash_seq(seq)
        refmolid = self.refmols[refmol]["molecule.id"]
        return {
            "name": sample_id,
            "header": header,
            "seqhash": seqhash,
            "sequence": seq,
            "refmol": refmol,
            "refmolid": refmolid,
            "translation_id": self.refmols[refmol]["translation.id"],
            "properties": self.get_properties(header),
        }

    def iter_fasta(self, *fnames: str) -> Iterator[Dict[str, Union[str, int]]]:
        """
        Iterate over FASTA files and extract sample information from each entry.

        Args:
            fnames (str): The paths of the FASTA files.

        Yields:
            dict: A dictionary containing the information for a single sample.
        """
        for fname in fnames:
            with sonarBasics.open_file_autodetect(fname) as handle, tqdm(
                desc="processing " + fname + "...",
                total=os.path.getsize(fname),
                unit="bytes",
                unit_scale=True,
                bar_format="{desc} {percentage:3.0f}% [{n_fmt}/{total_fmt}, {elapsed}<{remaining}, {rate_fmt}{postfix}]",
                disable=self.disable_progress,
            ) as pbar:
                seq = []
                header = None
                for line in handle:
                    pbar.update(len(line))
                    line = line.strip()
                    if line.startswith(">"):
                        if seq:
                            yield self.process_fasta_entry(header, "".join(seq))
                            seq = []
                        header = line[1:]
                    else:
                        seq.append(line)
                if seq:
                    yield self.process_fasta_entry(header, "".join(seq))

    def get_refmol(self, fasta_header: str) -> Optional[str]:
        """
        Get the reference molecule accession number from the FASTA header.

        Args:
            fasta_header (str): The header of the FASTA entry.

        Returns:
            str/None: The molecule accession number if found in the database, or None if not found int he dabase, or the default molecule accession number if no molecule accession found in the FASTA.
        """
        mol = self._molregex.search(fasta_header)
        if mol:
            try:
                return self.refmols[mol]["accession"]
            except Exception:
                return None
        return self.default_refmol_acc

    def get_refseq(self, refmol_acc: str) -> Optional[str]:
        """
        Get the reference sequence for a given reference molecule accession number.

        Args:
            refmol_acc (str): The reference molecule accession number.

        Returns:
            str/None: The reference sequence if found in the database, or None if not found.
        """
        try:
            return self.sources[refmol_acc]["sequence"]
        except Exception:
            return None

    def iter_cds(
        self, refmol_acc: str
    ) -> Iterator[Dict[str, Union[int, str, List[Tuple[int, int, int]]]]]:
        """
        Iterate over the coding sequences (CDS) for a given reference molecule.

        Args:
            refmol_acc (str): The reference molecule accession number.

        Yields:
            dict: A dictionary containing the information for a single coding sequence.
        """
        cds = {}
        prev_elem = None
        with sonarDbManager(self.db) as dbm:
            for row in dbm.get_annotation(
                reference_accession=refmol_acc,
                molecule_accession=refmol_acc,
                element_type="cds",
            ):
                if prev_elem is None:
                    prev_elem = row["element.id"]
                elif row["element.id"] != prev_elem:
                    yield cds
                    cds = {}
                    prev_elem = row["element.id"]
                if not cds:
                    cds = {
                        "id": row["element.id"],
                        "accession": row["element.accession"],
                        "symbol": row["element.symbol"],
                        "sequence": row["element.sequence"],
                        "ranges": [
                            (
                                row["elempart.start"],
                                row["elempart.end"],
                                row["elempart.strand"],
                            )
                        ],
                    }
                else:
                    cds["ranges"].append(
                        (
                            row["elempart.start"],
                            row["elempart.end"],
                            row["elempart.strand"],
                        )
                    )
        if cds:
            yield cds

    def get_refseq_id(self, refmol_acc: str) -> Optional[int]:
        """
        Get the reference sequence ID for a given reference molecule accession number.

        Args:
            refmol_acc (str): The reference molecule accession number.

        Returns:
            int/None: The reference sequence ID if found, or None if not found.
        """
        try:
            return self.sources[refmol_acc]["id"]
        except Exception:
            return None

    def get_refhash(self, refmol_acc: str) -> Optional[str]:
        """
        Get the sequence hash for a given reference molecule accession number.

        Args:
            refmol_acc (str): The reference molecule accession number.

        Returns:
            str/None: The sequence hash if found, or None if not found.
        """
        try:
            if "seqhash" not in self.sources[refmol_acc]:
                self.sources[refmol_acc]["seqhash"] = sonarBasics.hash_seq(
                    self.sources[refmol_acc]["sequence"]
                )
            return self.sources[refmol_acc]["seqhash"]
        except Exception:
            return None

    def get_properties(self, fasta_header: str) -> Dict[str, str]:
        """
        Get the properties from the FASTA header.

        Args:
            fasta_header (str): The header of the FASTA entry.

        Returns:
            dict: A dictionary of properties extracted from the header.
        """
        return {x.group(1): x.group(2) for x in self._propregex.finditer(fasta_header)}

    def add_fasta(
        self,
        *fnames: str,
        properties: DefaultDict[str, Dict[str, str]] = defaultdict(dict),
    ) -> None:
        """
        Add entries from a given FASTA file to the cache.

        Args:
            *fnames (str): Variable-length argument list of FASTA file names.
            properties (defaultdict): Additional properties to associate with the samples.

        Raises:
            OSError: If there is an error while reading or writing the cache files.
            ValueError: If the sample refers to an unknown reference molecule.
        """
        default_properties = {
            x: self.properties[x]["standard"] for x in self.properties
        }
        with sonarDbManager(self.db, readonly=False) as dbm:
            for fname in fnames:
                for data in self.iter_fasta(fname):
                    # check sample
                    data["sampleid"], seqhash = dbm.get_sample_data(data["name"])
                    data["sourceid"] = dbm.get_source(data["refmolid"])["id"]

                    # check properties
                    if data["sampleid"] is None:
                        props = default_properties.copy()
                        props.update(data["properties"])
                        props.update(properties[data["sampleid"]])
                        data["properties"] = props
                    elif not self.allow_updates:
                        continue
                    else:
                        data["properties"].update(properties[data["sampleid"]])

                    # check reference
                    data["refseq_id"] = self.get_refseq_id(data["refmol"])
                    self.write_checkref_log(data, data["refseq_id"])

                    # check alignment
                    data["algnid"] = dbm.get_alignment_id(
                        data["seqhash"], data["sourceid"]
                    )
                    data = self.add_data_files(
                        data, seqhash, dbm
                    )  # seqhash is provided for validation
                    del data["sequence"]
                    self.cache_sample(**data)

    def write_checkref_log(
        self, data: Dict[str, Any], refseq_id: Optional[int] = None
    ) -> None:
        """
        Write a log entry for a reference sequence check.

        Args:
            data (dict): The data dict of the sample.
            refseq_id (int): The reference sequence ID, optional.

        Raises:
            ValueError: If the reference sequence is unknown and `ignore_errors` is False (via self.log function).
        """
        if not refseq_id:
            if not self.ignore_errors:
                LOGGER.error(
                    "Fasta header refers to an unknown reference ('"
                    + data["header"]
                    + "').",
                )
                sys.exit(1)
            else:
                LOGGER.warning(
                    "Skipping '"
                    + data["name"]
                    + "' referring to an unknown reference ('"
                    + data["header"]
                    + "')."
                )

    def add_data_files(
        self, data: Dict[str, Any], seqhash: str, dbm: sonarDbManager
    ) -> Dict[str, Any]:
        """
        Assign additional data files to a given data dict of a sample.

        Args:
            data (dict): The data for the sample.
            seqhash (str): The sequence hash, optional.
            dbm (sonarDbManager): The sonarDbManager instance.

        Returns:
            (dict) The updated data dictionary.
        """
        if data["algnid"] is None:
            data["seqfile"] = self.cache_sequence(data["seqhash"], data["sequence"])
            data["reffile"] = self.cache_reference(
                data["refseq_id"], self.get_refseq(data["refmol"])
            )
            data["ttfile"] = self.cache_translation_table(data["translation_id"], dbm)
            data["liftfile"] = self.cache_lift(
                data["refseq_id"], data["refmol"], self.get_refseq(data["refmol"])
            )
            data["cdsfile"] = self.cache_cds(data["refseq_id"], data["refmol"])
            data["algnfile"] = self.get_algn_fname(
                data["seqhash"] + "@" + self.get_refhash(data["refmol"])
            )
            data["varfile"] = self.get_var_fname(
                data["seqhash"] + "@" + self.get_refhash(data["refmol"])
            )
        else:
            if data["seqhash"] != seqhash:
                data["seqfile"] = self.cache_sequence(data["seqhash"], data["sequence"])
            else:
                data["seqhash"] = None
                data["seqfile"] = None
            data["reffile"] = None
            data["ttfile"] = None
            data["liftfile"] = None
            data["cdsfile"] = None
            data["algnfile"] = None
            data["varfile"] = None

        return data

    def import_cached_samples(self) -> None:
        """
        Import the cached samples into the database.

        Raises:
            OSError: If there is an error while reading the cache files.
            ValueError: If the original sequence of a sample cannot be restored from the stored genomic profile.
        """
        refseqs = {}
        with sonarDbManager(self.db, readonly=False) as dbm:
            for sample_data in tqdm(
                self.iter_samples(),
                total=len(self._samplefiles),
                desc="importing samples...",
                unit="samples",
                bar_format="{desc} {percentage:3.0f}% [{n_fmt}/{total_fmt}, {elapsed}<{remaining}, {rate_fmt}{postfix}]",
                disable=self.disable_progress,
            ):
                try:
                    # nucleotide level import
                    if sample_data["seqhash"] is not None:
                        dbm.insert_sample(sample_data["name"], sample_data["seqhash"])
                        algnid = dbm.insert_alignment(
                            sample_data["seqhash"], sample_data["sourceid"]
                        )

                    if sample_data["var_file"] is not None:
                        with open(sample_data["var_file"], "r") as handle:
                            for line in handle:
                                if line == "//":
                                    break
                                vardat = line.strip("\r\n").split("\t")
                                dbm.insert_variant(
                                    algnid,
                                    element_id=vardat[0],
                                    ref=vardat[1],
                                    alt=vardat[2],
                                    start=vardat[3],
                                    end=vardat[4],
                                    frameshift=vardat[5],
                                    label=vardat[6],
                                )
                            if line != "//":
                                raise OSError(
                                    "cache error: corrupted file ("
                                    + sample_data["var_file"]
                                    + ")"
                                )

                    # paranoia test
                    if sample_data["seqhash"] is not None:
                        self.paranoid_test(refseqs, sample_data, dbm)

                except Exception as e:
                    LOGGER.exception(
                        "Unknown import error occurred. Debugging Information: %s", e
                    )
                    LOGGER.critical(
                        "Exiting the program due to an unrecoverable import error."
                    )
                    sys.exit(1)

    def paranoid_test(
        self, refseqs: Dict[str, str], sample_data: Dict[str, Any], dbm: sonarDbManager
    ) -> None:
        """
        Perform a test to verify the integrity of the imported sample.

        Args:
            refseqs (dict): A dictionary to cache reference sequences.
            sample_data (dict): The data for the sample.
            dbm (sonarDbManager): The sonarDbManager session to use.

        Raises:
            SystemExit: If the original sequence of a sample cannot be restored from the stored genomic profile.

        Returns:
            True
        """
        # Consider: change from refmolid (moleculeID) to sourceid
        # to support multiple references.
        try:
            seq = list(refseqs[sample_data["sourceid"]])
        except Exception:
            refseqs[sample_data["sourceid"]] = list(
                dbm.get_sequence(sample_data["sourceid"])
            )
            seq = list(refseqs[sample_data["sourceid"]])

        # sequence recovery based on nucleotide variations
        prefix = ""
        gaps = {".", " "}

        for vardata in dbm.iter_dna_variants(
            sample_data["name"], sample_data["sourceid"]
        ):
            # insert deletions
            if vardata["variant.alt"] in gaps:
                for i in range(vardata["variant.start"], vardata["variant.end"]):
                    seq[i] = ""
            # insert snp or insertion
            elif vardata["variant.start"] >= 0:
                if len(vardata["variant.alt"]) == 1:
                    seq[vardata["variant.start"]] = (
                        vardata["variant.alt"] + seq[vardata["variant.start"]][1:]
                    )
                else:
                    seq[vardata["variant.start"]] += vardata["variant.alt"][1:]
            # insert prefix (extension before start)
            else:
                prefix = vardata["variant.alt"]
        # assemble stored sequence
        seq = prefix + "".join(seq)
        # read original input sequence
        with open(sample_data["seq_file"], "r") as handle:
            orig_seq = handle.read()
        # compare
        if seq != orig_seq:
            aligner = sonarAligner()
            with TemporaryDirectory() as tempdir:
                qryfile = os.path.join(tempdir, "qry")
                reffile = os.path.join(tempdir, "ref")
                with open(qryfile, "w") as handle:
                    handle.write(seq)
                with open(reffile, "w") as handle:
                    handle.write(orig_seq)
                ref, qry, cigar = aligner.align(seq, orig_seq)
            with open("paranoid.alignment.fna", "w") as handle:
                handle.write(
                    ">original_"
                    + sample_data["name"]
                    + "\n"
                    + ref
                    + "\n>restored_"
                    + sample_data["name"]
                    + "\n"
                    + qry
                    + "\n"
                )
            LOGGER.critical(
                "Original sequence of sample "
                + sample_data["name"]
                + " cannot be restored from stored genomic profile for sample (see paranoid.alignment.fna):"
                + " "
                + cigar
            )
            LOGGER.error("Paranoid test failed.")
            sys.exit(1)
        return True


if __name__ == "__main__":
    pass
