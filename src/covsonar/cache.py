#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

import base64
from collections import defaultdict
import gzip
import hashlib
import lzma
import os
import pickle
from pickle import load as load_pickle
import pprint
import re
import shutil
import sys
from tempfile import mkdtemp
from tempfile import TemporaryDirectory

import pandas as pd
from tqdm import tqdm
import yaml

from .align import sonarAligner
from .basics import sonarBasics
from .dbm import sonarDBManager

pp = pprint.PrettyPrinter(indent=4)


class sonarCache:
    """ """

    def __init__(
        self,
        db=None,
        outdir=None,
        refacc=None,
        logfile=None,
        allow_updates=True,
        ignore_errors=False,
        temp=False,
        debug=False,
        disable_progress=False,
    ):
        self.db = db
        self.allow_updates = allow_updates
        self.debug = debug
        self.refacc = refacc
        self.temp = temp
        self.ignore_errors = ignore_errors
        self.disable_progress = disable_progress

        with sonarDBManager(self.db, debug=self.debug) as dbm:
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
        self._propregex = re.compile(
            r"\[(" + "|".join(self.properties.keys()) + r")=([^\[\]=]+)\]"
        )
        self._molregex = re.compile(r"\[molecule=([^\[\]=]+)\]")

        self.logfile = open(logfile, "w") if logfile else None
        self.basedir = (
            os.path.abspath(mkdtemp(prefix=".sonarCache_"))
            if not outdir
            else os.path.abspath(outdir)
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
        # os.makedirs(self.algn_dir, exist_ok=True)
        os.makedirs(self.var_dir, exist_ok=True)
        os.makedirs(self.sample_dir, exist_ok=True)
        self._samplefiles = set()
        self._samplefiles_to_profile = set()
        self._refs = set()
        self._lifts = set()
        self._tt = set()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        if os.path.isdir(self.basedir) and self.temp:
            shutil.rmtree(self.basedir)
        if self.logfile:
            self.logfile.close()
        if self.sample_yaml:
            self.sample_yaml.close()

    @staticmethod
    def slugify(string):
        return (
            base64.urlsafe_b64encode(string.encode("UTF-8")).decode("UTF-8").rstrip("=")
        )

    @staticmethod
    def deslugify(string):
        while len(string) % 3 != 0:
            string += "="
        return base64.urlsafe_b64decode(string).decode("utf-8")

    def log(self, msg, die=False, errtype="error"):
        if self.logfile:
            self.logfile.write(msg)
        elif not die:
            sys.stderr(msg)
        else:
            exit(errtype + ": " + msg)

    @staticmethod
    def write_pickle(fname, data):
        with open(fname, "wb") as handle:
            pickle.dump(data, handle)

    @staticmethod
    def read_pickle(fname):
        with open(fname, "rb") as handle:
            return pickle.load(handle, encoding="bytes")

    @staticmethod
    def pickle_collision(fname, data):
        if os.path.isfile(fname) and load_pickle(fname) != data:
            return True
        return False

    @staticmethod
    def file_collision(fname, data):
        with open(fname, "r") as handle:
            if handle.read() != data:
                return True
        return False

    def sample_collision(self, key, datadict):
        if sonarCache.read_pickle(key) != datadict:
            return True
        return False

    def write_smk_config(self):
        data = {
            "debug": self.debug,
            "sample_dir": self.sample_dir,
            "seq_dir": self.seq_dir,
            "algn_dir": self.algn_dir,
            "var_dir": self.var_dir,
        }

        with open(self.smk_config, "w") as handle:
            yaml.dump(data, handle)

    def get_seq_fname(self, seqhash):
        fn = self.slugify(seqhash)
        return os.path.join(self.seq_dir, fn[:2], fn + ".seq")

    def get_ref_fname(self, refid):
        return os.path.join(self.ref_dir, str(refid) + ".seq")

    def get_lift_fname(self, refid):
        return os.path.join(self.ref_dir, str(refid) + ".lift")

    def get_tt_fname(self, refid):
        return os.path.join(self.ref_dir, str(refid) + ".tt")

    def get_algn_fname(self, seqhash):
        fn = self.slugify(seqhash)
        return os.path.join(self.algn_dir, fn[:2], fn + ".algn")

    def get_var_fname(self, seqhash):
        fn = self.slugify(seqhash)
        return os.path.join(self.var_dir, fn[:2], fn + ".var")

    def get_sample_fname(self, sample_name):
        fn = self.slugify(hashlib.sha1(sample_name.encode("utf-8")).hexdigest())
        return os.path.join(self.sample_dir, fn[:2], fn + ".sample")

    def cache_sample(
        self,
        name,
        sampleid,
        seqhash,
        header,
        refmol,
        refmolid,
        sourceid,
        translation_id,
        algnid,
        seqfile,
        reffile,
        ttfile,
        algnfile,
        varfile,
        liftfile,
        properties,
    ):
        """
        The function takes in a bunch of arguments and returns a filename.
        :return: A list of dictionaries. Each dictionary contains the information for a single sample.
        """
        data = {
            "name": name,
            "sampleid": sampleid,
            "refmol": refmol,
            "refmolid": refmolid,
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
            "properties": properties,
        }
        fname = self.get_sample_fname(name)  # full path
        try:
            self.write_pickle(fname, data)
        except OSError:
            os.makedirs(os.path.dirname(fname), exist_ok=True)
            self.write_pickle(fname, data)

        self._samplefiles.add(fname)
        if algnid is None:
            self._samplefiles_to_profile.add(fname)
        return fname

    def iter_samples(self):
        for fname in self._samplefiles:
            yield self.read_pickle(fname)

    def cache_sequence(self, seqhash, sequence):
        fname = self.get_seq_fname(seqhash)
        if os.path.isfile(fname):
            if self.file_collision(fname, sequence):
                sys.exit(
                    "seqhash collision: sequences differ for seqhash " + seqhash + "."
                )
        else:
            try:
                with open(fname, "w") as handle:
                    handle.write(sequence)
            except OSError:
                os.makedirs(os.path.dirname(fname), exist_ok=True)
                with open(fname, "w") as handle:
                    handle.write(sequence)
        return fname

    def cache_reference(self, refid, sequence):
        fname = self.get_ref_fname(refid)
        if refid not in self._refs:
            with open(fname, "w") as handle:
                handle.write(sequence)
            self._refs.add(refid)
        return fname

    def cache_translation_table(self, translation_id, dbm):
        """
        If the translation table
        is not in the cache, it is retrieved from the database and written to a file

        :param translation_id: The id of the translation table
        :param dbm: the database manager
        :return: A file name.
        """
        fname = self.get_tt_fname(translation_id)  # write under /cache/ref/
        if translation_id not in self._tt:
            self.write_pickle(fname, dbm.get_translation_dict(translation_id))
            self._tt.add(translation_id)
        return fname

    def cache_lift(self, refid, refmol_acc, sequence):
        """
                        The function takes in a reference id, a reference molecule accession number,
                        and a reference sequence. It then checks to see if the reference molecule accession number is in the set of molecules that
                        have been cached. If it is not, it iterates through all of the coding sequences for that molecule and creates a
                        dataframe for each one.
        .
                        It then saves the dataframe to a pickle file and adds the reference molecule accession number to
                        the set of molecules that have been cached.
                        It then returns the name of the pickle file

        """
        fname = self.get_lift_fname(refid)
        rows = []
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
            for cds in self.iter_cds(refmol_acc):
                elemid = cds["id"]
                symbol = cds["symbol"]
                seq = cds["sequence"] + "*"
                codon = 0
                i = 0
                coords = []
                for rng in cds["ranges"]:
                    coords.extend(list(rng))
                while len(coords) % 3 != 0:
                    coords.append("")
                coords_len = len(coords)
                while len(seq) < coords_len / 3:
                    seq += "-"
                while len(sequence) < coords_len:
                    sequence += "-"
                for i, coord in enumerate(
                    [coords[x : x + 3] for x in range(0, len(coords), 3)]
                ):
                    codon = [sequence[coord[0]], sequence[coord[1]], sequence[coord[2]]]
                    rows.append(
                        [elemid] + coord + codon + codon + [symbol, i, seq[i].strip()]
                    )
                df = pd.DataFrame.from_records(rows, columns=cols, coerce_float=False)
                df = df.reindex(df.columns.tolist(), axis=1)
                df.to_pickle(fname)
                if self.debug:
                    df.to_csv(fname + ".csv")
            self._lifts.add(refmol_acc)
        return fname

    @staticmethod
    def open_file(fname, mode="r", compressed="auto", encoding=None):
        if not os.path.isfile(fname):
            sys.exit("input error: " + fname + " does not exist.")
        if compressed == "auto":
            compressed = os.path.splitext(fname)[1][1:]
        try:
            if compressed == "gz":
                return gzip.open(fname, mode + "t", encoding=encoding)
            if compressed == "xz":
                return lzma.open(fname, mode + "t", encoding=encoding)
            return open(fname, mode, encoding=encoding)
        except Exception:
            sys.exit("input error: " + fname + " cannot be opened.")

    def process_fasta_entry(self, header, seq):
        sample_id = header.replace("\t", " ").replace("|", " ").split(" ")[0]
        refmol = self.get_refmol(header)
        if not refmol:
            sys.exit(
                "input error: "
                + sample_id
                + " refers to an unknown reference molecule ("
                + refmol
                + ")."
            )
        seq = sonarBasics.harmonize(seq)
        seqhash = sonarBasics.hash(seq)
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

    def iter_fasta(self, *fnames):
        """
        This function iterates over the fasta files and returns a dictionary for each record
        """
        for fname in fnames:
            with self.open_file(fname) as handle, tqdm(
                desc="processing " + fname + "...",
                total=os.path.getsize(fname),
                unit="bytes",
                unit_scale=True,
                bar_format="{desc} {percentage:3.0f}% [{n_fmt}/{total_fmt}, {elapsed}<{remaining}, {rate_fmt}{postfix}]",
                disable=self.disable_progress,
            ) as pbar:
                seq = []
                for line in handle:
                    pbar.update(len(line))
                    line = line.strip()
                    if line.startswith(">"):
                        header = line[1:]
                        if seq:
                            yield self.process_fasta_entry(header, "".join(seq))
                            seq = []
                    else:
                        seq.append(line)
                if seq:
                    yield self.process_fasta_entry(header, "".join(seq))

    def get_refmol(self, fasta_header):
        mol = self._molregex.search(fasta_header)
        if not mol:
            try:
                return self.refmols[mol]["accession"]
            except Exception:
                None
        return self.default_refmol_acc

    def get_refseq(self, refmol_acc):
        try:
            return self.sources[refmol_acc]["sequence"]
        except Exception:
            return None

    def iter_cds(self, refmol_acc):
        cds = {}
        prev_elem = None
        with sonarDBManager(self.db, debug=self.debug) as dbm:
            for row in dbm.get_annotation(
                molecule_accession=refmol_acc, element_type="cds"
            ):
                if prev_elem is None:
                    prev_elem = row["element.id"]
                elif row["element.id"] != prev_elem:
                    yield cds
                    cds = {}
                    prev_elem = row["element.id"]
                if cds == {}:
                    cds = {
                        "id": row["element.id"],
                        "accession": row["element.accession"],
                        "symbol": row["element.symbol"],
                        "sequence": row["element.sequence"],
                        "ranges": [
                            range(
                                row["elempart.start"],
                                row["elempart.end"],
                                row["elempart.strand"],
                            )
                        ],
                    }
                else:
                    cds["ranges"].append(
                        range(
                            row["elempart.start"],
                            row["elempart.end"],
                            row["elempart.strand"],
                        )
                    )
        if cds:
            yield cds

    def get_refseq_id(self, refmol_acc):
        try:
            return self.sources[refmol_acc]["id"]
        except Exception:
            return None

    def get_refhash(self, refmol_acc):
        try:
            if "seqhash" not in self.sources[refmol_acc]:
                self.sources[refmol_acc]["seqhash"] = sonarBasics.hash(
                    self.sources[refmol_acc]["sequence"]
                )
            return self.sources[refmol_acc]["seqhash"]
        except Exception:
            return None

    def get_properties(self, fasta_header):
        return {x.group(1): x.group(2) for x in self._propregex.finditer(fasta_header)}

    def add_fasta(self, *fnames, propdict=defaultdict(dict)):  # noqa: C901
        default_properties = {
            x: self.properties[x]["standard"] for x in self.properties
        }
        with sonarDBManager(self.db, debug=self.debug) as dbm:
            for fname in fnames:
                for data in self.iter_fasta(fname):
                    # check sample
                    data["sampleid"], seqhash = dbm.get_sample_data(data["name"])
                    data["sourceid"] = dbm.get_source(data["refmolid"])["id"]

                    # check properties
                    if data["sampleid"] is None:
                        props = default_properties
                        for k, v in data["properties"].items():
                            props[k] = v
                        for k, v in propdict[data["sampleid"]].items():
                            props[k] = v
                        data["properties"] = props
                    elif not self.allow_updates:
                        continue
                    else:
                        for k, v in propdict[data["sampleid"]].items():
                            data["properties"][k] = v

                    # check reference
                    refseq_id = self.get_refseq_id(data["refmol"])
                    if not refseq_id:
                        if not self.ignore_errors:
                            self.log(
                                "fasta header refers to an unknown refrence ("
                                + data["header"]
                                + ")",
                                True,
                                "input error",
                            )
                        else:
                            self.log(
                                "skipping "
                                + data["name"]
                                + " referring to an unknown reference ("
                                + data["header"]
                                + ")"
                            )

                    # check alignment
                    data["algnid"] = dbm.get_alignment_id(data["seqhash"], refseq_id)
                    if data["algnid"] is None:
                        data["seqfile"] = self.cache_sequence(
                            data["seqhash"], data["sequence"]
                        )
                        data["reffile"] = self.cache_reference(
                            refseq_id, self.get_refseq(data["refmol"])
                        )
                        data["ttfile"] = self.cache_translation_table(
                            data["translation_id"], dbm
                        )
                        data["liftfile"] = self.cache_lift(
                            refseq_id, data["refmol"], self.get_refseq(data["refmol"])
                        )
                        data["algnfile"] = self.get_algn_fname(
                            data["seqhash"] + "@" + self.get_refhash(data["refmol"])
                        )
                        data["varfile"] = self.get_var_fname(
                            data["seqhash"] + "@" + self.get_refhash(data["refmol"])
                        )
                    elif data["seqhash"] != seqhash:
                        data["seqfile"] = self.cache_sequence(
                            data["seqhash"], data["sequence"]
                        )
                        data["reffile"] = None
                        data["ttfile"] = None
                        data["liftfile"] = None
                        data["algnfile"] = None
                        data["varfile"] = None
                    else:
                        data["seqhash"] = None
                        data["seqfile"] = None
                        data["reffile"] = None
                        data["ttfile"] = None
                        data["liftfile"] = None
                        data["algnfile"] = None
                        data["varfile"] = None
                    del data["sequence"]
                    self.cache_sample(**data)

    def import_cached_samples(self):  # noqa: C901
        refseqs = {}
        with sonarDBManager(self.db, readonly=False, debug=self.debug) as dbm:
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
                    if not sample_data["seqhash"] is None:
                        dbm.insert_sample(sample_data["name"], sample_data["seqhash"])
                        algnid = dbm.insert_alignment(
                            sample_data["seqhash"], sample_data["refmolid"]
                        )
                    if not sample_data["var_file"] is None:
                        with open(sample_data["var_file"], "r") as handle:
                            for line in handle:
                                if line == "//":
                                    break
                                vardat = line.strip("\r\n").split("\t")
                                dbm.insert_variant(
                                    algnid,
                                    vardat[4],
                                    vardat[0],
                                    vardat[3],
                                    vardat[1],
                                    vardat[2],
                                    vardat[5],
                                )
                            if line != "//":
                                sys.exit(
                                    "cache error: corrupted file ("
                                    + sample_data["var_file"]
                                    + ")"
                                )
                    # paranoia test
                    if not sample_data["seqhash"] is None:
                        try:
                            seq = list(refseqs[sample_data["refmolid"]])
                        except Exception:
                            refseqs[sample_data["refmolid"]] = list(
                                dbm.get_sequence(sample_data["refmolid"])
                            )
                            seq = list(refseqs[sample_data["refmolid"]])

                        prefix = ""
                        for vardata in dbm.iter_dna_variants(
                            sample_data["name"], sample_data["refmolid"]
                        ):
                            if vardata["variant.alt"] == " ":
                                for i in range(
                                    vardata["variant.start"], vardata["variant.end"]
                                ):
                                    seq[i] = ""
                            elif vardata["variant.start"] >= 0:
                                seq[vardata["variant.start"]] = vardata["variant.alt"]
                            else:
                                prefix = vardata["variant.alt"]
                        seq = prefix + "".join(seq)

                        with open(sample_data["seq_file"], "r") as handle:
                            orig_seq = handle.read()

                        if seq != orig_seq:
                            aligner = sonarAligner()
                            with TemporaryDirectory() as tempdir:
                                qryfile = os.path.join(tempdir, "qry")
                                reffile = os.path.join(tempdir, "ref")
                                with open(qryfile, "w") as handle:
                                    handle.write(">seq\n" + seq)
                                with open(reffile, "w") as handle:
                                    handle.write(">ref\n" + orig_seq)
                                qry, ref = aligner.align(qryfile, reffile)
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
                            sys.exit(
                                "import error: original sequence of sample "
                                + sample_data["name"]
                                + " cannot be restored from stored genomic profile for sample (see paranoid.alignment.fna)"
                            )

                except Exception as e:
                    print("\n------- Fatal Error ---------")
                    print("\nDebugging Information")
                    print(e)
                    pp.pprint(sample_data)
                    sys.exit("unknown import error")


if __name__ == "__main__":
    pass
