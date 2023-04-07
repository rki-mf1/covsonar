#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

import os
import pickle
import re
import sys

from Bio.Emboss.Applications import StretcherCommandline
import pandas as pd
import parasail


class sonarAligner(object):
    """
    this object performs a pairwise sequence alignment and provides/stores selected
    alignment functionalities/statistics.
    """

    def __init__(self):
        self.nuc_profile = []
        self.nuc_n_profile = []
        self.aa_profile = []
        self.aa_n_profile = []
        self.cigar_pattern = re.compile(r"(\d+)(\D)")

    def read_seqcache(self, fname):
        with open(fname, "r") as handle:
            seq = handle.readline().strip()
        return seq

    def align_global(self, qry, ref, gapopen=16, gapextend=4):
        """ """
        cline = StretcherCommandline(
            asequence=qry,
            bsequence=ref,
            gapopen=gapopen,
            gapextend=gapextend,
            outfile="stdout",
            aformat="fasta",
        )
        stdout, stderr = cline()
        s1 = stdout.find("\n") + 1
        e = stdout[1:].find(">") + 1
        s2 = stdout[e:].find("\n") + e
        qry = stdout[s1:e].replace("\n", "")
        ref = stdout[s2:].replace("\n", "")
        return qry, ref

    def align(self, qryseq, refseq, gapopen=16, gapextend=4):
        """ """
        result = parasail.sg_trace(
            qryseq, refseq, gapopen, gapextend, parasail.blosum62
        )
        return (
            result.traceback.ref,
            result.traceback.query,
            result.get_cigar().decode.decode(),
        )

    def extract_vars_from_cigar(  # noqa: C901
        self, qryseq, refseq, cigar, elemid, cds_file
    ):

        # get annotation for frameshift detection
        cds_df = pd.read_pickle(cds_file)
        cds_set = set(cds_df[cds_df["end"] == 0])

        # extract
        refpos = 0
        qrypos = 0
        prefix = False
        qrylen = len(qryseq)
        vars = []
        for match in self.cigar_pattern.finditer(cigar):
            vartype = match.group(2)
            varlen = int(match.group(1))
            # identical sites
            if vartype == "=":
                refpos += varlen
                qrypos += varlen
            # snp handling
            elif vartype == "X":
                for x in range(varlen):
                    ref = refseq[refpos]
                    alt = qryseq[qrypos]
                    vars.append(
                        (
                            ref,
                            str(refpos),
                            str(refpos + 1),
                            alt,
                            elemid,
                            ref + str(refpos + 1) + alt,
                            "0",
                        )
                    )
                    refpos += 1
                    qrypos += 1
            # deletion handling
            elif vartype == "D":
                if (
                    refpos == 0 and prefix is False
                ) or qrypos == qrylen:  # deletion at sequence terminus
                    vars.append(
                        (
                            refseq[refpos : refpos + varlen],
                            str(refpos),
                            str(refpos + varlen),
                            ".",
                            elemid,
                            "del:" + str(refpos + 1) + "-" + str(refpos + varlen),
                            "0",
                        )
                    )
                    refpos += varlen
                elif varlen == 1:  # 1bp deletion
                    vars.append(
                        (
                            refseq[refpos],
                            str(refpos),
                            str(refpos + 1),
                            " ",
                            elemid,
                            "del:" + str(refpos + 1),
                            self.detect_frameshifts(
                                refpos, refpos + 1, " ", cds_df, cds_set
                            ),
                        )
                    )
                    refpos += 1
                else:  # multi-bp inner deletion
                    vars.append(
                        (
                            refseq[refpos : refpos + varlen],
                            str(refpos),
                            str(refpos + varlen),
                            " ",
                            elemid,
                            "del:" + str(refpos + 1) + "-" + str(refpos + varlen),
                            self.detect_frameshifts(
                                refpos, refpos + varlen, " ", cds_df, cds_set
                            ),
                        )
                    )
                    refpos += varlen
            # insertion handling
            elif vartype == "I":
                if refpos == 0:  # to consider insertion before sequence start
                    ref = "."
                    alt = qryseq[:varlen]
                    prefix = True
                    fs = "0"
                else:
                    ref = refseq[refpos - 1]
                    alt = qryseq[qrypos - 1 : qrypos + varlen]
                    fs = self.detect_frameshifts(
                        refpos - 1, refpos, alt, cds_df, cds_set
                    )
                vars.append(
                    (
                        ref,
                        str(refpos - 1),
                        str(refpos),
                        alt,
                        elemid,
                        ref + str(refpos) + alt,
                        fs,
                    )
                )
                qrypos += varlen
            # unknown
            else:
                sys.exit(
                    "error: covSonar cannot interpret '"
                    + vartype
                    + "' in cigar string."
                )
        return vars

    def process_cached_sample(self, fname):
        """
        This function takes a sample file and processes it.
        create var file with NT and AA mutations
        """

        with open(fname, "rb") as handle:
            data = pickle.load(handle, encoding="bytes")

        if data["var_file"] is None:
            return True
        elif os.path.isfile(data["var_file"]):
            with open(data["var_file"], "r") as handle:
                for line in handle:
                    pass
            if line == "//":
                return True

        elemid = str(data["sourceid"])
        qryseq = self.read_seqcache(data["seq_file"])
        refseq = self.read_seqcache(data["ref_file"])
        _, __, cigar = self.align(qryseq, refseq)

        nuc_vars = [
            x
            for x in self.extract_vars_from_cigar(
                qryseq, refseq, cigar, elemid, data["cds_file"]
            )
        ]
        vars = "\n".join(["\t".join(x) for x in nuc_vars])
        if nuc_vars:
            # create AA mutation
            aa_vars = "\n".join(
                [
                    "\t".join(x)
                    for x in self.lift_vars(
                        nuc_vars, data["lift_file"], data["tt_file"]
                    )
                ]
            )
            if aa_vars:
                # concatenate to the same file of NT variants
                vars += "\n" + aa_vars
            vars += "\n"
        try:
            with open(data["var_file"], "w") as handle:
                handle.write(vars + "//")
        except OSError:
            os.makedirs(os.path.dirname(data["var_file"]), exist_ok=True)
            with open(data["var_file"], "w") as handle:
                handle.write(vars + "//")
        return True

    def translate(self, seq, tt):
        aa = []
        while len(seq) % 3 != 0:
            seq = seq[: len(seq) - 1]
        for codon in [seq[i : i + 3] for i in range(0, len(seq), 3)]:
            aa.append(tt[codon])
        return "".join(aa)

    def detect_frameshifts(self, start, end, alt, cds_df, coding_sites_set):
        # handling deletions
        if alt == " " or alt == "":
            coords = set(range(start, end))
            cds_df["gap"] = cds_df.apply(lambda x: 1 if x.pos in coords else 0, axis=1)
            groups = cds_df.groupby(
                [cds_df["elemid"], (cds_df["gap"].shift() != cds_df["gap"]).cumsum()]
            ).agg({"gap": sum})
            for _, row in groups.iterrows():
                if row["gap"] % 3 != 0:
                    return "1"
        # handling insertions
        elif (len(alt) - 1) % 3 > 0 and start in coding_sites_set:
            return "1"
        return "0"

    def lift_vars(self, nuc_vars, lift_file, tt_file):  # noqa: C901
        df = pd.read_pickle(lift_file)
        with open(tt_file, "rb") as handle:
            tt = pickle.load(handle, encoding="bytes")
        for nuc_var in nuc_vars:
            if nuc_var[3] == ".":
                continue  # ignore uncovered terminal regions
            for i in range(int(nuc_var[1]), int(nuc_var[2])):
                alt = "-" if nuc_var[3] == " " else nuc_var[3]
                df.loc[df["nucPos1"] == i, "alt1"] = alt
                df.loc[df["nucPos2"] == i, "alt2"] = alt
                df.loc[df["nucPos3"] == i, "alt3"] = alt

        df = df.loc[
            (df["ref1"] != df["alt1"])
            | (df["ref2"] != df["alt2"])
            | (df["ref3"] != df["alt3"])
        ]
        df["altAa"] = df.apply(
            lambda x: self.translate(x["alt1"] + x["alt2"] + x["alt3"], tt), axis=1
        )
        df = df.loc[df["aa"] != df["altAa"]]

        # snps or inserts
        for index, row in df.loc[(df["altAa"] != "-") & (df["altAa"] != "")].iterrows():
            pos = row["aaPos"] + 1
            label = row["aa"] + str(pos) + row["altAa"]
            yield row["aa"], str(pos - 1), str(pos), row["altAa"], str(
                row["elemid"]
            ), label, "0"

        # deletions
        prev_row = None
        for index, row in (
            df.loc[(df["altAa"] == "-")].sort_values(["elemid", "aaPos"]).iterrows()
        ):
            if prev_row is None:
                prev_row = row
            elif (
                prev_row["elemid"] == row["elemid"]
                and abs(prev_row["aaPos"] - row["aaPos"]) == 1
            ):
                prev_row["aa"] += row["aa"]
            else:
                start = prev_row["aaPos"]
                end = prev_row["aaPos"] + len(prev_row["aa"])
                if end - start == 1:
                    label = "del:" + str(start + 1)
                else:
                    label = "del:" + str(start + 1) + "-" + str(end)
                yield prev_row["aa"], str(start), str(end), " ", str(
                    prev_row["elemid"]
                ), label, "0"
                prev_row = row
        if prev_row is not None:
            start = prev_row["aaPos"]
            end = prev_row["aaPos"] + len(prev_row["aa"])
            if end - start == 1:
                label = "del:" + str(start + 1)
            else:
                label = "del:" + str(start + 1) + "-" + str(end)
            yield prev_row["aa"], str(start), str(end), " ", str(
                prev_row["elemid"]
            ), label, "0"
