#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

import os
import pickle
import re
import sys
from typing import Dict
from typing import Generator
from typing import List
from typing import Optional
from typing import Pattern
from typing import Tuple

import pandas as pd
import parasail

from covsonar.logging import LoggingConfigurator

# Initialize logger
LOGGER = LoggingConfigurator.get_logger()


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
        self._cigar_pattern: Pattern[str] = re.compile(r"(\d+)([M=XIDSHNP])")

    def read_seqcache(self, fname: str) -> str:
        """
        Reads a sequence cache file and returns the sequence stored in it.

        Args:
            fname (str): The filename of the sequence cache file to read.

        Returns:
            str: The sequence stored in the sequence cache file.

        Raises:
            FileNotFoundError: If the specified file does not exist.

        """
        with open(fname, "r") as handle:
            seq = handle.readline().strip()
        return seq

    def align(
        self, qryseq: str, refseq: str, gapopen: int = 16, gapextend: int = 4
    ) -> Tuple[str, str, str]:
        """
        Align two sequences using semi-global Smith-Waterman algorithm provided by PARASAIL and return the aligned sequences.

        Args:
            qryseq (str): Query sequence to align.
            refseq (str): Reference sequence to align against.
            gapopen (int, optional): Gap open penalty. Defaults to 16.
            gapextend (int, optional): Gap extension penalty. Defaults to 4.

        Returns:
            Tuple[str, str, str]: Tuple containing the aligned reference sequence,
                                  aligned query sequence, and CIGAR string.
        """
        result = parasail.sg_trace_striped_32(
            qryseq, refseq, gapopen, gapextend, parasail.blosum62
        )

        # Extract the aligned sequences and CIGAR string from the parasail result
        aligned_refseq = result.traceback.ref
        aligned_qryseq = result.traceback.query
        cigar_string = result.get_cigar().decode.decode()

        return aligned_refseq, aligned_qryseq, cigar_string

    def parse_cigar(  # noqa: C901
        self,
        elemid: int,
        reference_seq: str,
        query_seq: str,
        cigar: str,
        cds_df: Optional[pd.DataFrame] = None,
        del_alignment: Optional[str] = None,
    ) -> List[Tuple[int, str, str, int, int, int, str]]:
        """
        Extracts variant information from a reference and query sequence, given their aligned CIGAR string.

        Args:
            elemid (int): The element ID.
            reference_seq (str): The reference sequence.
            query_seq (str): The query sequence.
            cigar (str): The CIGAR string representing the alignment between the reference and query sequences.
            cds_set (set): set of coding positions (without end position), optional.
            del_alignment (str): how to align indels (left or right) or None for no alignment.

        Returns:
            list: A list of tuples containing the element ID, reference base, alternate base, start position, end position, frameshift indicator, and mutation label for each variant.

        Example:
            >>> sa = sonarAligner()
            >>> sa.parse_cigar(0, "AAGTCC", "AACTCG", "2=1X2=1X")
            [(0, 'G', 'C', 2, 3, 0, 'G3C'), (0, 'C', 'G', 5, 6, 0, 'C6G')]
        """

        def is_frameshift_insert(alt: str, ref_pos: int) -> int:
            """Checks if a snp or insert results in an framsehift."""
            nonlocal cds_set
            if ref_pos in cds_set and len(alt) % 3 != 0:
                return 1
            return 0

        def is_frameshift_del(ref_start: int, ref_end: int) -> int:
            """Checks if a snp or insert results in an framsehift."""
            nonlocal cds_df
            coords = set(range(ref_start, ref_end))
            cds_df["gap"] = cds_df["pos"].isin(coords).astype(int)
            groups = cds_df.groupby(
                [cds_df["elemid"], (cds_df["gap"].shift() != cds_df["gap"]).cumsum()]
            )["gap"].sum()
            return int(any(g % 3 != 0 for g in groups))

        def align_dels(
            direction: str,
            ref_pos: int,
            ref_end: int,
            ref_base: str,
            reference_seq: str,
        ) -> Tuple[int, int, str]:
            """Aligns deletions based on the specified mode."""
            if direction == "left":
                while (
                    ref_pos > 0
                    and reference_seq[ref_pos - 1] == reference_seq[ref_end - 1]
                ):
                    ref_pos -= 1
                    ref_end -= 1
                    ref_base = reference_seq[ref_pos:ref_end]
            elif direction == "right":
                while reference_seq[ref_pos + 1 :].startswith(ref_base):
                    ref_pos += 1
                    ref_end += 1
            elif direction is not None:
                LOGGER.error(
                    "'" + str(direction) + "' is an unknown indel alignment mode."
                )
                sys.exit(1)

            return ref_pos, ref_end, ref_base

        def handle_match(length: int):
            """Processes identical sites in query  and reference sequences."""
            nonlocal ref_pos, query_pos
            ref_pos += length
            query_pos += length

        def handle_snp(length: int):
            """Processes single nucleotide polymorphisms."""
            nonlocal ref_pos, query_pos, variants
            for _ in range(length):
                ref_base = reference_seq[ref_pos]
                alt_base = query_seq[query_pos]
                label = ref_base + str(ref_pos + 1) + alt_base
                variants.append(
                    (elemid, ref_base, alt_base, ref_pos, ref_pos + 1, 0, label)
                )
                ref_pos += 1
                query_pos += 1

        def handle_deletion(length: int):
            """Processes deletions in the query sequence."""
            nonlocal ref_pos, variants
            ref_end = ref_pos + length
            ref_base = reference_seq[ref_pos:ref_end]
            alt_base = " " if ref_pos != 0 else "."

            if del_alignment is not None:
                ref_pos, ref_end, ref_base = align_dels(
                    del_alignment, ref_pos, ref_end, ref_base, reference_seq
                )

            fs = is_frameshift_del(ref_pos, ref_end)

            if ref_end - ref_pos == 1:
                label = "del:" + str(ref_pos + 1)
            else:
                label = "del:" + str(ref_pos + 1) + "-" + str(ref_end)

            variants.append((elemid, ref_base, alt_base, ref_pos, ref_end, fs, label))
            ref_pos += length

        def handle_insertion(length: int):
            """Processes insertions in the query sequence."""
            nonlocal ref_pos, query_pos, variants
            ref_base = reference_seq[ref_pos - 1] if ref_pos > 0 else "."
            alt_base = (
                ref_base.replace(".", "") + query_seq[query_pos : query_pos + length]
            )

            fs = is_frameshift_insert(alt_base, ref_pos)

            label = ref_base + str(ref_pos) + alt_base
            variants.append(
                (elemid, ref_base, alt_base, ref_pos - 1, ref_pos, fs, label)
            )
            query_pos += length

        def handle_clipping(length: int):
            """Processes soft and hard clipping."""
            nonlocal query_pos
            query_pos += length

        def handle_skipped_bases(length: int):
            """Processes skipped bases in the reference sequence."""
            nonlocal ref_pos
            ref_pos += length

        operations = {
            "M": handle_match,
            "=": handle_match,
            "X": handle_snp,
            "D": handle_deletion,
            "I": handle_insertion,
            "S": handle_clipping,
            "H": handle_clipping,
            "N": handle_skipped_bases,
            "P": lambda length: (
                ref_pos,
                query_pos,
            ),  # Padding does not consume bases in either sequence, so do nothing
        }

        ref_pos = 0
        query_pos = 0
        variants = []
        cds_set = set() if cds_df is None else set(cds_df[cds_df["end"] == 0])

        # Parse the CIGAR string and process each operation
        for length, operation in self._cigar_pattern.findall(cigar):
            length = int(length)
            if operation in operations:
                operations[operation](length)
            else:
                LOGGER.critical("Unknown operator '" + operation + "' in cigar string")
                sys.exit(1)

        if cigar.endswith("D"):
            v = list(variants[-1])
            v[2] = "."
            variants[-1] = v

        return variants

    def process_cached_sample(self, fname: str) -> bool:
        """
        Processes a cached sample data file, aligns sequences, extracts and lifts variants and writes them to a file.

        Args:
            fname (str): the path to the file containing the sample data to process.

        Returns:
            bool: True if the operation is successful, else False.
        """

        # Load data from sample file
        with open(fname, "rb") as handle:
            data = pickle.load(handle, encoding="bytes")

        # If "var_file" is None, return True
        if data.get("var_file") is None:
            return True

        # Verify if content of existing var file is complete, if so return True
        if os.path.isfile(data["var_file"]):
            with open(data["var_file"], "r") as handle:
                for line in handle:
                    pass
                if line == "//":
                    return True

        # Align sequence data
        elemid = str(data["sourceid"])
        qryseq = self.read_seqcache(data["seq_file"])
        refseq = self.read_seqcache(data["ref_file"])
        _, _, cigar = self.align(qryseq, refseq)

        # get annotation for frameshift detection
        cds_df = pd.read_pickle(data["cds_file"])

        # Extract nucleotide vars
        nuc_vars = self.parse_cigar(elemid, refseq, qryseq, cigar, cds_df)

        if not nuc_vars:
            vars = "//"
        else:
            # Create AA mutation and concatenate to the same file of NT variants
            aa_vars = [
                x for x in self.lift_vars(nuc_vars, data["lift_file"], data["tt_file"])
            ]
            vars = (
                "\n".join(
                    "\t".join([str(y) for y in x]) for x in nuc_vars + sorted(aa_vars)
                )
                + "\n//"
            )

        # Write variants to the "var_file"
        os.makedirs(os.path.dirname(data["var_file"]), exist_ok=True)
        with open(data["var_file"], "w") as handle:
            handle.write(vars)

        return True

    def translate(self, seq: str, tt: int):
        """
        Translates a given nucleotide sequence into amino acid code using a given translation table

        Parameters
        ----------
        seq : str
            The nucleotide sequence.
        tt : int
            The ID of the tranlation table to use.

        Returns
        -------
        str
            The translated amino acid sequence
        """
        aa = []
        for codon in [seq[i : i + 3] for i in range(0, (len(seq) // 3) * 3, 3)]:
            aa.append(tt[codon])
            if tt[codon] == "*":
                break
        return "".join(aa)

    def update_nuc_positions(
        self,
        df: pd.DataFrame,
        nuc_vars: List[Tuple[int, str, str, int, int, int, str]],
    ) -> pd.DataFrame:
        """
        Update nucleotide positions with alternate nucleotides

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame with genomic information.
        nuc_vars : List[Tuple[int, str, str, int, int, int, str]]
            A list of tuples, where each tuple contains the element ID, reference allele, alterante allele, start position, end position, frameshift indicator, and mutation label.

        Returns
        -------
        pd.DataFrame
            Updated DataFrame with alternate nucleotides.
        """
        # Update nucleotide positions
        rng = ["1", "2", "3"]
        for nuc_var in nuc_vars:
            alt = "-" if nuc_var[2] == " " else nuc_var[2]
            if alt == ".":
                continue  # ignore uncovered terminal regions
            x, y = int(nuc_var[3]), int(nuc_var[4]) - 1
            for i in rng:
                df.loc[df["nucPos" + i].between(x, y), "alt" + i] = alt
        return df

    def filter_nuc_positions(
        self, df: pd.DataFrame, tt: Dict[str, str]
    ) -> pd.DataFrame:
        """
        Filter the dataframe rows based on the alternate and reference nucleotides and translate the alternate nucleotides.

        :param df: The dataframe to be filtered.
        :param tt: Codon translation table.
        :return: Filtered dataframe.
        """
        # Filter rows where alternate and reference nucleotides are different
        df = df.copy().loc[
            (df["ref1"] != df["alt1"])
            | (df["ref2"] != df["alt2"])
            | (df["ref3"] != df["alt3"])
        ]

        # Translate alternate nucleotides and create new column
        if not df.empty:
            df["altAa"] = df.apply(
                lambda x: self.translate(x["alt1"] + x["alt2"] + x["alt3"], tt), axis=1
            )

            # Further filter dataframe to only include rows where aa and altAa are different
            df = df.loc[df["aa"] != df["altAa"]]
        else:
            df["altAa"] = pd.Series(dtype="object")
        return df

    def lift_vars(
        self,
        nuc_vars: List[Tuple[int, str, str, int, int, int, str]],
        lift_file: str,
        tt_file: str,
    ) -> Generator[Tuple[int, str, str, int, int, int, str], None, None]:
        """
        Lift over nucleotide variants to protein (amino acid) variants.

        Parameters
        ----------
        nuc_vars : List[Tuple[int, str, str, int, int, int]]
            A list of tuples, where each tuple contains the element ID, reference base, alternate base, start position, end position, frameshift indicator, and mutation label for each variant.
        lift_file : str
            File path to a pickled DataFrame with genomic information.
        tt_file : str
            File path to a pickled translation table.

        Yields
        -------
        Tuple[int, str, str, int, int, int, str]]
            A tuple containing the element ID, reference base, alternate base, start position, end position, frameshift indicator, and mutation label for each variant.
        """

        def get_aa_sub_label(df_tuple):
            """Returns the label for amino acid substituions"""
            return df_tuple.aa + str(row.aaPos + 1) + df_tuple.altAa

        def get_aa_del_label(start, end):
            """Returns the label for amino acid deletions"""
            if end - start == 1:
                return "del:" + str(start + 1)
            return "del:" + str(start + 1) + "-" + str(end)

        df = pd.read_pickle(lift_file)
        with open(tt_file, "rb") as handle:
            tt = pickle.load(handle, encoding="bytes")

        # Update and filter nucleotide positions with alternate nucleotides
        pd.set_option("display.max_columns", None)
        pd.set_option("display.max_rows", None)
        pd.set_option("display.expand_frame_repr", False)

        df = self.update_nuc_positions(df, nuc_vars)
        df = self.filter_nuc_positions(df, tt)

        # handle snps or inserts
        filtered_df = df[
            (df["altAa"] != "-") & (df["altAa"] != " ") & (df["altAa"] != "")
        ]
        for row in filtered_df.itertuples():
            pos = row.aaPos
            label = get_aa_sub_label(row)
            yield row.elemid, row.aa, row.altAa, pos, pos + 1, 0, label

        # handle deletions
        filtered_df = df[
            (df["altAa"] == "-") | (df["altAa"] == " ") | (df["altAa"] == "")
        ].sort_values(["elemid", "aaPos"])
        prev_row = None
        for row in filtered_df.itertuples():
            if prev_row is None:
                prev_row = row
                prev_aa = [row.aa]
                prev_pos = row.aaPos
            elif prev_row.elemid == row.elemid and row.aaPos - prev_pos == 1:
                prev_aa.append(row.aa)
                prev_pos = row.aaPos
            else:
                start = prev_row.aaPos
                end = start + len(prev_aa)
                label = get_aa_del_label(start, end)
                yield prev_row.elemid, "".join(prev_aa), " ", start, end, 0, label
                prev_row = row
                prev_aa = [row.aa]
                prev_pos = row.aaPos

        if prev_row is not None:
            start = prev_row.aaPos
            end = start + len(prev_aa)
            label = get_aa_del_label(start, end)
            yield prev_row.elemid, "".join(prev_aa), " ", start, end, 0, label
