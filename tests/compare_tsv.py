#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

import argparse
import re
import sys

import pandas as pd


def check_format(file_path):
    with open(file_path, "r") as f:
        first_line = f.readline()
    if first_line.startswith("sample.name"):
        version = 1
        col_delim = first_line[11]
        val_delim = " "
    elif first_line.startswith("SAMPLE_NAME"):
        version = 2
        col_delim = first_line[11]
        val_delim = ","
    else:
        sys.exit(f"file format unknown for {file_path}")
    return version, col_delim, val_delim


def sort_delimited_values(s, delimiter):
    if pd.isna(s):
        return s
    items = s.split(delimiter)
    items.sort()
    return delimiter.join(items)


def replace_cds_symbols(value, mapping):
    mapping = {"ORF1a": "ORF1ab", "ORF1b": "ORF1ab"}
    for key, val in mapping.items():
        value = value.replace(key, val)
    return value


def main(args):  # noqa: C901
    err_code = 0
    # config cols
    profile_cols = ["GENOMIC_PROFILE", "PROTEOMIC_PROFILE", "FRAMESHIFT_MUTATIONS"]
    column_mapping = {
        "sample.name": "SAMPLE_NAME",
        "NUC_PROFILE": "GENOMIC_PROFILE",
        "AA_PROFILE": "PROTEOMIC_PROFILE",
        "FRAMESHIFTS": "FRAMESHIFT_MUTATIONS",
    }

    # checking format
    f1_version, f1_col_delim, f1_val_delim = check_format(args.file1)
    f2_version, f2_col_delim, f2_val_delim = check_format(args.file2)

    print("File 1:", args.file1, f"(cs {f1_version} format)")
    print("File 2:", args.file2, f"(cs {f2_version} format)")

    # Read the DataFrames with specified delimiters
    df1 = pd.read_csv(args.file1, delimiter=f1_col_delim)
    df2 = pd.read_csv(args.file2, delimiter=f2_col_delim)

    # Mapping columns v1 -> v2
    if f1_version == 1:
        df1.rename(columns=column_mapping, inplace=True)
    if f2_version == 1:
        df2.rename(columns=column_mapping, inplace=True)

    # Check column names
    missing_columns_df1 = [col for col in df2.columns if col not in df1.columns]
    missing_columns_df2 = [col for col in df1.columns if col not in df2.columns]

    if missing_columns_df1:
        print(f"Warning: Columns {missing_columns_df1} are missing in file 1.")
    if missing_columns_df2:
        print(f"Warning: Columns {missing_columns_df2} are missing in file 2.")
    if missing_columns_df1 or missing_columns_df2:
        print("Note: skipping missing columns when comparing.")

    common_cols = df1.columns.intersection(df2.columns)
    if len(common_cols) == 0:
        print("No common columns to compare.")
        return

    # Check data types
    mismatched_dtypes = []
    for col in common_cols:
        if df1[col].dtype != df2[col].dtype:
            mismatched_dtypes.append((col, df1[col].dtype, df2[col].dtype))

    if mismatched_dtypes:
        err_code = 1
        print("Data types do not match for the following columns:")
        for col, dtype1, dtype2 in mismatched_dtypes:
            print(f"  Column: {col}, File 1: {dtype1}, File 2: {dtype2}")
    else:
        print("Data types match for all common columns.")

    # Check for SAMPLE_NAME column
    if "SAMPLE_NAME" not in df1.columns or "SAMPLE_NAME" not in df2.columns:
        print(
            "SAMPLE_NAME column not found in one or both files. Skipping value checks."
        )
        sys.exit(1)

    # Check 'SAMPLE_NAME' column values and order
    if df1["SAMPLE_NAME"].equals(df2["SAMPLE_NAME"]):
        print("The SAMPLE_NAME columns contain the same values in the same order.")
    else:
        if set(df1["SAMPLE_NAME"]) == set(df2["SAMPLE_NAME"]):
            print(
                "The SAMPLE_NAME columns contain the same values but in a different order."
            )
        else:
            print("The SAMPLE_NAME columns contain different values.")

    # Perform an inner join on 'SAMPLE_NAME'
    merged_df = pd.merge(df1, df2, on="SAMPLE_NAME", suffixes=("_file1", "_file2"))

    if merged_df.empty:
        print("No common SAMPLE_NAME values found to compare.")
        sys.exit(1)

    # check properties
    discrepancies = []
    for _, row in merged_df.iterrows():
        sample_name = row["SAMPLE_NAME"]
        for col in merged_df.columns:
            orig_col = col.replace("_file1", "").replace("_file2", "")
            if col == "SAMPLE_NAME" or orig_col in profile_cols:
                continue
            if "_file1" in col or "_file2" in col:
                val1 = row[f"{orig_col}_file1"]
                val2 = row[f"{orig_col}_file2"]

                if pd.notna(val1) and pd.notna(val2) and val1 != val2:
                    discrepancies.append(
                        f"\n  SAMPLE_NAME: {sample_name}"
                        + f"\n  Column: {orig_col}"
                        + f"\n  Value_in_File1: {val1}"
                        + f"\n  Value_in_File2: {val2}"
                    )

    if not discrepancies:
        print("The property columns contain the same profile information.")
    else:
        err_code = 1
        print("The property columns contain different values:\n")
        print("\n".join(discrepancies) + "\n\n")

    # check profiles
    for field in profile_cols:
        discrepancies = []
        for _, row in merged_df.iterrows():
            sample_name = row["SAMPLE_NAME"]
            file1_values = (
                str(row[f"{field}_file1"]) if pd.notna(row[f"{field}_file1"]) else ""
            )
            file2_values = (
                str(row[f"{field}_file2"]) if pd.notna(row[f"{field}_file2"]) else ""
            )
            file1_set = set(re.split("[, ]", file1_values)) if file1_values else set()
            file2_set = set(re.split("[, ]", file2_values)) if file2_values else set()

            if field == "PROTEOMIC_PROFILE":
                if f1_version == "output version 1":
                    file1_set = set(replace_cds_symbols(x) for x in file1_set)
                if f2_version == "output version 1":
                    file2_set = set(replace_cds_symbols(x) for x in file2_set)

            missing_in_file1 = file2_set - file1_set
            missing_in_file2 = file1_set - file2_set

            if missing_in_file1 or missing_in_file2:
                discrepancies.append(
                    f"  SAMPLE_NAME {sample_name}\n  Missing_in_File1': "
                    + ", ".join(missing_in_file1)
                    + "\n  Missing_in_File2': "
                    + ", ".join(missing_in_file2)
                )
        if not discrepancies:
            print(f"The {field} columns contain the same profile information.")
        else:
            err_code = 1
            print(f"The {field} columns contain different values:\n")
            print("\n\n".join(discrepancies) + "\n\n")

        if err_code != 0:
            sys.exit(err_code)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare two CSV files.")
    parser.add_argument("--file1", required=True, help="Path to the first CSV file.")
    parser.add_argument("--file2", required=True, help="Path to the second CSV file.")

    args = parser.parse_args()
    main(args)
