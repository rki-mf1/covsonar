#!/usr/bin/python
# -*- coding: utf-8 -*-
from contextlib import ExitStack
import gzip
import math
from multiprocessing import Pool
import os
from os import getpid
import shutil
import subprocess
from tempfile import mkdtemp
import traceback
import warnings

from lib.sonardb import sonarDBManager
import numpy as np
import pandas as pd
from tqdm import tqdm

warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)


def create_fix_vcf_header(ref):
    header = "##fileformat=VCFv4.2\n##CreatedBy=covSonarV1.1.4\n##reference=" + ref
    format = '\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
    info = '\n##INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">'
    info = (
        info
        + '\n##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">'
    )
    info = (
        info
        + '\n##INFO=<ID=TYPE,Number=.,Type=String,Description="Mutation type e.g., SNP,INS and DEL">\n'
    )
    note = "##Note_1='Currently we ignore DEL type'\n"
    note = (
        note
        + "##Note_2='This VCF file is genereted by using var2vcf with betaV2, if you find any bugs, then please write a bug report to us'\n"
    )
    # column = "\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+sample_id+"\n"
    return header + format + info + note


def bgzip(filename):
    """Call bgzip to compress a file."""
    cmd = ["bgzip", "-f", filename]
    with subprocess.Popen(
        cmd, encoding="utf8", stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    ) as process:
        try:
            stdout, stderr = process.communicate(cmd)

        except subprocess.TimeoutExpired:
            process.kill()
            stdout, stderr = process.communicate()
            raise subprocess.TimeoutExpired(output=stdout, stderr=stderr)
        except Exception:
            process.kill()
            raise
    return


def tabix_index(filename):
    """Call tabix to create an index for a bgzip-compressed file."""
    cmd = ["tabix", "-p", "vcf", filename + ".gz"]
    with subprocess.Popen(
        cmd, encoding="utf8", stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    ) as process:
        try:
            stdout, stderr = process.communicate(cmd)
        except subprocess.TimeoutExpired:
            process.kill()
            stdout, stderr = process.communicate()
            raise subprocess.TimeoutExpired(output=stdout, stderr=stderr)
        except Exception:
            process.kill()
            raise
    return


def bcftool_index(filename):
    """Call tabix to create an index for a bgzip-compressed file."""
    cmd = ["bcftools", "index", filename]
    with subprocess.Popen(
        cmd, encoding="utf8", stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    ) as process:
        try:
            stdout, stderr = process.communicate(cmd)
        except subprocess.TimeoutExpired:
            process.kill()
            stdout, stderr = process.communicate()
            raise subprocess.TimeoutExpired(output=stdout, stderr=stderr)
        except Exception:
            process.kill()
            raise
    return


def calculate_AC_AN(final_df):
    # calculate AC AN
    # order-preserving index between POS and INFO AC
    # e.g. a.POS;b.POS a.AC,b.AC
    for row in final_df.itertuples():
        unique, counts = np.unique(
            np.asarray(row[11:]), return_counts=True
        )  # row[11:] means we start from sample ID column
        # for unique, counts in zip(unique, counts):
        AN = 0
        AC = ""
        for idx, val in enumerate(unique):
            if val == ".":  # ignore it
                continue
            else:
                _AC = counts[idx]
                AN = AN + _AC
                AC = str(_AC) if not AC else AC + "," + str(_AC)
        final_df.at[row.Index, "INFO"] = "AN=" + str(AN) + ";AC=" + AC
    return final_df


def _check_variant_type(ref, alt):
    if len(ref) == len(alt) and len(alt) == 1:  # SNP
        return "SNP"
    elif len(ref) < len(alt) and len(alt) > 0:  # INS
        if ref == alt[0]:
            return "INS"
        else:
            return "INDEL"
    elif len(ref) > len(alt) and len(ref) > 0:  # DEL
        return "DEL"
    else:
        print("Unknown:", ref, alt)
        return "Unknown"


def create_vcf(rows_grouped, tmp_dirname, refdescr):
    refdescr = refdescr.split()[0].replace(">", "")
    process_id = str(getpid())
    # print(process_id+" Start")
    # iterate over each group
    final_snp_df = pd.DataFrame(
        {
            "#CHROM": pd.Series(dtype="str"),
            "POS": pd.Series(dtype="int"),
            "ID": pd.Series(dtype="str"),
            "REF": pd.Series(dtype="str"),
            "ALT": pd.Series(dtype="str"),
            "FILTER": pd.Series(dtype="str"),
            "QUAL": pd.Series(dtype="str"),
            "INFO": pd.Series(dtype="str"),
            "FORMAT": pd.Series(dtype="str"),
            "TYPE": pd.Series(dtype="str"),
        },
        index=np.arange(1, 29904),
    )
    final_snp_df["POS"] = np.arange(1, 29904)
    final_snp_df[
        "#CHROM"
    ] = refdescr  # duplicated line here did a trick for assigning the value
    final_snp_df["ID"] = "."
    final_snp_df["REF"] = "."
    final_snp_df["ALT"] = "."
    final_snp_df["FILTER"] = "."
    final_snp_df["QUAL"] = "."
    final_snp_df["INFO"] = "."
    final_snp_df["FORMAT"] = "GT"
    final_snp_df["TYPE"] = "."

    final_indel_df = pd.DataFrame(
        {
            "#CHROM": pd.Series(dtype="str"),
            "POS": pd.Series(dtype="int"),
            "ID": pd.Series(dtype="str"),
            "REF": pd.Series(dtype="str"),
            "ALT": pd.Series(dtype="str"),
            "FILTER": pd.Series(dtype="str"),
            "QUAL": pd.Series(dtype="str"),
            "INFO": pd.Series(dtype="str"),
            "FORMAT": pd.Series(dtype="str"),
            "TYPE": pd.Series(dtype="str"),
        }
    )
    vcf_filename = process_id + ".vcf"
    full_path = os.path.join(tmp_dirname, vcf_filename)
    with open(full_path, "w") as f:
        # print("Create VCF file:",full_path)
        f.write(create_fix_vcf_header(refdescr))
        for group_name, df_ in tqdm(rows_grouped, mininterval=0.5):
            try:
                # final_snp_df[group_name] = "." # init value
                # final_indel_df[group_name] = "." # init value
                for row in df_.itertuples():
                    _id = (
                        getattr(row, "ref")
                        + str(getattr(row, "start"))
                        + getattr(row, "alt")
                    )
                    _type = _check_variant_type(
                        getattr(row, "ref"), getattr(row, "alt")
                    )

                    if _type == "SNP":
                        # find  ID
                        index_start_postion = getattr(row, "start")
                        selected_row = final_snp_df.loc[
                            index_start_postion
                        ]  # return scalar instead DF
                        if selected_row.ID == ".":  # A
                            final_snp_df.at[index_start_postion, "ID"] = _id
                            final_snp_df.at[index_start_postion, group_name] = "1"
                            final_snp_df.at[index_start_postion, "REF"] = getattr(
                                row, "ref"
                            )
                            final_snp_df.at[index_start_postion, "ALT"] = getattr(
                                row, "alt"
                            )
                            final_snp_df.at[index_start_postion, "TYPE"] = "SNP"

                        elif selected_row.ID != ".":  # B
                            if (
                                selected_row.ID == _id
                            ):  # only one ID exists and just update GT of a sample
                                final_snp_df.at[index_start_postion, group_name] = "1"
                            else:
                                splited_final_id_list = selected_row.ID.split(";")
                                totel_len = len(splited_final_id_list)
                                for new_GT, splited_final_id in enumerate(
                                    splited_final_id_list, start=1
                                ):
                                    if splited_final_id == _id:  # Found the exist one
                                        final_snp_df.at[
                                            index_start_postion, group_name
                                        ] = str(new_GT)
                                        break
                                    elif (
                                        totel_len == new_GT
                                    ):  # cannot find the same ID ,so we append the new one to the string
                                        final_snp_df.at[index_start_postion, "ID"] = (
                                            final_snp_df.at[index_start_postion, "ID"]
                                            + ";"
                                            + _id
                                        )
                                        # with new GT number
                                        final_snp_df.at[
                                            index_start_postion, group_name
                                        ] = str(new_GT + 1)
                                        # appends new alt
                                        final_snp_df.at[index_start_postion, "ALT"] = (
                                            final_snp_df.at[index_start_postion, "ALT"]
                                            + ","
                                            + getattr(row, "alt")
                                        )
                        else:
                            # get index
                            print(selected_row)
                            raise ValueError("Something went wrong")

                    # elif(_type == 'DEL'):
                    #    continue
                    elif (
                        _type == "INS" or _type == "INDEL" or _type == "DEL"
                    ):  # C      or _type == 'DEL'
                        selected_row = final_indel_df[
                            (final_indel_df["POS"] == getattr(row, "start"))
                            & (final_indel_df["TYPE"] == _type)
                        ]
                        if len(selected_row) == 0:  # D, always insert new records
                            new_row = {
                                "#CHROM": refdescr,
                                "ID": _id,
                                "POS": int(getattr(row, "start")),
                                "REF": getattr(row, "ref"),
                                "ALT": getattr(row, "alt"),
                                "TYPE": _type,
                                "FILTER": ".",
                                "QUAL": ".",
                                "INFO": ".",
                                "FORMAT": "GT",
                                group_name: "1",
                            }
                            final_indel_df = final_indel_df.append(
                                new_row, ignore_index=True
                            )
                        elif (
                            len(selected_row) > 0
                        ):  # E , found more than 0 which means we can update or insert it
                            for _E_rows in selected_row.itertuples():
                                index_postion = _E_rows.Index
                                splited_IDs_list = _E_rows.ID.split(";")
                                totel_len = len(splited_IDs_list)
                                for _E_new_GT, splited_ID in enumerate(
                                    splited_IDs_list, start=1
                                ):
                                    # print('index',index_postion,'test', splited_ID, _E_new_GT)
                                    if splited_ID == _id:  # Found the exist one
                                        # print('Found the exist one', splited_final_id)
                                        final_indel_df.at[
                                            index_postion, group_name
                                        ] = str(_E_new_GT)
                                        break
                                    elif (
                                        totel_len == _E_new_GT
                                    ):  # cannot find the same ID ,so we append the new one to the string
                                        final_indel_df.at[index_postion, "ID"] = (
                                            final_indel_df.at[index_postion, "ID"]
                                            + ";"
                                            + _id
                                        )
                                        # with new GT number
                                        final_indel_df.at[
                                            index_postion, group_name
                                        ] = str(_E_new_GT + 1)
                                        # appends new alt
                                        final_indel_df.at[index_postion, "ALT"] = (
                                            final_indel_df.at[index_postion, "ALT"]
                                            + ","
                                            + getattr(row, "alt")
                                        )

                    else:
                        print("Skip this ", _id)  # right now skip for Unknown
                        continue
            except ValueError as err:
                print(err.args)
                raise
            except Exception as e:
                print("An exception occurred at...")
                print(group_name)
                print(row)
                print(selected_row)
                print(traceback.format_exc())
                raise e

        final_df = pd.concat([final_snp_df, final_indel_df], axis=0, ignore_index=True)
        final_df = final_df.drop(final_df[final_df.ID == "."].index)
        final_df.replace(np.nan, ".", inplace=True)
        final_df["ALT"].replace("", ".", inplace=True)  # for deletion
        final_df = final_df.sort_values(["POS"], ascending=True)
        final_df = calculate_AC_AN(final_df)
        final_df["INFO"] = final_df["INFO"] + ";TYPE=" + final_df["TYPE"]
        final_df = final_df.drop(columns=["TYPE"])
        final_df.to_csv(f, sep="\t", encoding="utf-8", index=False)

    bgzip(full_path)
    tabix_index(full_path)
    return full_path + ".gz"
    # print(process_id+" Finish")


def parallelize_dataframe(df, tmp_dirname, num_cores, refdescr, func):
    _tmp_lis = np.array_split(df, num_cores)
    counter = 0
    zip_items = [
        (_tmp_lis[i], tmp_dirname, refdescr) for i in range(len(_tmp_lis))
    ]  # same order
    # with Pool(processes=num_cores) as pool:
    #    res = pool.starmap(func, zip_items)
    pool = Pool(num_cores)
    full_paht_list = pool.starmap(func, zip_items)

    # finish all tasks
    pool.close()
    pool.join()
    # print("Tmp result: ", full_paht_list)
    return full_paht_list


def export2VCF(
    db_path,
    include_acc,
    include_dates,
    output,
    num_cores,
    refdescr="No ref is provided",
):
    print("----- You are using sonartoVCF_V2 --------")
    print("WARNING: the function is still experimental/not fully implemented.")
    print("Prepare export2VCF workspace for", num_cores, "cpu")
    with ExitStack() as stack:
        dbm = stack.enter_context(sonarDBManager(db_path))

        where_clause = []
        where_vals = []

        if include_acc:
            where_clause.append(
                dbm.get_metadata_in_condition("accession", *include_acc)
            )
            where_vals.extend(include_acc)
        if include_dates:
            where_clause.append(dbm.get_metadata_date_condition("date", *include_dates))

        # print(where_clause)
        fields = "accession, start, end, alt, ref "
        if where_clause:
            sql = (
                "SELECT "
                + fields
                + " FROM dna_view WHERE "
                + " AND ".join(where_clause)
                + ";"
            )
        else:
            sql = "SELECT " + fields + " FROM dna_view;"

        # print("query: " + sql)
        # print("vals: ", where_vals)
        ##############################
        print("Start Bigquery...")
        rows = pd.read_sql(sql, dbm.connection, params=where_vals)
        print("Return:", len(rows), " records")
        track_vcf = []
        count = 0
        if not rows.empty:
            # rows.to_pickle("dummy.pkl")
            tmp_dirname = mkdtemp(prefix=".sonarCache_")
            # vcf_path=os.path.join(tmp_dirname,)

            # create fasta_id

            # rows['CHROM'] = chrom_id
            # rows['QUAL'] = '.'
            # rows['FILTER'] = '.'
            # rows['INFO'] = '.'
            # rows['FORMAT'] = 'GT'
            # POS or start position: The reference position, with the 1st base is position 1 not 0 , but in covsonar use 0 as the 1st position
            # so we should + 1
            # http://samtools.github.io/hts-specs/VCFv4.2.pdf
            rows["start"] = rows["start"] + 1
            # rows['end'] = rows['end']+1

            rows = rows.loc[
                (1 <= rows["start"]) & (rows["start"] <= 29903)
            ]  # filter out

            rows["alt"] = rows["alt"].replace("", np.nan)  # remove Deletion
            # rows['start'] = rows['start'].replace('', np.nan) # remove Insertion
            rows = rows.dropna(axis=0, subset=["alt"])
            rows_grouped = rows.groupby("accession")
            print("With :", len(rows_grouped), " accessions")
            # split data and write each ACC into individual VCF file.
            print("Start Divide and Conquer ...")
            track_vcf = parallelize_dataframe(
                rows_grouped, tmp_dirname, num_cores, refdescr, create_vcf
            )

            # bundle all vcf together
            print("Integrate all VCFs ...")
            divide_merge_vcf(track_vcf, output, num_cores)

            if os.path.isdir(tmp_dirname):
                shutil.rmtree(tmp_dirname)

            print("Finish! compress final result (gz):")


def divide_merge_vcf(list_track_vcf, global_output, num_cores):
    chunk = 500
    list_length = math.ceil(len(list_track_vcf) / chunk)  # try to merge every
    print("size:", list_length)
    first_create_ = True
    second_create_ = True
    tmp_dirname = mkdtemp(
        dir="/home/kongkitimanonk/SCRATCH_NOBAK/CovSonar1/workdir_covsonar/test-vcf/",
        prefix=".final.sonarCache_",
    )
    # we can tweak performance by using U at Bcftools for piping between bcftools subcommands (future work)
    bar = tqdm(range(list_length), desc="Create Global VCF:")
    merge_type = "b"
    for i in bar:
        _vcfs = " ".join(list_track_vcf[chunk * i : chunk * i + chunk])

        if len(list_track_vcf) == 1:
            tmp_output = list_track_vcf[i].replace(".gz", "")
            continue

        if i == list_length - 1:
            merge_type = "v"
            # print('final merge')

        if first_create_:
            tmp_output = os.path.join(tmp_dirname, "vcf.2")
            cmd = "bcftools merge {} -o {} -O{} --threads {}".format(
                _vcfs, tmp_output, merge_type, num_cores
            )
            with subprocess.Popen(cmd, encoding="utf8", shell=True) as process:
                stdout, stderr = process.communicate(cmd)
            # bgzip(tmp_output)
            # tabix_index(tmp_output)
            bcftool_index(tmp_output)
            first_create_ = False
            second_create_ = True
        elif second_create_:
            _vcfs = _vcfs + " " + os.path.join(tmp_dirname, "vcf.2")
            tmp_output = os.path.join(tmp_dirname, "vcf.3")

            cmd = "bcftools merge {} -o {} -O{} --threads {}".format(
                _vcfs, tmp_output, merge_type, num_cores
            )
            with subprocess.Popen(cmd, encoding="utf8", shell=True) as process:
                stdout, stderr = process.communicate(cmd)
            # bgzip(tmp_output)
            # tabix_index(tmp_output)
            bcftool_index(tmp_output)
            second_create_ = False
        else:
            _vcfs = _vcfs + " " + os.path.join(tmp_dirname, "vcf.3")
            tmp_output = os.path.join(tmp_dirname, "vcf.2")

            cmd = "bcftools merge {} -o {} -O{} --threads {}".format(
                _vcfs, tmp_output, merge_type, num_cores
            )
            with subprocess.Popen(cmd, encoding="utf8", shell=True) as process:
                stdout, stderr = process.communicate(cmd)
            # bgzip( tmp_output)
            # tabix_index(tmp_output)
            bcftool_index(tmp_output)
            second_create_ = True

        if merge_type == "v":
            bgzip(tmp_output)
            tabix_index(tmp_output)
    tmp_output = clean_stranger_things(tmp_output + ".gz", tmp_dirname)
    shutil.copy(tmp_output, global_output + ".gz")
    # sys.exit(" exist.")
    print("Clean workspace ...")
    if os.path.isdir(tmp_dirname):
        shutil.rmtree(tmp_dirname)


def clean_stranger_things(path_to_vcfgz, tmp_dirname):
    print("Clean strange things in vcf ...")
    output_path_file = os.path.join(tmp_dirname, "vcf.final.gz")
    with gzip.open(path_to_vcfgz, "rt") as f:
        with gzip.open(output_path_file, "wt") as output_file:
            for line in f:
                if line.startswith("#"):
                    if "bcftools_mergeCommand" in line:
                        continue
                    else:
                        output_file.write(line)
                else:
                    rows = line.split("\t")
                    ### fix duplicate
                    ID = rows[2]
                    rows[2] = ";".join(set(ID.split(";")))

                    output_file.write("\t".join(rows))
    return output_path_file
