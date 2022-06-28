#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
from random import sample
from sonar.dbm import sonarDBManager
from contextlib import ExitStack
from tempfile import mkstemp, mkdtemp
import shutil
import pandas as pd
import numpy as np
import subprocess
from os import getpid
from multiprocessing import Pool
import warnings
import math
from tqdm import tqdm
import gzip
from multiprocessing import Pool

warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)


def create_fix_vcf_header(ref, sample_id):
    header = "##fileformat=VCFv4.2\n##poweredby=CovSonarV2\n##reference=" + ref
    format = '\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
    info = '\n##INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">'
    info = (
        info
        + '\n##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">\n'
    )
    note = "##Note_1='Currently, Deletion site is represented as ' ' (one whitespace)\n"
    column = (
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_id + "\n"
    )
    return header + format + info + note + column


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


def create_vcf(rows_grouped, tmp_dirname, refdescr, _pos):
    process_id = str(getpid())
    # print(process_id+" Start")
    # iterate over each group
    # position=_pos,bar_format='{l_bar}{bar:10}{r_bar}{bar:-2b}'
    for group_name, df_group in tqdm(rows_grouped, mininterval=0.5):
        # print("Create VCF file:",group_name)
        vcf_filename = group_name + ".vcf"
        full_path = os.path.join(tmp_dirname, vcf_filename)

        with open(full_path, "w") as f:
            f.write(create_fix_vcf_header(refdescr, group_name))
            df_group = df_group.sort_values(by="start", ascending=True)
            # replace null to .
            # df_group['ref'] = df_group['ref'].replace('', '.') # for insertion
            # df_group["alt"] = df_group["alt"].replace(" ", ".")  # for deletion
            # df_group["alt"] = df_group["alt"].replace(" ", "")  # remove Deletion
            # df_group['alt'] = df_group['alt'].replace('', np.nan) # remove Deletion
            # df_group = df_group.dropna(axis=0, subset=["alt"])  # remove Deletion
            for index, row in df_group.iterrows():
                id = row["ref"] + str(row["start"]) + row["alt"]

                f.write(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        row["CHROM"],
                        row["start"],
                        id,
                        row["ref"],
                        row["alt"],
                        row["QUAL"],
                        row["FILTER"],
                        row["INFO"],
                        row["FORMAT"],
                        "1",
                    )
                )
        bgzip(full_path)
        tabix_index(full_path)

    # print(process_id+" Finish")


def parallelize_dataframe(df, tmp_dir, num_cores, refdescr, func):
    _tmp_lis = np.array_split(df, num_cores)
    counter = 0
    zip_items = [
        (_tmp_lis[i], tmp_dir, refdescr, i + 2) for i in range(len(_tmp_lis))
    ]  # same order
    # with Pool(processes=num_cores) as pool:
    #    res = pool.starmap(func, zip_items)
    pool = Pool(num_cores)
    pool.starmap(func, zip_items)

    # finish all tasks
    pool.close()
    pool.join()


def export2VCF(db_path, output, num_cores=4, properties=None, acc=None):
    print("----- You are using covSonar-to-VCF V1 --------")
    print("Prepare export2VCF workspace for", num_cores, "cpu")
    with ExitStack() as stack:
        dbm = stack.enter_context(sonarDBManager(db_path))
        fields = '  "sample.name" as accession, "reference.accession" as ref_name , "variant.start" as start, "variant.end" as end, "variant.alt" as alt, "variant.ref" as ref '

        # collecting sqls for metadata-based filtering
        if properties or acc:
            if acc:
                _tmp_sample = dbm.get_sample_ids(acc)

            if properties:
                print("Calculate Properties...")
                print(properties)
                _tmp_prop = dbm.query_metadata_forVCF(properties)
                
            if properties and acc:
                sample_id_list = [value for value in _tmp_sample if value in _tmp_prop]
            else:
                if properties:
                    sample_id_list = _tmp_prop
                elif acc:
                    sample_id_list = _tmp_sample
                else:
                    print("ERROR: Something went wrong...")
                    raise
            sql = (
                "SELECT "
                + fields
                + ' FROM variantView WHERE "element.type"="source" '
                + ' AND "sample.id" IN ({0})'.format(", ".join(sample_id_list))
                + ";"
            )
        else:
            print("Warning: you query the whole database...")
            sql = (
                "SELECT "
                + fields
                + ' FROM variantView  WHERE "element.type"="source"  ;'
            )
        # print(sql)
        ##############################
        print("Start Bigquery...")
        rows = pd.read_sql(sql, dbm.connection)
        print("Return:", len(rows), " records")
        track_vcf = []
        count = 0
        if not rows.empty:
            tmp_dirname = mkdtemp(prefix=".sonarCache_")
            # vcf_path=os.path.join(tmp_dirname,)
            # create fasta_id
            if len(rows.groupby("ref_name").size()) == 1:
                refdescr = rows.iloc[0]["ref_name"]
            else:
                print("Error: Found reference more than 1")
                refdescr = "No ref is provided"
                raise

            chrom_id = refdescr.split()[0].replace(">", "")
            rows["CHROM"] = chrom_id
            rows["QUAL"] = "."
            rows["FILTER"] = "."
            rows["INFO"] = "AC=1;AN=1"
            rows["FORMAT"] = "GT"
            # POS or start position: The reference position, with the 1st base is position 1 not 0 , but in covsonar use 0 as the 1st position
            # so we should + 1
            # http://samtools.github.io/hts-specs/VCFv4.2.pdf
            rows["start"] = rows["start"] + 1
            rows_grouped = rows.groupby("accession")
            for group_name, df_group in rows_grouped:
                group_name = os.path.join(tmp_dirname, group_name + ".vcf.gz")
                full_path = os.path.join(tmp_dirname, group_name)
                track_vcf.append(full_path)
            print("With :", len(track_vcf), " accessions")
            # split data and write each ACC into individual VCF file.
            print("Start Divide and Conquer ...")
            parallelize_dataframe(
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
    tmp_dirname = mkdtemp(prefix=".final.sonarCache_")
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
            third_create_ = True
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
            third_create_ = True
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
            third_create_ = False

        if merge_type == "v":
            bgzip(tmp_output)
            tabix_index(tmp_output)
    tmp_output = clean_stranger_things(tmp_output + ".gz", tmp_dirname)
    shutil.copy(tmp_output, global_output + ".gz")

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
