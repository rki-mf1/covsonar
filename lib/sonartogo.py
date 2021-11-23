#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
from lib.sonardb import sonarDBManager
from contextlib import ExitStack
from more_itertools import consecutive_groups, split_when
from tempfile import mkstemp, mkdtemp
import shutil
import pandas as pd
import numpy as np
import subprocess
from os import getpid
from multiprocessing import Pool
import warnings
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning) 

num_partitions = 20 #number of partitions to split dataframe
num_cores = 20 #number of cores on your machine


def create_fix_vcf_header(ref,sample_id):
    header = "##fileformat=VCFv4.2\n##poweredby=CovSonarV1\n##reference="+ref
    format = '\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
    info = '\n##INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">'
    info = info+'\n##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">'

    column = "\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+sample_id+"\n"
    return header+format+info+column

from multiprocessing import Pool


def bgzip(filename):
    """Call bgzip to compress a file."""
    cmd = ['bgzip', '-f', filename]
    with subprocess.Popen(cmd, encoding='utf8', stdout=subprocess.PIPE, stderr=subprocess.STDOUT) as process:
        try:
            stdout, stderr = process.communicate(cmd)

        except subprocess.TimeoutExpired:
            process.kill()
            stdout, stderr = process.communicate()
            raise subprocess.TimeoutExpired( output=stdout, stderr=stderr)
        except Exception:
            process.kill()
            raise
    return

def tabix_index(filename):
    """Call tabix to create an index for a bgzip-compressed file."""
    cmd = ['tabix', '-p', 'vcf', filename+'.gz']
    with subprocess.Popen(cmd, encoding='utf8', stdout=subprocess.PIPE, stderr=subprocess.STDOUT) as process:
        try:
            stdout, stderr = process.communicate(cmd)
        except subprocess.TimeoutExpired:
            process.kill()
            stdout, stderr = process.communicate()
            raise subprocess.TimeoutExpired( output=stdout, stderr=stderr)
        except Exception:
            process.kill()
            raise
    return        

def create_vcf(rows_grouped, tmp_dirname, refdescr):
    process_id =str(getpid())
    # print(process_id+" Start")
    # iterate over each group
    for group_name, df_group in rows_grouped:
        #print("Create VCF file:",group_name)
        vcf_filename =group_name+'.vcf'
        full_path = os.path.join(tmp_dirname,vcf_filename)
        with open(full_path, 'w') as f:
            f.write(create_fix_vcf_header(refdescr,group_name))
            df_group = df_group.sort_values(by='start', ascending=True)
            # replace null to .
            df_group['alt'] = df_group['alt'].replace('', np.nan)
            df_group = df_group .dropna(axis=0, subset=['alt'])
            for index, row in df_group.iterrows():
                id = row['ref']+str(row['start'])+row['alt']
                f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(row['CHROM'], row['start'], id, 
                                        row['ref'], row['alt'],row['QUAL'],row['FILTER'],
                                        row['INFO'],row['FORMAT'],"1"))
        bgzip(full_path)
        tabix_index(full_path) 

    # print(process_id+" Finish")



def parallelize_dataframe(df, tmp_dir,refdescr, func):
    _tmp_lis = np.array_split(df, num_partitions)

    zip_items = [(d,tmp_dir,refdescr) for d in _tmp_lis] # same order
    #with Pool(processes=num_cores) as pool:
    #    res = pool.starmap(func, zip_items)
    pool = Pool(num_cores)
    pool.starmap(func, zip_items)

	# finish all tasks
    pool.close()
    pool.join()

def export2VCF(
    db_path,
    include_acc,
    include_dates,
    output,
    refdescr="No ref is provided"
    ):

    with ExitStack() as stack:
        dbm = stack.enter_context(sonarDBManager(db_path))
        #print(dbm.cursor)
        # for i in tqdm(range(len(fnames)), desc = msg, disable=disable_progressbar):
        # ....self.import_genome(**self.process_fasta(fnames[i]), dbm=dbm)
        # creating where condition

        where_clause = []
        where_vals = []

        if include_acc:
            where_clause.append(dbm.get_metadata_in_condition('accession'
                                , *include_acc))
            where_vals.extend(include_acc)
        if include_dates:
            where_clause.append(dbm.get_metadata_date_condition('date',
                                *include_dates))
        #print(where_clause)
        fields = 'accession, start, end, alt, ref '
        if where_clause:
            sql =  "SELECT " + fields + " FROM dna_view WHERE " + " AND ".join(where_clause) + ";"
        else:
            sql = "SELECT " + fields + " FROM dna_view;"

        # rows = dbm.execute(sql, [include_acc]).fetchall()
        rows = pd.read_sql(sql,
                   dbm.connection,params=where_vals)
        #print(rows['alt'].unique())
    if not rows.empty:
        tmp_dirname = mkdtemp(prefix=".sonarCache_")
        
        # vcf_path=os.path.join(tmp_dirname,)

        # create fasta_id
        chrom_id = refdescr.split()[0].replace(">", "")
        rows['CHROM'] = chrom_id
        rows['QUAL'] = '.'
        rows['FILTER'] = '.'
        rows['INFO'] = 'AC=1;AN=1'
        rows['FORMAT'] = 'GT'
        rows_grouped = rows.groupby('accession')
        
        # split data
        parallelize_dataframe(rows_grouped, tmp_dirname, refdescr, create_vcf)
        
        # bundle all vcf together 
        print("Finish all tasks, start to merge vcfs")
        cmd = "bcftools merge {}/*.vcf.gz -o {} -O v --threads 20".format(tmp_dirname, output)
        print(cmd)
        with subprocess.Popen(cmd, encoding='utf8', shell=True, stdout=subprocess.PIPE) as process:
            try:
                stdout, stderr = process.communicate(cmd)
            except subprocess.TimeoutExpired:
                process.kill()
                stdout, stderr = process.communicate()
                raise subprocess.TimeoutExpired( output=stdout, stderr=stderr)
            except Exception:
                process.kill()
                raise
        if os.path.isdir(tmp_dirname):
            shutil.rmtree(tmp_dirname)
    else:
        print("Noting to export")

        return None
