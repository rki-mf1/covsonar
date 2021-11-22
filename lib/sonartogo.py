#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
import re
import sys
import argparse
from lib.sonardb import sonarDBManager
from contextlib import ExitStack
from more_itertools import consecutive_groups, split_when
from tempfile import mkstemp, mkdtemp
import vcfpy
import shutil
import pandas as pd

def create_fix_vcf_header(ref,sample_id):
    header = "##fileformat=VCFv4.2\n##poweredby=covsonar\n##reference="+ref
    format = '\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
    info = '\n##INFO=<ID=AC,Number=1,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">'
    info = info+'\n##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">'

    column = "\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+sample_id+"\n"
    return header+format+info+column

def export2VCF(
    db_path,
    include_acc,
    include_dates,
    export_all=False,
    refdescr="No ref is provided"
    ):

    with ExitStack() as stack:
        dbm = stack.enter_context(sonarDBManager(db_path))
        print(dbm.cursor)
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
        print(where_clause)
        fields = 'accession, start, end, alt, ref '
        if where_clause:
            sql =  "SELECT " + fields + " FROM dna_view WHERE " + " AND ".join(where_clause) + ";"
        else:
            sql = "SELECT " + fields + " FROM dna_view;"

        # rows = dbm.execute(sql, [include_acc]).fetchall()
        rows = pd.read_sql(sql,
                   dbm.connection,params=where_vals)
        print(rows)
    if not rows.empty:
        #tmp_dirname = mkdtemp(prefix=".sonarCache_")
        tmp_dirname="../workdir_covsonar/vcf-test/"
        # vcf_path=os.path.join(tmp_dirname,)

        # create fasta_id
        chrom_id = refdescr.split()[0].replace(">", "")
        rows['CHROM'] = chrom_id
        rows['QUAL'] = '.'
        rows['FILTER'] = '.'
        rows['INFO'] = 'AC=1;AN=1'
        rows['FORMAT'] = 'GT'
        rows_grouped = rows.groupby('accession')

        # iterate over each group
        for group_name, df_group in rows_grouped:
            print("Create VCF file:",group_name)
            vcf_filename =group_name+'.vcf'
            with open(tmp_dirname+vcf_filename, 'w') as f:
                f.write(create_fix_vcf_header(refdescr,group_name))
                df_group = df_group.sort_values(by='start', ascending=True)
                for index, row in df_group.iterrows():
                    id = row['ref']+str(row['start'])+row['alt']
                    f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(row['CHROM'], row['start'], id, 
                                            row['ref'], row['alt'],row['QUAL'],row['FILTER'],
                                            row['INFO'],row['FORMAT'],"1"))
                    # Open input, add FILTER header, and open output file
                    
        # bundle all vcf together 
        #if os.path.isdir(tmp_dirname):
        #    shutil.rmtree(tmp_dirname
    else:
        print("Noting to export")

        return None
