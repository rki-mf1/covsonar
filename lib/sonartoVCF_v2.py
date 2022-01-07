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
import math
from tqdm import tqdm
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning) 
import traceback

#num_partitions = 20 # number of partitions to split dataframe
#num_cores = 20 # number of cores on your machine

def create_fix_vcf_header(ref):
    header = "##fileformat=VCFv4.2\n##CreatedBy=covSonarV1.1.2\n##reference="+ref
    format = '\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
    info = '\n##INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">'
    info = info+'\n##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">\n'
    note =  "##Note='Currently ignore INDEL'\n"
    # column = "\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+sample_id+"\n"
    return header+format+info+note

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

def bcftool_index(filename):
    """Call tabix to create an index for a bgzip-compressed file."""
    cmd = ['bcftools', 'index', filename]
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


def calculate_AC_AN(final_df):
    # calculate AC AN
    # order-preserving index between POS and INFO AC
    # e.g. a.POS;b.POS a.AC,b.AC
    for row in final_df.itertuples():
        # print('POS '+str(row.POS))
        unique, counts = np.unique(np.asarray(row[10:]), return_counts=True) # row[10:] means we start from sample ID column
        # for unique, counts in zip(unique, counts):
        AN=0
        AC=''
        for idx, val in enumerate(unique):
            if(val == '.'): # ignore it 
                continue
            else:
                _AC = counts[idx]
                AN = AN +_AC
                AC = str(_AC) if not AC else AC+','+str(_AC)
        # print('AN='+str(AN)+';AC='+AC)
        final_df.at[row.POS, 'INFO'] = 'AN='+str(AN)+';AC='+AC
    return final_df


def create_vcf(rows_grouped, tmp_dirname, refdescr):

    process_id =str(getpid())
    # print(process_id+" Start")
    # iterate over each group
    final_df = pd.DataFrame()
    final_df.index = np.arange(1, 29904)
    final_df['#CHROM'] = refdescr.split()[0].replace(">", "")
    final_df['POS'] = np.arange(1, 29904)
    final_df['ID'] = '.'
    final_df['REF'] = '.'
    final_df['ALT'] = '.'
    final_df['FILTER'] = '.'
    final_df['QUAL'] = '.'
    final_df['INFO'] = '.'
    final_df['FORMAT'] = 'GT'
    final_df['POS'] = final_df['POS'].astype('int16') # int16 bit max value 32767
    vcf_filename =process_id+'.vcf'
    full_path = os.path.join(tmp_dirname,vcf_filename)
    with open(full_path, 'w') as f:
        # print("Create VCF file:",full_path)
        f.write(create_fix_vcf_header(refdescr))
        for group_name, df_ in tqdm(rows_grouped, mininterval=0.5):
            final_df[group_name] = "."
            for row in df_.itertuples():
                try:
                    if( getattr(row, 'start') < 1 or getattr(row, 'start') > 29903):
                        continue
                    selected_final_row = final_df.loc[getattr(row, 'start')]
                    index_start_postion = getattr(row, 'start')
                    id = getattr(row, 'ref')+str(getattr(row, 'start'))+getattr(row, 'alt')

                    if(selected_final_row.ID=='.'):
                        final_df.at[index_start_postion, 'ID'] = id
                        final_df.at[index_start_postion, group_name] = 1
                        final_df.at[index_start_postion, 'REF'] = getattr(row, 'ref')
                        final_df.at[index_start_postion, 'ALT'] = getattr(row, 'alt')
                    else: # update 
                        # check ref and alt 
                        if(selected_final_row.ID==id): # only one ID exists
                            final_df.at[index_start_postion, group_name] = 1
                        else:
                            splited_final_id_list = selected_final_row.ID.split(";")
                            totel_len = len(splited_final_id_list)
                            for new_GT, splited_final_id in enumerate(splited_final_id_list, start=1):
                                if(splited_final_id==id): # find the exist one
                                    final_df.at[index_start_postion, group_name] = new_GT
                                    break
                                elif(totel_len == new_GT): # cannot find the same id ,so we append the new one to the string
                                    final_df.at[index_start_postion, 'ID'] = final_df.at[index_start_postion,'ID']+ ';'+id
                                    # with new GT number
                                    final_df.at[index_start_postion, group_name] = new_GT + 1
                                    # appends new alt 
                                    final_df.at[index_start_postion, 'ALT'] = final_df.at[index_start_postion,'ALT']+ ','+getattr(row, 'alt')
                except Exception as e:
                    print("An exception occurred at...") 
                    print(group_name)
                    print(row)
                    print(traceback.format_exc())
                    continue

        final_df = calculate_AC_AN(final_df)
        final_df = final_df.drop(final_df[final_df.ID=='.'].index)   
        final_df.to_csv(f, sep='\t', encoding='utf-8', index=False)
    bgzip(full_path)
    tabix_index(full_path) 
    return full_path+'.gz'
    # print(process_id+" Finish")



def parallelize_dataframe(df, tmp_dirname, num_cores,refdescr, func):
    _tmp_lis = np.array_split(df, num_cores)
    counter=0
    zip_items = [(_tmp_lis[i],tmp_dirname,refdescr) for i in range(len(_tmp_lis))] # same order
    #with Pool(processes=num_cores) as pool:
    #    res = pool.starmap(func, zip_items)
    pool = Pool(num_cores)
    full_paht_list = pool.starmap(func, zip_items)

	# finish all tasks
    pool.close()
    pool.join()
    #print("Tmp result: ", full_paht_list)
    return full_paht_list

def export2VCF(
    db_path,
    include_acc,
    include_dates,
    output,
    num_cores,
    refdescr="No ref is provided"
    ):
    print('----- You are using sonartoVCF_V2 --------')
    print('WARNING: the function is still experimental/not fully implemented.')
    print('Prepare export2VCF workspace for',num_cores,'cpu')
    with ExitStack() as stack:
        dbm = stack.enter_context(sonarDBManager(db_path))

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

        # print("query: " + sql)
        # print("vals: ", where_vals)
        ##############################
        print('Start Bigquery...')
        rows = pd.read_sql(sql,
                   dbm.connection,params=where_vals)
        print('Return:', len(rows), ' records')
        track_vcf = []
        count = 0
        if not rows.empty:
            tmp_dirname = mkdtemp( prefix=".sonarCache_")
            # vcf_path=os.path.join(tmp_dirname,)

            # create fasta_id
  
            #rows['CHROM'] = chrom_id
            #rows['QUAL'] = '.'
            #rows['FILTER'] = '.'
            #rows['INFO'] = '.'
            #rows['FORMAT'] = 'GT'
            # POS or start position: The reference position, with the 1st base is position 1 not 0 , but in covsonar use 0 as the 1st position
            # so we should + 1
            # http://samtools.github.io/hts-specs/VCFv4.2.pdf
            rows['start'] = rows['start'] + 1
            rows['end'] = rows['end']+1
            rows['alt'] = rows['alt'].replace('', np.nan) # remove Deletion
            # rows['start'] = rows['start'].replace('', np.nan) # remove Insertion
            rows = rows.dropna(axis=0, subset=['alt'])
            rows_grouped = rows.groupby('accession')
            print('With :', len(rows_grouped), ' accessions')
            # split data and write each ACC into individual VCF file.
            print('Start Divide and Conquer ...')
            track_vcf = parallelize_dataframe(rows_grouped, tmp_dirname, num_cores, refdescr, create_vcf)
            
            # bundle all vcf together 
            print('Integrate all VCFs ...')
            divide_merge_vcf(track_vcf, output, num_cores)


            if os.path.isdir(tmp_dirname):
                shutil.rmtree(tmp_dirname)

            print("Finish! compress final result (gz):")
                
def divide_merge_vcf(list_track_vcf, global_output, num_cores):
    chunk=500
    list_length = math.ceil(len(list_track_vcf)/chunk) # try to merge every
    print('size:', list_length)
    first_create_ = True 
    second_create_ = True
    tmp_dirname = mkdtemp( prefix=".final.sonarCache_")
    # we can tweak performance by using U at Bcftools for piping between bcftools subcommands (future work)
    bar = tqdm(range(list_length), desc="Create Global VCF:")
    merge_type='b'
    for i in bar:
        _vcfs = " ".join(list_track_vcf[chunk*i:chunk*i+chunk])


        if(len(list_track_vcf)==1):
            tmp_output = list_track_vcf[i].replace('.gz', '')
            continue

        if(i == list_length-1):
            merge_type='v'
            #print('final merge')


        if first_create_:
            tmp_output = os.path.join(tmp_dirname,'vcf.2' )
            cmd = "bcftools merge {} -o {} -O{} --threads {}".format(_vcfs,tmp_output, merge_type, num_cores)    
            with subprocess.Popen(cmd, encoding='utf8', shell=True) as process:
                stdout, stderr = process.communicate(cmd)
            #bgzip(tmp_output)
            #tabix_index(tmp_output)  
            bcftool_index(tmp_output)
            first_create_ = False
            second_create_ = True
            third_create_ = True
        elif second_create_:
            _vcfs = _vcfs +' '+ os.path.join(tmp_dirname,'vcf.2' )
            tmp_output = os.path.join(tmp_dirname,'vcf.3' )

            cmd = "bcftools merge {} -o {} -O{} --threads {}".format(_vcfs,  tmp_output, merge_type, num_cores)
            with subprocess.Popen(cmd, encoding='utf8', shell=True) as process:
                stdout, stderr = process.communicate(cmd)
            #bgzip(tmp_output)
            #tabix_index(tmp_output) 
            bcftool_index(tmp_output)
            second_create_ = False
            third_create_ = True
        else:
            _vcfs = _vcfs +' '+ os.path.join(tmp_dirname,'vcf.3' )
            tmp_output = os.path.join(tmp_dirname,'vcf.2' )

            cmd = "bcftools merge {} -o {} -O{} --threads {}".format(_vcfs, tmp_output, merge_type, num_cores)  
            with subprocess.Popen(cmd, encoding='utf8', shell=True) as process:
                stdout, stderr = process.communicate(cmd)
            #bgzip( tmp_output)
            #tabix_index(tmp_output) 
            bcftool_index(tmp_output)
            second_create_ = True
            third_create_ = False

        if(merge_type =='v'):
            bgzip(tmp_output)
            tabix_index(tmp_output)  

    shutil.copy( tmp_output + '.gz', global_output+ '.gz')
    
    print('Clean workspace ...')
    if os.path.isdir(tmp_dirname):
        shutil.rmtree(tmp_dirname)   
