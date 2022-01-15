#!/usr/bin/python
# Maintainer: KongkitimanonK
# Script version 1 (copy from VOCAL)
# The method originally came from  https://github.com/cov-lineages/pango-designation.
# We just adapt and change some parts to be used in covsonar, vocal etc.
import os 
import pandas as pd
from tempfile import mkstemp, mkdtemp
import json
import requests



class Aliasor:
    def __init__(self, alias_file):

        aliases = pd.read_json(alias_file)

        self.alias_dict = {}
        for column in aliases.columns:
            if column.startswith('X'):
                self.alias_dict[column] = column
            else:
                self.alias_dict[column] = aliases[column][0]

        self.alias_dict['A'] = 'A'
        self.alias_dict['B'] = 'B'

        self.realias_dict = {v: k for k, v in self.alias_dict.items()}

    def compress(self,name):
        name_split = name.split('.')
        #print(name_split)
        if len(name_split) < 5:
            return name
        letter = self.realias_dict[".".join(name_split[0:4])]
        if len(name_split) == 5:
            #print('len5:'+letter + '.' + name_split[4])
            return letter + '.' + name_split[4]
        else:
            #print('len6:'+letter + '.' + ".".join(name_split[4:]))
            return letter + '.' + ".".join(name_split[4:])

    def uncompress(self,name):
        name_split = name.split('.')
        #print(name_split)
        letter = name_split[0]
        unaliased = self.alias_dict[letter]
        if len(name_split) == 1:
            return name
        if len(name_split) == 2:
            #print('len2:'+unaliased + '.' + name_split[1])
            return unaliased + '.' + name_split[1]
        else:
            #print('len3:'+unaliased + '.' + ".".join(name_split[1:]))
            return unaliased + '.' + ".".join(name_split[1:])
        
        
def lts(lineage):
    items = []
    for item in lineage.split("."):
        item_string = str(item)
        items.append((5-len(item))*"0" + item_string)
    return "".join(items)

def download_source(tmp_dir):
    lineages_url='https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineages.csv'
    alias_key_url='https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json'
    lineag=os.path.join(tmp_dir,'lineags.csv')
    alias_key=os.path.join(tmp_dir,'alias_key.json')
    print('Download lineags')
    url_content = requests.get(lineages_url).content
    csv_file = open(lineag, 'wb')
    csv_file.write(url_content)
    csv_file.close()
    print('Download alias_key')
    items = requests.get(alias_key_url)
    data = items.json()
    with open( alias_key , 'w') as f:
        json.dump(data, f)
    return alias_key,lineag

def process_lineage(alias_key_path, lineages_path, output):
    print('Calculate all lineages')
    aliasor = Aliasor(alias_key_path)
    df_lineages = pd.read_csv(lineages_path)
    lineages = df_lineages.lineage.unique()
    #%%
    uncompressed_lineages = []
    sorted_lineages = []
    print('Calculate parent-child relationship')
    for ch in map(aliasor.uncompress, lineages):
        uncompressed_lineages.append(ch)
    uncompressed_lineages.sort(key=lts)
    for ch in map(aliasor.compress, uncompressed_lineages):
        sorted_lineages.append(ch)
    #%%
    print('To output')
    df = pd.DataFrame()
    for _id in lineages:
        sub_lineage_char = aliasor.realias_dict.get(_id)
        sub_lineage_list = []
        
        for name_ in sorted_lineages:
            letter = name_.split('.')[0]
            if sub_lineage_char == letter:
                sub_lineage_list.append(name_)
        
        sub_lineage_list = list(filter((_id).__ne__, sub_lineage_list))
        if(len(sub_lineage_list)):
            df = df.append({'lineage': _id, 'sublineage': ",".join(sub_lineage_list)}, ignore_index=True)
        else:
            df = df.append({'lineage': _id, 'sublineage': 'none'}, ignore_index=True)
    df.to_csv(output, sep="\t", index=False)

