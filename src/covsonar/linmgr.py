#!/usr/bin/python
# Maintainer: KongkitimanonK
# Script version 1 (copy from VOCAL)
# The method originally came from  https://github.com/cov-lineages/pango-designation.
# We just adapt and change some parts to be used in covsonar, vocal etc.

import json
import os
import shutil
from tempfile import mkdtemp

import pandas as pd
import requests


class Aliasor:
    def __init__(self, alias_file):
        aliases = pd.read_json(alias_file)
        self.alias_dict = {
            x: x if x.startswith("X") else aliases[x][0] for x in aliases.columns
        }
        self.alias_dict["A"] = "A"
        self.alias_dict["B"] = "B"
        self.realias_dict = {v: k for k, v in self.alias_dict.items()}

    def compress(self, name):
        name_split = name.split(".")
        if len(name_split) < 5:
            return name
        letter = self.realias_dict[".".join(name_split[0:4])]
        if len(name_split) == 5:
            return letter + "." + name_split[4]
        else:
            return letter + "." + ".".join(name_split[4:])

    def uncompress(self, name):
        name_split = name.split(".")
        letter = name_split[0]
        unaliased = self.alias_dict[letter]
        if len(name_split) == 1:
            return name
        if len(name_split) == 2:
            return unaliased + "." + name_split[1]
        else:
            return unaliased + "." + ".".join(name_split[1:])


class sonarLinmgr:
    def __init__(self, tmpdir=None):
        self._tmpdir = mkdtemp(prefix=".tmp_sonarLinmgr_")
        self._linurl = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineages.csv"
        self._aliurl = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json"
        self.lineage_file = os.path.join(self._tmpdir, "lineages.csv")
        self.alias_file = os.path.join(self._tmpdir, "alias.json")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        # print('Clean env.')
        shutil.rmtree(self._tmpdir)

    # download
    def download_lineage_data(self):
        url_content = requests.get(self._linurl)
        with open(self.lineage_file, "wb") as handle:
            handle.write(url_content.content)
        items = requests.get(self._aliurl)
        with open(self.alias_file, "w") as handle:
            json.dump(items.json(), handle)

    # processing
    @staticmethod
    def lts(lineage):
        items = []
        for item in lineage.split("."):
            item_string = str(item)
            items.append((5 - len(item)) * "0" + item_string)
        return "".join(items)

    def process_lineage_data(self, output_file):
        # handle duplicate values
        with open(self.alias_file) as f:
            # load json objects to dictionaries
            data_dict = json.load(f)
        for k, v in data_dict.items():
            if type(v) is list:
                data_dict[k] = list(set(v))
        # rewrite the json
        with open(self.alias_file, "w") as nf:
            json.dump(data_dict, nf)

        aliasor = Aliasor(self.alias_file)
        df_lineages = pd.read_csv(self.lineage_file)
        lineages = df_lineages.lineage.unique()

        uncompressed_lineages = []
        sorted_lineages = []

        # Calculating parent-child relationship
        uncompressed_lineages = [x for x in map(aliasor.uncompress, lineages)]
        uncompressed_lineages.sort(key=sonarLinmgr.lts)
        sorted_lineages = [x for x in map(aliasor.compress, uncompressed_lineages)]

        datadict = []
        for _id in lineages:
            sub_lineage_char = aliasor.realias_dict.get(_id)
            sub_lineage_list = [
                x for x in sorted_lineages if x.split(".")[0] == sub_lineage_char
            ]
            sub_lineage_list = list(filter((_id).__ne__, sub_lineage_list))
            if len(sub_lineage_list):
                datadict.append(
                    {"lineage": _id, "sublineage": ",".join(sub_lineage_list)}
                )
            else:
                datadict.append({"lineage": _id, "sublineage": "none"})

        return pd.DataFrame(
            datadict
        )  # pd.DataFrame(datadict).to_csv(output_file, sep="	", index=False)

    def update_lineage_data(self, output_file):
        self.download_lineage_data()
        df = self.process_lineage_data(output_file)
        return df
