#!/usr/bin/python
# Maintainer: KongkitimanonK
# The method originally came from  https://github.com/cov-lineages/pango-designation.
# We just adapt and change some parts to be used in covsonar, vocal etc.

import importlib.resources
import json
import logging
import os
import shutil
from tempfile import mkdtemp

import pandas as pd
import requests

try:  # noqa: C901
    from pango_aliasor.aliasor import Aliasor


except ModuleNotFoundError:  # pragma: no cover
    logging.WARN(
        "Dependency `pango_aliasor` missing, please install using `pip install pango_aliasor`"
    )
    logging.WARN("Fall back to original Aliasor...")

    class Aliasor:
        def __init__(self, alias_file=None):
            if alias_file is None:
                with importlib.resources.open_text(
                    "pango_designation", "alias_key.json"
                ) as file:
                    file = json.load(file)

            else:
                with open(alias_file) as file:
                    file = json.load(file)

            self.alias_dict = {}
            for column in file.keys():
                if type(file[column]) is list or file[column] == "":
                    self.alias_dict[column] = column
                else:
                    self.alias_dict[column] = file[column]

            self.realias_dict = {v: k for k, v in self.alias_dict.items()}

        def compress(self, name):
            name_split = name.split(".")
            levels = len(name_split) - 1
            num_indirections = (levels - 1) // 3
            if num_indirections <= 0:
                return name
            alias = ".".join(name_split[0 : (3 * num_indirections + 1)])
            ending = ".".join(name_split[(3 * num_indirections + 1) :])
            return self.realias_dict[alias] + "." + ending

        def uncompress(self, name):
            name_split = name.split(".")
            letter = name_split[0]
            try:
                unaliased = self.alias_dict[letter]
            except KeyError:
                return name
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

    def process_lineage_data(self):
        aliasor = Aliasor(self.alias_file)
        df_lineages = pd.read_csv(self.lineage_file)
        lineages = df_lineages.lineage.unique()

        uncompressed_lineages = []
        sorted_lineages = []

        # Calculating parent-child relationship
        cleanedlineages = [x for x in lineages if str(x) != "nan"]
        uncompressed_lineages = list(map(aliasor.uncompress, cleanedlineages))
        uncompressed_lineages.sort(key=sonarLinmgr.lts)
        sorted_lineages = list(map(aliasor.compress, uncompressed_lineages))

        _final_list = []
        for _id in sorted_lineages:
            alias_lineage_char = aliasor.uncompress(_id)
            sub_lineage_list = []
            row_dict = {}
            # print(_id, '=',alias_lineage_char)

            for name_ in uncompressed_lineages:  # fetch all lineage again
                root = ""
                for index, letter in enumerate(name_.split(".")):
                    if index != 0:
                        letter = root + "." + letter
                    root = letter
                    if letter == alias_lineage_char:
                        sub_lineage_list.append(aliasor.compress(name_))
            # remove root lineage
            sub_lineage_list.remove(_id)
            if len(sub_lineage_list) > 0:
                row_dict["lineage"] = _id
                row_dict["sublineage"] = ",".join(sub_lineage_list)
            else:
                row_dict["lineage"] = _id
                row_dict["sublineage"] = "none"
            _final_list.append(row_dict)

        df = pd.DataFrame.from_dict(_final_list, orient="columns")

        return df.sort_values(by=["lineage"])

    def update_lineage_data(self):
        self.download_lineage_data()
        df = self.process_lineage_data()
        return df
