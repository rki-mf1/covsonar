#!/usr/bin/python
# Maintainer: KongkitimanonK
# The method originally came from  https://github.com/cov-lineages/pango-designation.
# We just adapt and change some parts to be used in covsonar, vocal etc.

import importlib.resources
import json
import os
import shutil
from tempfile import mkdtemp
from typing import Optional

import pandas as pd
import requests


class Aliasor:
    """
    Aliasor class is used to handle aliases.

    Attributes:
        alias_dict (dict): Alias dictionary.
        realias_dict (dict): Reverse alias dictionary.
    """

    def __init__(self, alias_file: Optional[str] = None):
        """
        Aliasor Constructor.

        Args:
            alias_file (str): File containing the alias information.
        """

        # Load the alias file
        if alias_file is None:
            with importlib.resources.open_text(
                "pango_designation", "alias_key.json"
            ) as file:
                alias_data = json.load(file)
        else:
            with open(alias_file) as file:
                alias_data = json.load(file)

        # Create the alias and realias dictionaries
        self.alias_dict = {
            column: (
                alias_data[column]
                if type(alias_data[column]) is not list and alias_data[column] != ""
                else column
            )
            for column in alias_data.keys()
        }
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
    """
    sonarLinmgr class is used to manage the lineages.

    Attributes:
        _tmpdir (str): Temporary directory.
        _linurl (str): URL of the lineages file.
        _aliurl (str): URL of the alias file.
        lineage_file (str): Local path of the lineages file.
        alias_file (str): Local path of the alias file.
    """

    def __init__(self, tmpdir: Optional[str] = None):
        """
        sonarLinmgr Constructor.

        Args:
            tmpdir (str): Temporary directory.
        """
        self._tmpdir = mkdtemp(prefix=".tmp_sonarLinmgr_")
        self._linurl = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineages.csv"
        self._aliurl = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json"
        self.lineage_file = os.path.join(self._tmpdir, "lineages.csv")
        self.alias_file = os.path.join(self._tmpdir, "alias.json")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        shutil.rmtree(self._tmpdir)

    def download_lineage_data(self):
        """
        Download lineage and alias data.
        """
        # Download and write the lineage file
        url_content = requests.get(self._linurl)
        with open(self.lineage_file, "wb") as handle:
            handle.write(url_content.content)

        # Download and write the alias file
        items = requests.get(self._aliurl)
        with open(self.alias_file, "w") as handle:
            json.dump(items.json(), handle)

    @staticmethod
    def lts(lineage: str) -> str:
        """
        Convert lineage into a sortable format.

        Args:
            lineage (str): The lineage.

        Returns:
            str: The lineage in sortable format.
        """
        items = []
        for item in lineage.split("."):
            item_string = str(item)
            items.append((5 - len(item)) * "0" + item_string)
        return "".join(items)

    def process_lineage_data(self) -> pd.DataFrame:
        """
        Process the lineage data.

        Returns:
            pd.DataFrame: The dataframe with processed data.
        """
        aliasor = Aliasor(self.alias_file)
        df_lineages = pd.read_csv(self.lineage_file, usecols=[0, 1])
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

    def update_lineage_data(self) -> pd.DataFrame:
        """
        Update the lineage data.

        Returns:
            pd.DataFrame: The dataframe with updated data.
        """
        self.download_lineage_data()
        df = self.process_lineage_data()
        return df
