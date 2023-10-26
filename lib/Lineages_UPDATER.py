#!/usr/bin/python
# Maintainer: KongkitimanonK
# The method originally came from
# https://github.com/cov-lineages/pango-designation.
# We just adapt and change some parts to be used in covsonar, vocal etc.
import json
import os
import sys

import pandas as pd
import requests

# warnings.simplefilter(action='ignore', category=FutureWarning)

# Due to the https://github.com/cov-lineages/pango-designation/issues/853
# Now the new repo of Aliasor is here https://github.com/corneliusroemer/pango_aliasor
# try:
#    from pango_aliasor.aliasor import Aliasor
# except ModuleNotFoundError as e:
#    print("Dependency `pango_aliasor` missing, please install using `pip install pango_aliasor`")
#    raise e


# In covSonar 1 we will copy the class from https://github.com/corneliusroemer/pango_aliasor/blob/main/src/pango_aliasor/aliasor.py
# but for V.2 we will switch using the package to "pip install pango_aliasor"
class Aliasor:
    def __init__(self, alias_file=None):
        import json

        if alias_file is None:
            import importlib.resources

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
        # A function that prints the output to the screen.
        if pd.isna(name):
            return ""
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


def lts(lineage):
    items = []
    for item in lineage.split("."):
        item_string = str(item)
        items.append((5 - len(item)) * "0" + item_string)
    return "".join(items)


def download_source(tmp_dir):
    lineages_url = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineages.csv"
    alias_key_url = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json"
    lineag = os.path.join(tmp_dir, "lineags.csv")
    alias_key = os.path.join(tmp_dir, "alias_key.json")
    print("Download lineages", file=sys.stderr)
    url_content = requests.get(lineages_url).content
    csv_file = open(lineag, "wb")
    csv_file.write(url_content)
    csv_file.close()
    print("Download alias_key", file=sys.stderr)
    items = requests.get(alias_key_url)
    data = items.json()
    with open(alias_key, "w") as f:
        json.dump(data, f)
    return alias_key, lineag


def process_lineage(alias_key_path, lineages_path, output):
    """
    # handle duplicate values
    with open(alias_key_path) as f:
                # load json objects to dictionaries
        data_dict = json.load(f)

    for k, v in data_dict.items():
        if type(v) is list:
            data_dict[k] = list(set(v))
        # rewrite the json
    with open(alias_key_path ,'w') as nf:
        json.dump(data_dict, nf)
    """
    print("Create all lineages", file=sys.stderr)
    aliasor = Aliasor(alias_key_path)
    df_lineages = pd.read_csv(lineages_path)
    lineages = df_lineages.lineage.unique()
    # create all available lineages.
    uncompressed_lineages = []
    sorted_lineages = []

    uncompressed_lineages = list(map(aliasor.uncompress, lineages))
    uncompressed_lineages.sort(key=lts)

    sorted_lineages = list(map(aliasor.compress, uncompressed_lineages))
    print("Calculate parent-child relationship", file=sys.stderr)
    # -- Fill the sub child --
    # the current method is working, but it is not good in term of performance
    # the algoritm use 3 loops
    _final_list = []
    for _id in sorted_lineages:
        alias_lineage_char = aliasor.uncompress(_id)
        sub_lineage_list = []
        row_dict = {}
        # print(_id, '=',alias_lineage_char)  # check alias name
        for name_ in uncompressed_lineages:  # fetch all lineage again
            for index, letter in enumerate(name_.split(".")):
                if index != 0:
                    letter = root + "." + letter
                    root = letter
                else:
                    root = letter
                if letter == alias_lineage_char:
                    sub_lineage_list.append(aliasor.compress(name_))

        sub_lineage_list.remove(
            aliasor.compress(alias_lineage_char)
        )  # remove root lineage
        if len(sub_lineage_list) > 0:
            row_dict["lineage"] = _id
            row_dict["sublineage"] = ",".join(sub_lineage_list)
        else:
            row_dict["lineage"] = _id
            row_dict["sublineage"] = "none"
        _final_list.append(row_dict)
    print("Write output:", output, file=sys.stderr)
    df = pd.DataFrame.from_dict(_final_list, orient="columns")
    df = df[df.lineage != ""]
    df.sort_values(by=["lineage"]).to_csv(output, sep="\t", index=False)
