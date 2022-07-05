#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

import argparse
from collections import defaultdict
import os
import sys
from tempfile import mkstemp
from textwrap import fill

from mpire import WorkerPool
from tabulate import tabulate
from tqdm import tqdm

from . import sonardb
from .align import sonarAligner
from .basics import sonarBasics
from .cache import sonarCache
from .dbm import sonarDBManager
from .linmgr import sonarLinmgr


class arg_namespace(object):
    pass


def parse_args(args):
    """
    setting and handling command line arguments
    """

    VERSION = sonarBasics.get_version()

    # preparations
    user_namespace = arg_namespace()
    parser = argparse.ArgumentParser(prog="sonar", description="covsonar " + VERSION)
    subparsers = parser.add_subparsers(
        help="detect, store, and screen for mutations in genomic sequences"
    )
    subparsers.dest = "tool"
    subparsers.required = True

    # parser components
    # parser component: db input
    general_parser = argparse.ArgumentParser(add_help=False)
    general_parser.add_argument(
        "--db", metavar="FILE", help="sonar database to use", type=str, required=True
    )
    general_parser.add_argument(
        "--debug",
        help="activate debugging mode showing all sqllite queries on screen",
        action="store_true",
    )

    # parser component: output
    output_parser = argparse.ArgumentParser(add_help=False)
    output_parser.add_argument(
        "-o",
        "--out",
        metavar="FILE",
        help="write output to file",
        type=str,
        default=None,
    )

    # parser component: sample input
    sample_parser = argparse.ArgumentParser(add_help=False)
    sample_parser.add_argument(
        "--sample",
        metavar="STR",
        help="sample accession(s) to consider ",
        type=str,
        nargs="+",
        default=[],
    )
    sample_parser.add_argument(
        "--sample-file",
        metavar="FILE",
        help="file containing sample accession(s) to consider (one per line)",
        type=str,
        nargs="+",
        default=[],
    )

    # parser component: property name
    prop_parser = argparse.ArgumentParser(add_help=False)
    prop_parser.add_argument(
        "--name", metavar="STR", help="name of sample property", type=str, required=True
    )

    # parser component: reference
    ref_parser = argparse.ArgumentParser(add_help=False)
    ref_parser.add_argument(
        "-r",
        "--reference",
        metavar="STR",
        help="reference accession",
        type=str,
        default=None,
    )

    # parser component: threads
    thread_parser = argparse.ArgumentParser(add_help=False)
    thread_parser.add_argument(
        "-t",
        "--threads",
        metavar="INT",
        help="number of threads to use (default: 1)",
        type=int,
        default=1,
    )

    # main parser
    # setup parser
    parser_setup = subparsers.add_parser(
        "setup", parents=[general_parser], help="Setup a new database."
    )
    parser_setup.add_argument(
        "-a",
        "--auto-create",
        help="Auto create important properties",
        action="store_true",
    )
    parser_setup.add_argument(
        "--gbk",
        metavar="FILE",
        help="genbank file of the reference genome (default MN908947.3 is used as reference)",
        type=str,
        default=None,
    )

    # import parser
    parser_import = subparsers.add_parser(
        "import",
        parents=[general_parser, thread_parser],
        help="Import genome sequences and sample information to the database.",
    )
    parser_import.add_argument(
        "--fasta",
        help="fasta file containing genome sequences to import",
        type=str,
        nargs="+",
        default=None,
    )
    parser_import.add_argument(
        "--tsv",
        help="tab-delimited file containing sample properties to import",
        type=str,
        nargs="+",
        default=None,
    )
    parser_import.add_argument(
        "--cols",
        help="define column names for sample properties (if different from property name)",
        type=str,
        nargs="+",
        default=[],
    )
    parser_import.add_argument(
        "--no-autodetect",
        help="do not auto-detect of columns for sample properties based on property names",
        action="store_true",
    )
    parser_import.add_argument(
        "--no-update",
        help="skip samples already existing in the database",
        action="store_true",
    )
    parser_import.add_argument(
        "--cache",
        metavar="DIR",
        help="directory for chaching data (default: a temporary directory is created)",
        type=str,
        default=None,
    )
    parser_import_g1 = parser.add_mutually_exclusive_group()
    parser_import_g1.add_argument(
        "--no-progress",
        "-p",
        help="don't show progress bars while importing",
        action="store_true",
    )

    # view-prop parser
    subparsers.add_parser(
        "list-prop",
        parents=[general_parser],
        help="View sample properties added to the database.",
    )

    # add-prop parser
    parser_addprops = subparsers.add_parser(
        "add-prop",
        parents=[general_parser, prop_parser],
        help="Add a sample property to the database.",
    )
    parser_addprops.add_argument(
        "--descr",
        metavar="STR",
        help="description of the new property",
        type=str,
        required=True,
    )
    parser_addprops.add_argument(
        "--dtype",
        metavar="STR",
        help="data type of the new property",
        type=str,
        choices=["integer", "float", "text", "date", "zip"],
        required=True,
    )
    parser_addprops.add_argument(
        "--qtype",
        metavar="STR",
        help="query type of the new property",
        type=str,
        choices=["numeric", "float", "date", "text", "zip"],
        default=None,
    )
    parser_addprops.add_argument(
        "--default",
        metavar="VAR",
        help="default value of the new property (none by default)",
        type=str,
        default=None,
    )

    # delete-prop parser
    parser_delprops = subparsers.add_parser(
        "delete-prop",
        parents=[general_parser, prop_parser],
        help="Delete a sample property from the database.",
    )
    parser_delprops.add_argument(
        "--force", help="force property to be deleted", action="store_true"
    )

    # match parser
    parser_match = subparsers.add_parser(
        "match",
        parents=[sample_parser, output_parser, general_parser],
        help="get mutations profiles for given accessions.",
    )
    parser_match.add_argument(
        "--profile",
        "-p",
        metavar="STR",
        help="match genomes sharing the given mutation profile",
        type=str,
        action="append",
        nargs="+",
        default=[],
    )
    parser_match_format = parser_match.add_mutually_exclusive_group()
    parser_match_format.add_argument(
        "--count", help="count instead of listing matching genomes", action="store_true"
    )
    parser_match_format.add_argument(
        "--format",
        help="output format (default: tsv)",
        choices=["csv", "tsv", "vcf"],
        default="tsv",
    )
    parser_match.add_argument(
        "--with-sublineage",
        metavar="STR",
        help="recursively get all sublineages from a given lineage. ",
        type=str,
        default=None,
    )

    # delete parser
    parser_delete = subparsers.add_parser(
        "delete",
        parents=[ref_parser, sample_parser, general_parser],
        help="delete one or more samples from the database.",
    )
    parser_delete.add_argument(
        "--aligned",
        help="ise aligned form (deletions indicated by - and insertions by lower-case letters)",
        action="store_true",
    )

    # restore parser
    parser_restore = subparsers.add_parser(
        "restore",
        parents=[ref_parser, sample_parser, general_parser],
        help="restore sequence(s) from the database.",
    )
    parser_restore.add_argument(
        "--aligned",
        help="ise aligned form (deletions indicated by - and insertions by lower-case letters)",
        action="store_true",
    )

    # info parser
    subparsers.add_parser(
        "info", parents=[general_parser], help="show software and database info"
    )

    # optimize parser
    subparsers.add_parser(
        "optimize", parents=[general_parser], help="optimizes the database."
    )

    # dev parser
    subparsers.add_parser("dev", parents=[general_parser])

    # db-upgrade parser
    subparsers.add_parser(
        "db-upgrade",
        parents=[general_parser],
        help="upgrade a database to the latest version",
    )

    # update-lineage-info parser
    subparsers.add_parser(
        "update-lineage-info",
        parents=[general_parser],
        help="download latest lineage information",
    )

    # version parser
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="covsonar " + VERSION,
        help="Show program's version number and exit.",
    )

    # register known arguments
    pargs = parser.parse_known_args(args=args, namespace=user_namespace)

    # register dynamic arguments
    # register additonal arguments from database-specific sample property names when matching genomes
    if pargs[0].tool == "match" and hasattr(pargs[0], "db"):
        with sonarDBManager(pargs[0].db, readonly=True, debug=pargs[0].debug) as dbm:
            for prop in dbm.properties.values():
                if prop["datatype"] == "integer":
                    t = int
                elif prop["datatype"] == "float":
                    t = float
                else:
                    t = str
                if pargs[0].tool == "match":
                    parser_match.add_argument(
                        "--" + prop["name"],
                        type=t,
                        nargs="+",
                        default=argparse.SUPPRESS,
                    )

    # return
    return parser.parse_args(args=args, namespace=user_namespace)


class sonar:
    def __init__(self, db, gff=None, debug=False):
        self.dbfile = db if db else mkstemp()[1]
        self.db = sonardb.sonarDB(self.dbfile)
        self.gff = gff
        self.debug = debug

    def show_system_info(self):
        """
        Prints out the version of the database, the name of the reference genome, the length of the
        reference genome, the names of the annotated proteins, and the translation table used
        """
        print("sonarDB version:       ", self.db.get_version())
        print("reference genome:      ", self.db.refdescr)
        print("reference length:      ", str(len(self.db.refseq)) + "bp")
        print("annotated proteins:    ", ", ".join(self.db.refgffObj.symbols))
        print("used translation table:", self.db.translation_table)

    def show_db_info(self):
        with sonardb.sonarDBManager(self.dbfile, readonly=True) as dbm:
            print("database path:             ", dbm.dbfile)
            print("database version:          ", dbm.get_db_version())
            print("database size:             ", self.get_db_size())
            g = dbm.count_genomes()
            print("genomes:                   ", g)
            print("unique sequences:          ", dbm.count_sequences())
            print("labs:                      ", dbm.count_labs())
            print("earliest genome import:    ", dbm.get_earliest_import())
            print("latest genome import:      ", dbm.get_latest_import())
            print("earliest sampling date:    ", dbm.get_earliest_date())
            print("latest sampling date:      ", dbm.get_latest_date())
            print("metadata:          ")
            fields = sorted(
                [
                    "lab",
                    "source",
                    "collection",
                    "technology",
                    "platform",
                    "chemistry",
                    "software",
                    "software_version",
                    "material",
                    "ct",
                    "gisaid",
                    "ena",
                    "lineage",
                    "zip",
                    "date",
                ]
            )
            maxlen = max([len(x) for x in fields])
            for field in fields:
                if g == 0:
                    c = 0
                    p = 0
                else:
                    c = dbm.count_metadata(field)
                    p = c / g * 100
                spacer = " " * (maxlen - len(field))
                print("   " + field + " information:" + spacer, f"{c} ({p:.{2}f}%)")


def check_file(fname):
    if not os.path.isfile(fname):
        sys.exit("iput error: " + fname + " is not a valid file.")


def main(args):  # noqa: C901
    # process arguments
    # args = parse_args()

    if hasattr(args, "db") and args.db:
        if args.tool != "setup" and args.db is not None and not os.path.isfile(args.db):
            sys.exit("input error: database does not exist.")
            # check_db_compatibility

    # set debugging mode
    if hasattr(args, "debug") and args.debug:
        debug = True
    else:
        debug = False

    # tool procedures
    # setup, db-upgrade
    if args.tool == "setup":
        if args.gbk:
            check_file(args.gbk)
        sonarBasics.setup_db(
            args.db, args.auto_create, reference_gb=args.gbk, debug=debug
        )
    elif args.tool == "db-upgrade":
        input("Warning: Backup db file before upgrading, Press Enter to continue...")
        sonarDBManager.upgrade_db(args.db)
    else:
        with sonarDBManager(args.db, readonly=True) as dbm:
            dbm.check_db_compatibility()
    # other than the above
    # import
    if args.tool == "import":
        if args.no_update:
            update = False
            print("import mode: skipping existing samples")
            if args.fasta is None:
                print("Nothing to import.")
                exit(0)
        else:
            update = True
            print("import mode: updating existing samples")

        if args.tsv is None and args.fasta is None:
            print("Nothing to import.")
            exit(0)

        # prop handling
        with sonarDBManager(args.db, readonly=True) as dbm:
            db_properties = set(dbm.properties.keys())
            db_properties.add("sample")

        colnames = {} if args.no_autodetect else {x: x for x in db_properties}
        if args.cols is not None:
            for x in args.cols:
                if x.count("=") != 1:
                    sys.exit(
                        "input error: "
                        + x
                        + " is not a valid sample property assignment."
                    )
                k, v = x.split("=")
                if k not in db_properties:
                    sys.exit(
                        "input error: sample property "
                        + k
                        + " is unknown to the selected database. Use list-props to see all valid properties."
                    )
                colnames[k] = v

            if "sample" not in colnames:
                sys.exit("input error: a sample column has to be assigned.")

        properties = defaultdict(dict)

        if args.tsv is not None:
            for tsv in args.tsv:
                with open(tsv, "r") as handle, tqdm(
                    desc="processing " + tsv + "...",
                    total=os.path.getsize(tsv),
                    unit="bytes",
                    unit_scale=True,
                    bar_format="{desc} {percentage:3.0f}% [{n_fmt}/{total_fmt}, {elapsed}<{remaining}, {rate_fmt}{postfix}]",
                    disable=args.no_progress,
                ) as pbar:
                    line = handle.readline()
                    pbar.update(len(line))
                    fields = line.strip("\r\n").split("\t")
                    tsv_cols = {}
                    print()
                    for x in sorted(colnames.keys()):
                        c = fields.count(colnames[x])
                        if c == 1:
                            tsv_cols[x] = fields.index(colnames[x])
                            print("  " + x + " <- " + colnames[x])
                        elif c > 1:
                            sys.exit(
                                "error: " + colnames[x] + " is not an unique column."
                            )
                    if "sample" not in tsv_cols:
                        sys.exit(
                            "error: tsv file does not contain required sample column."
                        )
                    elif len(tsv_cols) == 1:
                        sys.exit(
                            "input error: tsv does not provide any informative column."
                        )
                    for line in handle:
                        pbar.update(len(line))
                        fields = line.strip("\r\n").split("\t")
                        sample = fields[tsv_cols["sample"]]
                        for x in tsv_cols:
                            if x == "sample":
                                continue
                            properties[sample][x] = fields[tsv_cols[x]]

        # setup cache
        temp = False if args.cache else True
        cache = sonarCache(
            args.db,
            outdir=args.cache,
            logfile="import.log",
            allow_updates=update,
            temp=temp,
            debug=args.debug,
            disable_progress=args.no_progress,
        )

        # importing sequences
        if args.fasta is not None:
            cache.add_fasta(*args.fasta, propdict=properties)

            aligner = sonarAligner()
            l = len(cache._samplefiles_to_profile)
            with WorkerPool(n_jobs=args.threads, start_method="fork") as pool, tqdm(
                desc="profiling sequences...",
                total=l,
                unit="seqs",
                bar_format="{desc} {percentage:3.0f}% [{n_fmt}/{total_fmt}, {elapsed}<{remaining}, {rate_fmt}{postfix}]",
                disable=args.no_progress,
            ) as pbar:
                for _ in pool.imap_unordered(
                    aligner.process_cached_sample, cache._samplefiles_to_profile
                ):
                    pbar.update(1)

            cache.import_cached_samples()

        # importing properties
        if args.tsv:
            with sonarDBManager(args.db, readonly=False, debug=args.debug) as dbm:
                for sample_name in tqdm(
                    properties,
                    desc="import data ...",
                    total=len(properties),
                    unit="samples",
                    bar_format="{desc} {percentage:3.0f}% [{n_fmt}/{total_fmt}, {elapsed}<{remaining}, {rate_fmt}{postfix}]",
                    disable=args.no_progress,
                ):
                    sample_id = dbm.get_sample_id(sample_name)
                    if not sample_id:
                        continue
                    for property_name, value in properties[sample_name].items():
                        dbm.insert_property(sample_id, property_name, value)

    # view-prop
    elif args.tool == "list-prop":
        with sonarDBManager(args.db, debug=debug) as dbm:
            if not dbm.properties:
                print("*** no properties ***")
                exit(0)

            cols = [
                "name",
                "argument",
                "description",
                "data type",
                "query type",
                "standard value",
            ]
            rows = []
            for prop in sorted(dbm.properties.keys()):
                dt = (
                    dbm.properties[prop]["datatype"]
                    if dbm.properties[prop]["datatype"] != "float"
                    else "decimal number"
                )
                rows.append([])
                rows[-1].append(prop)
                rows[-1].append("--" + prop)
                rows[-1].append(fill(dbm.properties[prop]["description"], width=25))
                rows[-1].append(dt)
                rows[-1].append(dbm.properties[prop]["querytype"])
                rows[-1].append(dbm.properties[prop]["standard"])

            print(tabulate(rows, headers=cols, tablefmt="fancy_grid"))
            print()
            print("DATE FORMAT")
            print("dates must comply with the following format: YYYY-MM-DD")
            print()
            print("OPERATORS")
            print(
                "integer, floating point and decimal number data types support the following operators prefixed directly to the respective value without spaces:"
            )
            print("  > larger than (e.g. >1)")
            print("  < smaller than (e.g. <1)")
            print("  >= larger than or equal to (e.g. >=1)")
            print("  <= smaller than or equal to (e.g. <=1)")
            print("  != different than (e.g. !=1)")
            print()
            print("RANGES")
            print(
                "integer, floating point and date data types support ranges defined by two values directly connected by a colon (:) with no space between them:"
            )
            print("  e.g. 1:10 (between 1 and 10)")
            print("  e.g. 2021-01-01:2021-12-31 (between 1st Jan and 31st Dec of 2021)")
            print()

    # add-prop
    elif args.tool == "add-prop":
        with sonarDBManager(args.db, readonly=False, debug=args.debug) as dbm:
            if args.qtype is None:
                if args.dtype == "integer":
                    args.qtype = "numeric"
                elif args.dtype == "float":
                    args.qtype = "numeric"
                elif args.dtype == "text":
                    args.qtype = "text"
                elif args.dtype == "date":
                    args.qtype = "date"
                elif args.dtype == "zip":
                    args.qtype = "zip"
            dbm.add_property(
                args.name, args.dtype, args.qtype, args.descr, args.default
            )
        print("Inserted successfully:", args.name)

    # delprop
    elif args.tool == "delete-prop":
        with sonarDBManager(args.db, readonly=False, debug=debug) as dbm:
            if args.name not in dbm.properties:
                sys.exit("input error: unknown property.")
            a = dbm.count_property(args.name)
            b = dbm.count_property(args.name, ignore_standard=True)
            if args.force:
                decision = "YES"
            else:
                print(
                    "WARNING: There are",
                    a,
                    "samples with content for this property. Amongst those,",
                    b,
                    "samples do not share the default value of this property.",
                )
                decision = ""
                while decision not in ("YES", "no"):
                    decision = input(
                        "Do you really want to delete this property? [YES/no]: "
                    )
            if decision == "YES":
                dbm.delete_property(args.name)
                print("property deleted.")
            else:
                print("property not deleted.")

    # delete
    elif args.tool == "delete":
        samples = set([x.strip() for x in args.sample])
        for file in args.sample_file:
            check_file(file)
            with sonarBasics.open_file(file, compressed="auto") as handle:
                for line in handle:
                    samples.add(line.strip())
        if len(samples) == 0:
            print("Nothing to delete.")
        else:
            sonarBasics.delete(args.db, *samples, debug=args.debug)

    # restore
    elif args.tool == "restore":
        samples = set([x.strip() for x in args.sample])
        for file in args.sample_file:
            check_file(file)
            with sonarBasics.open_file(file, compressed="auto") as handle:
                for line in handle:
                    samples.add(line.strip())
        if len(samples) == 0:
            print("Nothing to restore.")
        else:
            sonarBasics.restore(
                args.db, *samples, aligned=args.aligned, debug=args.debug
            )

    # update-lineage-info
    elif args.tool == "update-lineage-info":
        print("Start to update parent-child relationship")
        fname = os.path.join(
            os.path.join(os.path.dirname(os.path.realpath(__file__))),
            "data/lineage.all.tsv",
        )
        # obj = sonarLinmgr()
        # lin_df = obj.update_lineage_data(fname)
        with sonarLinmgr() as obj:
            lin_df = obj.update_lineage_data(fname)

        with sonarDBManager(args.db, readonly=False, debug=args.debug) as dbm:
            dbm.add_update_lineage(lin_df)
            print("Update has been successfully")

    # info
    elif args.tool == "info":
        sonarBasics.show_db_info(args.db)

    # match
    elif args.tool == "match":
        props = {}
        reserved_props = {}

        with sonarDBManager(args.db, readonly=False, debug=args.debug) as dbm:
            for pname in dbm.properties:
                if hasattr(args, pname):
                    props[pname] = getattr(args, pname)
            if args.with_sublineage:
                if args.with_sublineage in dbm.properties:
                    reserved_props["with_sublineage"] = args.with_sublineage
                else:
                    sys.exit(
                        "input error: with-sublineage value is mismatch to the available properties"
                    )

        # for reserved keywords
        reserved_key = ["sample"]
        for pname in reserved_key:
            if hasattr(args, pname):
                if pname == "sample" and len(getattr(args, pname)) > 0:
                    # reserved_props[pname] = set([x.strip() for x in args.sample])
                    reserved_props = sonarBasics.set_key(
                        reserved_props, pname, getattr(args, pname)
                    )
                    # reserved_props[pname] = getattr(args, pname)

        # Support file upload
        if args.sample_file:
            for sample_file in args.sample_file:
                check_file(sample_file)
                with sonarBasics.open_file(sample_file, compressed="auto") as file:
                    for line in file:
                        reserved_props = sonarBasics.set_key(
                            reserved_props, "sample", line.strip()
                        )
        format = "count" if args.count else args.format

        sonarBasics.match(
            args.db,
            args.profile,
            reserved_props,
            props,
            outfile=args.out,
            debug=args.debug,
            format=format,
        )

    # optimize
    if args.tool == "optimize":
        with sonarDBManager(args.db, debug=args.debug) as dbm:
            dbm.optimize(args.db)

    # dev
    if args.tool == "dev":
        print("***dev mode***")
        with sonarDBManager(args.db, debug=debug) as dbm:
            for feature in dbm.get_annotation():
                print(())
    # Finished successfully
    return 0


def run():
    parsed_args = parse_args(sys.argv[1:])
    main(parsed_args)


if __name__ == "__main__":
    run()
