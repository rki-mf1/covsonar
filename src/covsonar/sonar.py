#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

import argparse
import os
import sys
from textwrap import fill

from tabulate import tabulate

from src.covsonar import logging
from src.covsonar.basics import sonarBasics
from src.covsonar.cache import sonarCache  # noqa: F401
from src.covsonar.dbm import sonarDBManager
from src.covsonar.linmgr import sonarLinmgr


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
        "--out-column",
        help="select output columns to the output file (support csv and tsv)",
        type=str,
        default="all",
    )
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
    parser_match.add_argument(
        "--showNX",
        help="include any variant character in the result (X for AA and N for NT)",
        action="store_true",
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
                # if prop["datatype"] == "integer":
                #    t = int
                # elif prop["datatype"] == "float":
                #    t = float
                # else:
                #    t = str
                if pargs[0].tool == "match":
                    parser_match.add_argument(
                        "--" + prop["name"],
                        # type=t,
                        nargs="+",
                        # default=argparse.SUPPRESS,
                    )

    # return
    return parser.parse_args(args=args, namespace=user_namespace)


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
        print("WARNING: Backup db file before upgrading")
        decision = ""
        while decision not in ("YES", "no"):
            decision = input("Do you really want to perform this action? [YES/no]: ")
        if decision == "YES":
            sonarDBManager.upgrade_db(args.db)
        else:
            logging.info("No operation is performed")
    else:
        with sonarDBManager(args.db, readonly=True) as dbm:
            dbm.check_db_compatibility()
    # other than the above
    # import
    if args.tool == "import":
        sonarBasics.import_data(
            db=args.db,
            fasta=args.fasta,
            tsv=args.tsv,
            cols=args.cols,
            cachedir=args.cache,
            autodetect=not args.no_autodetect,
            progress=not args.no_progress,
            update=not args.no_update,
            threads=args.threads,
            debug=args.debug,
        )

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
        logging.info("Inserted successfully: %s", args.name)

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
                logging.warning(
                    "There are"
                    " %d"
                    " samples with content for this property. Amongst those,"
                    " %d"
                    " samples do not share the default value of this property." % (a, b)
                )
                decision = ""
                while decision not in ("YES", "no"):
                    decision = input(
                        "Do you really want to delete this property? [YES/no]: "
                    )
            if decision == "YES":
                dbm.delete_property(args.name)
                logging.info("property deleted.")
            else:
                logging.info("property not deleted.")

    # delete
    elif args.tool == "delete":
        samples = set([x.strip() for x in args.sample])
        for file in args.sample_file:
            check_file(file)
            with sonarBasics.open_file(file, compressed="auto") as handle:
                for line in handle:
                    samples.add(line.strip())
        if len(samples) == 0:
            logging.info("Nothing to delete.")
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
            logging.info("Nothing to restore.")
        else:
            sonarBasics.restore(
                args.db, *samples, aligned=args.aligned, debug=args.debug
            )

    # update-lineage-info
    elif args.tool == "update-lineage-info":
        logging.info("Start to update parent-child relationship...")
        # fname = os.path.join(
        #    os.path.join(os.path.dirname(os.path.realpath(__file__))),
        #    "data/lineage.all.tsv",
        # )
        # obj = sonarLinmgr()
        # lin_df = obj.update_lineage_data(fname)
        with sonarLinmgr() as obj:
            lin_df = obj.update_lineage_data()

        with sonarDBManager(args.db, readonly=False, debug=args.debug) as dbm:
            dbm.add_update_lineage(lin_df)
            logging.info("Update has been successfully")

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
                        "input error: "
                        + args.with_sublineage
                        + " is mismatch to the available properties"
                    )
            # check column output and property name
            if args.out_column != "all":
                out_column = args.out_column.strip()
                out_column_list = out_column.split(",")
                check = all(item in dbm.properties for item in out_column_list)
                if check:
                    # sample.name is fixed
                    valid_output_column = out_column_list + ["sample.name"]
                else:
                    sys.exit(
                        "input error: "
                        + str(out_column_list)
                        + " one or more given name mismatch the available properties"
                    )
            else:
                valid_output_column = "all"

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
            output_column=valid_output_column,
            debug=args.debug,
            format=format,
            showNX=args.showNX,
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
                print(feature)
    # Finished successfully
    return 0


def run():
    parsed_args = parse_args(sys.argv[1:])
    main(parsed_args)


if __name__ == "__main__":
    run()
