#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

import argparse
import os
import sys
from textwrap import fill
from typing import Optional

from tabulate import tabulate

from covsonar.basics import sonarBasics
from covsonar.cache import sonarCache  # noqa: F401
from covsonar.dbm import sonarDbManager
from covsonar.linmgr import sonarLinmgr
from covsonar.logging import LoggingConfigurator
from covsonar.utils import sonarUtils


# Constants
VERSION = sonarBasics.get_version()

# Initialize logger
LOGGER = LoggingConfigurator.get_logger()


class args_namespace:
    """An empty class for storing command-line arguments as object attributes."""

    pass


def parse_args(args=None):
    """
    Parse command-line arguments using argparse.ArgumentParser.

    Args:
        args (list): List of command-line arguments. Default is None.
        namespace (argparse.Namespace): An existing namespace to populate with parsed arguments. Default is None.

    Returns:
        argparse.Namespace: Namespace containing parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        prog="sonar",
        description="covSonar2: Integrated Genome Information System for Pathogen Surveillance.",
    )

    # Create parent parsers for common arguments and options for the command-line interface
    database_parser = create_parser_database()
    output_parser = create_parser_output()
    sample_parser = create_parser_sample()
    property_parser = create_parser_property()
    reference_parser = create_parser_reference()
    thread_parser = create_parser_thread()

    # Create all subparsers for the command-line interface
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers, _ = create_subparser_setup(
        subparsers, database_parser, reference_parser
    )
    subparsers, _ = create_subparser_import(subparsers, database_parser, thread_parser)
    subparsers, _ = create_subparser_list_prop(subparsers, database_parser)
    subparsers, _ = create_subparser_add_prop(
        subparsers, database_parser, property_parser
    )
    subparsers, _ = create_subparser_delete_prop(
        subparsers, database_parser, property_parser
    )
    subparsers, _ = create_subparser_delete(subparsers, database_parser, sample_parser)
    subparsers, subparser_match = create_subparser_match(
        subparsers, database_parser, reference_parser, sample_parser, output_parser
    )
    subparsers, _ = create_subparser_restore(
        subparsers, database_parser, sample_parser, output_parser
    )
    subparsers, _ = create_subparser_info(subparsers, database_parser)
    subparsers, _ = create_subparser_optimize(subparsers, database_parser)
    subparsers, _ = create_subparser_db_upgrade(subparsers, database_parser)
    subparsers, _ = create_subparser_update_pangolin(subparsers, database_parser)
    subparsers, _ = create_subparser_direct_query(
        subparsers, database_parser, output_parser
    )

    # add version option
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="covsonar " + VERSION,
        help="show program's version number and exit",
    )

    # add database-specific properties to match subparser
    user_namespace = args_namespace()
    known_args, _ = parser.parse_known_args(args=args, namespace=user_namespace)
    if is_match_selected(known_args):
        LoggingConfigurator(debug=known_args.debug)
        with sonarDbManager(known_args.db, readonly=True) as db_manager:
            for property in db_manager.properties.values():
                subparser_match.add_argument(
                    "--" + property["name"], type=str, nargs="+"
                )

    return parser.parse_args(args=args, namespace=user_namespace)


def create_parser_database() -> argparse.ArgumentParser:
    """Creates a 'database' parent parser with common arguments and options for the command-line interface.

    Returns:
        argparse.ArgumentParser: The created 'database' parent parser.
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "--db", metavar="FILE", help="path to Sonar database", type=str, required=True
    )
    parser.add_argument(
        "--debug",
        help="activate debugging mode and show all SQLite queries on screen",
        action="store_true",
    )
    return parser


def create_parser_output() -> argparse.ArgumentParser:
    """Creates an 'output' parent parser with common arguments and options for the command-line interface.

    Returns:
        argparse.ArgumentParser: The created 'output' parent parser.
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "-o",
        "--out",
        metavar="FILE",
        help="write output file (existing files will be overwritten!)",
        type=str,
        default=None,
    )
    return parser


def create_parser_sample() -> argparse.ArgumentParser:
    """Creates a 'sample' parent parser with common arguments and options for the command-line interface.

    Returns:
        argparse.ArgumentParser: The created 'sample' parent parser.
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "--sample",
        metavar="STR",
        help="sample accession(s) to consider",
        type=str,
        nargs="+",
        default=[],
    )
    parser.add_argument(
        "--sample-file",
        metavar="FILE",
        help="file containing sample accession(s) to consider (one per line)",
        type=str,
        nargs="+",
        default=[],
    )
    return parser


def create_parser_property() -> argparse.ArgumentParser:
    """Creates a 'property' parent parser with common arguments and options for the command-line interface.

    Returns:
        argparse.ArgumentParser: The created 'property' parent parser.
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "--name", metavar="STR", help="property name", type=str, required=True
    )
    return parser


def create_parser_reference() -> argparse.ArgumentParser:
    """Creates a 'reference' parent parser with common arguments and options for the command-line interface.

    Returns:
        argparse.ArgumentParser: The created 'reference' parent parser.
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "-r",
        "--reference",
        metavar="STR",
        help="reference accession",
        type=str,
        default=None,
    )
    return parser


def create_parser_thread() -> argparse.ArgumentParser:
    """Creates a 'thread' parent parser with common arguments and options for the command-line interface.

    Returns:
        argparse.ArgumentParser: The created 'thread' parent parser.
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "-t",
        "--threads",
        metavar="INT",
        help="number of threads to use (default: 1)",
        type=int,
        default=1,
    )
    return parser


def create_subparser_setup(
    subparsers: argparse._SubParsersAction, *parent_parsers: argparse.ArgumentParser
) -> argparse.ArgumentParser:
    """
    Creates a 'setup' subparser with command-specific arguments and options for the command-line interface.

    Args:
        subparsers (argparse._SubParsersAction): ArgumentParser object to attach the 'setup' subparser to.
        parent_parsers (argparse.ArgumentParser): ArgumentParser objects providing common arguments and options.

    Returns:
        argparse.ArgumentParser: The created 'setup' subparser.
    """
    parser = subparsers.add_parser(
        "setup", help="set up a new database", parents=parent_parsers
    )
    parser.add_argument(
        "--default-props",
        help="add commonly used properties to the new database",
        action="store_true",
    )
    parser.add_argument(
        "--gbk",
        metavar="GBK_FILE",
        help="path to GenBank reference file ( MN908947.3 is used as default reference)",
        type=str,
        default=None,
    )
    return subparsers, parser


def create_subparser_import(
    subparsers: argparse._SubParsersAction, *parent_parsers: argparse.ArgumentParser
) -> argparse.ArgumentParser:
    """
    Creates an 'import' subparser with command-specific arguments and options for the command-line interface.

    Args:
        subparsers (argparse._SubParsersAction): ArgumentParser object to attach the 'import' subparser to.
        parent_parsers (argparse.ArgumentParser): ArgumentParser objects providing common arguments and options.

    Returns:
        argparse.ArgumentParser: The created 'import' subparser.
    """
    parser = subparsers.add_parser(
        "import",
        help="import genome sequences and sample properties into the database",
        parents=parent_parsers,
    )
    parser.add_argument(
        "--fasta",
        metavar="FASTA_FILE",
        help="file containing genome sequences to import",
        type=str,
        nargs="+",
        default=[],
    )
    parser.add_argument(
        "--tsv",
        metavar="TSV_FILE",
        help="tab-delimited file containing sample properties to import",
        type=str,
        nargs="+",
        default=[],
    )
    parser.add_argument(
        "--csv",
        metavar="CSV_FILE",
        help="comma-delimited file containing sample properties to import",
        type=str,
        nargs="+",
        default=[],
    )
    parser.add_argument(
        "--cols",
        metavar="PROP=COL",
        help="assign column names used in the provided TSV/CSV file to the matching property names provided by the database in the form PROP=COL (e.g. SAMPLE=GenomeID)",
        type=str,
        nargs="+",
        default=[],
    )
    parser.add_argument(
        "--auto-link",
        help="automatically link TSV/CSV columns with database fields based on identical names",
        action="store_true",
    )
    parser.add_argument(
        "--no-update",
        help="do not update sequences or properties of samples already existing in the database",
        action="store_true",
    )
    parser.add_argument(
        "--cache",
        metavar="DIR",
        help="directory for caching data (default: a temporary directory is created)",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--no-progress",
        help="do not show progress bars while importing",
        action="store_true",
    )
    return subparsers, parser


def create_subparser_list_prop(
    subparsers: argparse._SubParsersAction, *parent_parsers: argparse.ArgumentParser
) -> argparse.ArgumentParser:
    """
    Creates a 'list-prop' subparser with command-specific arguments and options for the command-line interface.

    Args:
        subparsers (argparse._SubParsersAction): ArgumentParser object to attach the 'list-prop' subparser to.
        parent_parsers (argparse.ArgumentParser): ArgumentParser objects providing common arguments and options.

    Returns:
        argparse.ArgumentParser: The created 'list-prop' subparser.
    """
    parser = subparsers.add_parser(
        "list-prop",
        help="view sample properties added to the database",
        parents=parent_parsers,
    )
    return subparsers, parser


def create_subparser_add_prop(
    subparsers: argparse._SubParsersAction, *parent_parsers: argparse.ArgumentParser
) -> argparse.ArgumentParser:
    """
    Creates an 'add-prop' subparser with command-specific arguments and options for the command-line interface.

    Args:
        subparsers (argparse._SubParsersAction): ArgumentParser object to attach the 'add-prop' subparser to.
        parent_parsers (argparse.ArgumentParser): ArgumentParser objects providing common arguments and options.

    Returns:
        argparse.ArgumentParser: The created 'add-prop' subparser.
    """
    parser = subparsers.add_parser(
        "add-prop", help="add a property to the database", parents=parent_parsers
    )
    parser.add_argument(
        "--descr",
        metavar="STR",
        help="a short description of the property",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--dtype",
        metavar="STR",
        help="the data type of the property",
        type=str,
        choices=["integer", "float", "text", "date", "zip", "pango"],
        required=True,
    )
    parser.add_argument(
        "--qtype",
        metavar="STR",
        help="the query type of the property",
        type=str,
        choices=["numeric", "float", "text", "date", "zip", "pango"],
        default=None,
    )
    parser.add_argument(
        "--default",
        metavar="VAR",
        help="the default value of the property (none by default)",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--subject",
        metavar="VAR",
        help="choose between sample or variant property (by default: sample)",
        choices=["sample", "variant"],
        default="sample",
    )
    return subparsers, parser


def create_subparser_delete_prop(
    subparsers: argparse._SubParsersAction, *parent_parsers: argparse.ArgumentParser
) -> argparse.ArgumentParser:
    """
    Creates an 'delete-prop' subparser with command-specific arguments and options for the command-line interface.

    Args:
        subparsers (argparse._SubParsersAction): ArgumentParser object to attach the 'delete-prop' subparser to.
        parent_parsers (argparse.ArgumentParser): ArgumentParser objects providing common arguments and options.

    Returns:
        argparse.ArgumentParser: The created 'delete-prop' subparser.
    """
    parser = subparsers.add_parser(
        "delete-prop", help="delete a property to the database", parents=parent_parsers
    )
    parser.add_argument(
        "--force",
        help="skip user confirmation and force property to be deleted",
        action="store_true",
    )
    return subparsers, parser


def create_subparser_delete(
    subparsers: argparse._SubParsersAction, *parent_parsers: argparse.ArgumentParser
) -> argparse.ArgumentParser:
    """
    Creates an 'delete' subparser with command-specific arguments and options for the command-line interface.

    Args:
        subparsers (argparse._SubParsersAction): ArgumentParser object to attach the 'delete' subparser to.
        parent_parsers (argparse.ArgumentParser): ArgumentParser objects providing common arguments and options.

    Returns:
        argparse.ArgumentParser: The created 'delete' subparser.
    """
    parser = subparsers.add_parser(
        "delete", help="delete samples from the database", parents=parent_parsers
    )
    return subparsers, parser


def create_subparser_match(
    subparsers: argparse._SubParsersAction, *parent_parsers: argparse.ArgumentParser
) -> argparse.ArgumentParser:
    """
    Creates a 'match' subparser with command-specific arguments and options for the command-line interface.

    Args:
        subparsers (argparse._SubParsersAction): An ArgumentParser object to attach the 'match' subparser to.
        parent_parsers (argparse.ArgumentParser): A list of ArgumentParser objects providing common arguments and options.

    Returns:
        argparse.ArgumentParser: The created 'match' subparser.
    """
    parser = subparsers.add_parser(
        "match",
        help="match samples based on mutation profiles and/or properties",
        parents=parent_parsers,
    )
    parser.add_argument(
        "--profile",
        "-p",
        metavar="STR",
        help="match genomes sharing the given mutation profile",
        type=str,
        action="append",
        nargs="+",
        default=[],
    )
    parser.add_argument(
        "--showNX",
        help="include non-informative polymorphisms in resulting mutation profiles (X for AA and N for NT)",
        action="store_true",
    )
    parser.add_argument(
        "--frameshifts-only",
        help="match only mutation profiles with frameshift mutations",
        action="store_true",
    )
    parser.add_argument(
        "--out-cols",
        metavar="STR",
        help="define output columns for CSV and TSV files (by default all available columns are shown)",
        type=str,
        nargs="+",
        default=[],
    )
    mutually_exclusive_group = parser.add_mutually_exclusive_group()
    mutually_exclusive_group.add_argument(
        "--count", help="count matching genomes only", action="store_true"
    )
    mutually_exclusive_group.add_argument(
        "--format",
        help="output format (default: tsv)",
        choices=["csv", "tsv", "vcf"],
        default="tsv",
    )
    return subparsers, parser


def create_subparser_restore(
    subparsers: argparse._SubParsersAction, *parent_parsers: argparse.ArgumentParser
) -> argparse.ArgumentParser:
    """
    Creates a 'restore' subparser with command-specific arguments and options for the command-line interface.

    Args:
        subparsers (argparse._SubParsersAction): An ArgumentParser object to attach the 'restore' subparser to.
        parent_parsers (argparse.ArgumentParser): A list of ArgumentParser objects providing common arguments and options.

    Returns:
        argparse.ArgumentParser: The created 'restore' subparser.
    """
    parser = subparsers.add_parser(
        "restore",
        help="restore genome sequences from the database",
        parents=parent_parsers,
    )
    parser.add_argument(
        "--aligned",
        help='use pairwise aligned form (deletions indicated by "-" and insertions by lowercase letters)',
        action="store_true",
    )
    return subparsers, parser


def create_subparser_info(
    subparsers: argparse._SubParsersAction, *parent_parsers: argparse.ArgumentParser
) -> argparse.ArgumentParser:
    """
    Creates a 'info' subparser with command-specific arguments and options for the command-line interface.

    Args:
        subparsers (argparse._SubParsersAction): An ArgumentParser object to attach the 'info' subparser to.
        parent_parsers (argparse.ArgumentParser): A list of ArgumentParser objects providing common arguments and options.

    Returns:
        argparse.ArgumentParser: The created 'info' subparser.
    """
    parser = subparsers.add_parser(
        "info",
        help="show detailed information on a given database",
        parents=parent_parsers,
    )
    parser.add_argument(
        "--detailed",
        help="show numbers of stored mutations (dependent on database size this might take a while to process)",
        action="store_true",
    )
    return subparsers, parser


def create_subparser_optimize(
    subparsers: argparse._SubParsersAction, *parent_parsers: argparse.ArgumentParser
) -> argparse.ArgumentParser:
    """
    Creates an 'optimize' subparser with command-specific arguments and options for the command-line interface.

    Args:
        subparsers (argparse._SubParsersAction): ArgumentParser object to attach the 'optimize' subparser to.
        parent_parsers (argparse.ArgumentParser): ArgumentParser objects providing common arguments and options.

    Returns:
        argparse.ArgumentParser: The created 'optimize' subparser.
    """
    parser = subparsers.add_parser(
        "optimize", help="optimize database", parents=parent_parsers
    )
    parser.add_argument(
        "--tempdir",
        help="custom temporrary direcxtory (default: None)",
        type=str,
        default=None,
    )
    return subparsers, parser


def create_subparser_db_upgrade(
    subparsers: argparse._SubParsersAction, *parent_parsers: argparse.ArgumentParser
) -> argparse.ArgumentParser:
    """
    Creates an 'db-upgrade' subparser with command-specific arguments and options for the command-line interface.

    Args:
        subparsers (argparse._SubParsersAction): ArgumentParser object to attach the 'db-upgrade' subparser to.
        parent_parsers (argparse.ArgumentParser): ArgumentParser objects providing common arguments and options.

    Returns:
        argparse.ArgumentParser: The created 'db-upgrade' subparser.
    """
    parser = subparsers.add_parser(
        "db-upgrade",
        help="upgrade the database to the latest version",
        parents=parent_parsers,
    )
    return subparsers, parser


def create_subparser_update_pangolin(
    subparsers: argparse._SubParsersAction, *parent_parsers: argparse.ArgumentParser
) -> argparse.ArgumentParser:
    """
    Creates an 'update-lineages' subparser with command-specific arguments and options for the command-line interface.

    Args:
        subparsers (argparse._SubParsersAction): ArgumentParser object to attach the 'update-lineages' subparser to.
        parent_parsers (argparse.ArgumentParser): ArgumentParser objects providing common arguments and options.

    Returns:
        argparse.ArgumentParser: The created 'update-lineages' subparser.
    """
    parser = subparsers.add_parser(
        "update-lineages",
        help="download latest pangolin information",
        parents=parent_parsers,
    )
    parser.add_argument(
        "--alias-key",
        help="Pangolin alias_key.json file (default: auto download from GitHub)",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--lineages",
        help="Pangolin lineages.csv file (default: auto download from GitHub)",
        type=str,
        default=None,
    )
    return subparsers, parser


def create_subparser_direct_query(
    subparsers: argparse._SubParsersAction, *parent_parsers: argparse.ArgumentParser
) -> argparse.ArgumentParser:
    """
    Creates an 'direct-query' subparser with command-specific arguments and options for the command-line interface.

    Args:
        subparsers (argparse._SubParsersAction): ArgumentParser object to attach the 'direct-query' subparser to.
        parent_parsers (argparse.ArgumentParser): ArgumentParser objects providing common arguments and options.

    Returns:
        argparse.ArgumentParser: The created 'direct-query' subparser.
    """
    parser = subparsers.add_parser(
        "direct-query",
        help="connect read-only to the database for direct queries",
        parents=parent_parsers,
    )
    parser.add_argument(
        "--sql",
        help="sqlite query",
        type=str,
        required=True,
    )
    return subparsers, parser


def is_match_selected(namespace: Optional[argparse.Namespace] = None) -> bool:
    """
    Checks if the 'match' command is selected and the 'db' attribute is present in the arguments.

    Args:
        namespace: Namespace object for storing argument values (default: None)

    Returns:
        True if 'match' command is selected and 'db' attribute is present, False otherwise
    """
    # Check if the 'match' command is selected and the 'db' attribute is present
    match_selected = namespace.command == "match" and hasattr(namespace, "db")

    return match_selected


def check_file(fname, exit_on_fail=True):
    """
    Check if a given file path exists.

    Args:
        fname (string): The name and path to an existing file.
        exit_on_fail (boolean): Whether to exit the script if the file doesn't exist. Default is True.

    Returns:
        True if the file exists, False otherwise.
    """
    if not os.path.isfile(fname):
        if exit_on_fail:
            sys.exit("Error: The file '" + fname + "' does not exist.")
        return False
    return True


def handle_setup(args: argparse.Namespace):
    """
    Handle database setup.

    Args:
        args (argparse.Namespace): Parsed command line arguments.

    Raises:
        FileNotFoundError: If the specified GenBank file is not found.
    """
    if args.gbk:
        check_file(args.gbk)
    sonarUtils.setup_db(args.db, args.default_props, reference_gb=args.gbk)


def handle_db_upgrade(args: argparse.Namespace):
    """
    Handle database upgrade, prompting the user to confirm the action.

    Args:
        args (argparse.Namespace): Parsed command line arguments.
    """
    LOGGER.warning("Backup your database before upgrading")
    decision = ""
    while decision not in ("YES", "no"):
        decision = input("Do you really want to perform this action? [YES/no]: ")
    if decision == "YES":
        sonarDbManager.upgrade_db(args.db)
    else:
        LOGGER.info("No operation is performed")


def handle_import(args: argparse.Namespace):
    """
    Handle data import.

    Args:
        args (argparse.Namespace): Parsed command line arguments.
    """
    sonarUtils.import_data(
        db=args.db,
        fasta=args.fasta,
        csv_files=args.csv,
        tsv_files=args.tsv,
        prop_links=args.cols,
        cachedir=args.cache,
        autolink=args.auto_link,
        progress=not args.no_progress,
        update=not args.no_update,
        threads=args.threads,
    )


def handle_list_prop(args: argparse.Namespace):
    """
    Handle listing properties from the database.
    This function retrieves all properties stored in the database, sorts them by
    name, and formats the output as a table. The table includes columns such as
    name, argument, subject, description, data type, query type, and standard value.
    The formatted table is then printed to the console.

    Args:
        args (argparse.Namespace): Parsed command line arguments.
    """
    with sonarDbManager(args.db) as db_manager:
        if not db_manager.properties:
            print("*** no properties ***")
        else:
            cols = [
                "name",
                "argument",
                "subject",
                "description",
                "data type",
                "query type",
                "standard value",
            ]
            rows = []
            for prop in sorted(db_manager.properties.keys()):
                dt = (
                    db_manager.properties[prop]["datatype"]
                    if db_manager.properties[prop]["datatype"] != "float"
                    else "decimal number"
                )
                rows.append([])
                rows[-1].append(prop)
                rows[-1].append("--" + prop)
                rows[-1].append("sample")
                rows[-1].append(
                    fill(db_manager.properties[prop]["description"], width=25)
                )
                rows[-1].append(dt)
                rows[-1].append(db_manager.properties[prop]["querytype"])
                rows[-1].append(db_manager.properties[prop]["standard"])

            print(tabulate(rows, headers=cols, tablefmt="fancy_grid"))


def handle_add_prop(args: argparse.Namespace):
    """
    Handle adding a new property to the database.

    Args:
        args (argparse.Namespace): Parsed command line arguments.
    """
    with sonarDbManager(args.db, readonly=False) as db_manager:
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
            elif args.dtype == "pango":
                args.dtype = "text"
                args.qtype = "pango"
        db_manager.add_property(
            args.name,
            args.dtype,
            args.qtype,
            args.descr,
            args.subject,
            args.default,
        )
    LOGGER.info("Inserted successfully: %s", args.name)


def handle_delete_prop(args: argparse.Namespace):
    """
    Handle deleting an existing property from the database.

    This function removes a specified property from the database. If the 'force'
    option is not set, the user is prompted to confirm the deletion, especially
    when there are samples with non-default values for the property.

    Args:
    args (argparse.Namespace): Parsed command line arguments containing the
    property name to be deleted and the 'force' flag.

    Raises:
    SystemExit: If the specified property name is not found in the database.
    """
    with sonarDbManager(args.db, readonly=False) as db_manager:
        if args.name not in db_manager.properties:
            LOGGER.error("Unknown property.")
            sys.exit(1)

        non_default_values_count = db_manager.count_property(
            args.name, ignore_standard=True
        )

        if args.force:
            decision = "YES"
        else:
            LOGGER.warning(
                f"There are {non_default_values_count} samples with non-default values for this property."
            )
            decision = ""
            while decision not in ("YES", "no"):
                decision = input(
                    "Do you really want to delete this property? [YES/no]: "
                )

        if decision.lower() == "yes":
            db_manager.delete_property(args.name)
            LOGGER.info("Property deleted.")
        else:
            LOGGER.info("Property not deleted.")


def handle_delete(args: argparse.Namespace):
    """
    Handle deleting a sample from the database.

    This function removes specified samples from the database. The samples can be
    provided as command line arguments or within a file. If the 'force' option is
    not set, the user is prompted to confirm the deletion, especially when there
    are samples with non-default values for any property.

    Args:
        args (argparse.Namespace): Parsed command line arguments containing the
                                        sample names to be deleted and the 'force' flag.

    Raises:
        FileNotFoundError: If the specified sample file is not found.
        SystemExit: If none of the specified samples are found in the database.
    """
    samples = set([x.strip() for x in args.sample])
    for fname in args.sample_file:
        check_file(fname)
        with sonarBasics.open_file_autodetect(fname) as handle:
            for line in handle:
                samples.add(line.strip())
    if len(samples) == 0:
        LOGGER.info("Nothing to delete.")
    else:
        sonarUtils.delete_sample(args.db, samples)


def handle_restore(args: argparse.Namespace):
    """
    Handle restoring genome sequences from the database.

    This function retrieves specified genome sequences from the database and
    writes them to an output file. The samples can be provided as command line
    arguments or within a file. Optionally, the output can contain aligned
    sequences.

    Args:
        args (argparse.Namespace): Parsed command line arguments containing the
                                        sample names to be restored and optional flags
                                        such as 'aligned' and 'out'.
    Raises:
        FileNotFoundError: If the specified sample file is not found.
    """
    samples = set([x.strip() for x in args.sample])
    for file in args.sample_file:
        check_file(file)
        with sonarBasics.open_file_autodetect(file) as handle:
            for line in handle:
                samples.add(line.strip())
    sonarUtils.restore_seq(
        args.db,
        samples,
        aligned=args.aligned,
        outfile=args.out,
    )


def handle_update_pangolin(args: argparse.Namespace):
    """
    Handle update pangolin information.

    This function updates the parent-child relationship of pangolin lineages
    stored in the database by obtaining the latest lineage data.

    Args:
        args (argparse.Namespace): Parsed command line arguments.

    Raises:
        FileNotFoundError: If any required files are not found.
    """
    with sonarLinmgr() as lineage_manager:
        lineage_data = lineage_manager.update_lineage_data(
            args.alias_key, args.lineages
        )

    with sonarDbManager(args.db, readonly=False) as db_manager:
        db_manager.add_update_lineage(lineage_data)


def handle_info(args: argparse.Namespace):
    """
    Handle displaying database information.

    This function retrieves and displays the database information, including
    various statistics and metadata about the stored data.

    Args:
        args (argparse.Namespace): Parsed command line arguments.

    Raises:
        FileNotFoundError: If the database file is not found.
    """
    sonarUtils.show_db_info(args.db, args.detailed)


def handle_match(args: argparse.Namespace):
    """
    Handle profile and property matching.

    This function matches samples based on their properties and a given profile.
    The samples can be provided as command line arguments or within a file.
    The output can be customized by selecting specific output columns.

    Args:
        args (argparse.Namespace): Parsed command line arguments containing the
                                        profile, sample names, output format, and other
                                        optional flags.

    Raises:
        FileNotFoundError: If any required files are not found.
        SystemExit: If any unknown output columns are selected.
    """
    properties = {}
    with sonarDbManager(args.db, readonly=False) as db_manager:
        for property_name in db_manager.properties:
            if hasattr(args, property_name):
                properties[property_name] = getattr(args, property_name)

        # Check output column names and property names
        if len(args.out_cols) > 0:
            if all(item in db_manager.properties for item in args.out_cols):
                # sample.name is fixed for output
                args.out_cols = args.out_cols + ["sample.name"]
            else:
                LOGGER.error(
                    "Unknown output column(s) selected. Please select from: "
                    + ", ".join(db_manager.properties)
                )
                sys.exit(1)

    # Sample name handling
    samples = set(args.sample)
    if args.sample_file:
        for sample_file in args.sample_file:
            check_file(sample_file)
            with sonarBasics.open_file_autodetect(sample_file) as file:
                for line in file:
                    samples.add(line.strip())

    # Set output format
    output_format = "count" if args.count else args.format

    # Perform matching
    sonarUtils.match(
        db=args.db,
        profiles=args.profile,
        reference=args.reference,
        samples=samples,
        properties=properties,
        outfile=args.out,
        output_column=args.out_cols,
        format=output_format,
        showNX=args.showNX,
        frameshifts_only=args.frameshifts_only,
    )


def handle_optimize(args: argparse.Namespace):
    """
    Handle database optimization.

    This function optimizes the database by performing operations that
    improve the performance and reduce the storage space.

    Args:
        arguments (argparse.Namespace): Parsed command line arguments.

    Raises:
        FileNotFoundError: If the database file is not found.
    """
    with sonarDbManager(args.db) as db_manager:
        db_manager.optimize(args.db, tmpdir=args.tempdir)


def handle_direct_query(args: argparse.Namespace):
    """
    Handle read-only SQLite queries.

    This function allows the user to execute read-only SQLite queries
    directly on the database. The results can be written to an output file.

    Args:
        arguments (argparse.Namespace): Parsed command line arguments containing
                                        the SQL query and optional output file.

    Raises:
        FileNotFoundError: If the database file is not found.
    """
    if len(args.sql) > 1 and (
        (args.sql.startswith('"') and args.sql.endswith('"'))
        or (args.sql.startswith("'") and args.sql.endswith("'"))
    ):
        sql = args.sql[1:-1]
    else:
        sql = args.sql
    sonarUtils.direct_query(args.db, query=sql, outfile=args.out)


def execute_commands(args):  # noqa: C901
    """
    Execute the appropriate function based on the provided command.
    This function determines which command was provided as an argument and
    calls the corresponding function to handle the command. It also ensures
    the database compatibility before executing the command.

    Args:
        args (argparse.Namespace): Parsed command line arguments.
    """
    if args.command == "setup":
        handle_setup(args)
    elif args.command == "db-upgrade":
        handle_db_upgrade(args)
    else:
        with sonarDbManager(args.db, readonly=True) as db_manager:
            db_manager.check_db_compatibility()

    if args.command == "import":
        handle_import(args)
    elif args.command == "list-prop":
        handle_list_prop(args)
    elif args.command == "add-prop":
        handle_add_prop(args)
    elif args.command == "delete-prop":
        handle_delete_prop(args)
    elif args.command == "delete":
        handle_delete(args)
    elif args.command == "restore":
        handle_restore(args)
    elif args.command == "update-lineages":
        handle_update_pangolin(args)
    elif args.command == "match":
        handle_match(args)
    elif args.command == "info":
        handle_info(args)
    elif args.command == "optimize":
        handle_optimize(args)
    if args.command == "direct-query":
        handle_direct_query(args)


def main(args: Optional[argparse.Namespace] = None) -> int:
    """
    The main function that handles the execution of different commands.

    Args:
        args (Optional[argparse.Namespace]): Namespace containing parsed command-line arguments.
            If None, the function will parse the arguments itself.

    Returns:
        int: Returns 0 if finished successfully.
    """
    # process arguments
    if not args:
        args = parse_args(sys.argv[1:])

    # Set debugging mode
    if hasattr(args, "debug") and args.debug:
        debug = True
    else:
        debug = False

    LoggingConfigurator(debug=debug)

    # Check database
    if hasattr(args, "db") and args.db:
        if (
            args.command != "setup"
            and args.db is not None
            and not os.path.isfile(args.db)
        ):
            LOGGER.error("The database does not exists")
            sys.exit(1)

    # other commands
    execute_commands(args)

    # Finished successfully
    return 0


if __name__ == "__main__":
    main()
