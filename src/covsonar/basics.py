#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# author: Stephan Fuchs (Robert Koch Institute, MF1, fuchss@rki.de)

# DEPENDENCIES
from contextlib import contextmanager
import gzip
from hashlib import sha256
import logging
import lzma
import os
import sys
from typing import Union
import zipfile

from Bio.Seq import Seq
import magic

from covsonar import __version__


# CLASS
class sonarBasics:
    """
    A class providing basic operations to covSonar modules.
    """

    # VERSION CONTROL
    @staticmethod
    def get_version() -> str:
        """
        Retrieves the version of the covSonar package.
        """
        return __version__

    # LOGGING
    @staticmethod
    def set_logging_config() -> None:
        logging.basicConfig(
            format="%(asctime)s %(levelname)-4s: %(message)s",
            level=logging.WARNING,
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        logging.getLogger("requests").setLevel(logging.WARNING)

    # FILE HANDLING
    @staticmethod
    @contextmanager
    def open_file_autodetect(file_path: str, mode: str = "r"):
        """
        Opens a file with automatic packing detection.

        Args:
            file_path: The path of the file to open.
            mode: The mode in which to open the file. Default is 'r' (read mode).

        Returns:
            A context manager yielding a file object.
        """
        # Use the magic library to identify the file type
        file_type = magic.from_file(file_path, mime=True)

        if file_type == "application/x-xz":
            file_obj = lzma.open(file_path, mode + "t")  # xz
        elif file_type == "application/gzip":
            file_obj = gzip.open(file_path, mode + "t")  # gz
        elif file_type == "application/zip":
            zip_file = zipfile.ZipFile(file_path, mode)  # zip
            # Assumes there's one file in the ZIP, adjust as necessary
            file_obj = zip_file.open(zip_file.namelist()[0], mode)
        elif file_type == "text/plain":  # plain
            file_obj = open(file_path, mode)
        else:
            raise ValueError(f"Unsupported file type: {file_type}")

        try:
            yield file_obj
        finally:
            file_obj.close()
            if file_type == "application/zip":
                zip_file.close()

    @staticmethod
    def _files_exist(*files: str) -> bool:
        """Check if files exits."""
        for fname in files:
            if not os.path.isfile(fname):
                return False
        return True

    @staticmethod
    @contextmanager
    def out_autodetect(outfile=None):
        """
        Open a file if the 'outfile' is provided.
        If not, use standard output.

        Args:
            outfile: File path to the output file. If None, use standard output.
        """
        if outfile is not None:
            f = open(outfile, "w")
        else:
            f = sys.stdout
        try:
            yield f
        finally:
            if f is not sys.stdout:
                f.close()

    # SEQUENCE HANDLING
    def harmonize_seq(seq: str) -> str:
        """
        Harmonizes the input sequence.

        This function trims leading and trailing white spaces, converts the sequence to upper case and
        replaces all occurrences of "U" with "T". It's usually used to standardize the format of a DNA
        or RNA sequence.

        Args:
            seq (str): The input sequence as a string.

        Returns:
            str: The harmonized sequence.
        """
        try:
            return seq.strip().upper().replace("U", "T")
        except AttributeError as e:
            raise ValueError(
                f"Invalid input, expected a string, got {type(seq).__name__}"
            ) from e

    @staticmethod
    def hash_seq(seq: Union[str, Seq]) -> str:
        """
        Generate a hash from a sequence.

        Args:
            seq: The sequence to hash. This can either be a string or a Seq object from BioPython.

        Returns:
            The SHA-256 hash of the sequence.

        Notes:
            The SHA-256 hash algorithm is used as it provides a good balance
            between performance and collision resistance.
        """
        # If the input is a BioPython Seq object, convert it to a string.
        if isinstance(seq, Seq):
            seq = str(seq)

        return sha256(seq.upper().encode()).hexdigest()
