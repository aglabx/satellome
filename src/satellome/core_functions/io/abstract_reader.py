#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 04.06.2011
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com
"""
Abstract base classes and utilities for file and folder I/O operations.

Provides flexible file handling with automatic compression detection (.gz, .bz2),
abstract base classes for building custom readers/writers, and convenient
folder iteration utilities. Designed for genomic data processing workflows.

Classes:
    WiseOpener: Context manager for opening regular, gzip, and bz2 files
    AbstractFileIO: Base class for file I/O with batch and streaming operations
    AbstractFolderIO: Base class for iterating files in folders with regex filtering

Folder Shortcuts:
    sc_iter_filepath_folder: Iterate file paths in folder
    sc_iter_filename_folder: Iterate filenames (no paths)
    sc_iter_folders: Iterate subfolder names
    sc_iter_path_name_folder: Iterate (filename, path) tuples
    sc_iter_filedata_folder: Iterate file contents
    sc_move_files: Move files matching pattern

Processing Shortcuts:
    sc_process_file: Process single file in-place
    sc_process_folder: Process all files in folder in-place
    sc_process_folder_to_other: Process folder to output folder
    read_pickle_file: Load pickled data

Key Features:
    - Automatic compression handling (.gz, .bz2) via WiseOpener
    - Streaming and batch file processing
    - Regex-based file filtering in folders
    - MongoDB integration (deprecated, for legacy support)
    - Functional processing with arbitrary functions
    - Safe atomic file operations

Example:
    >>> # Automatic compression detection
    >>> with WiseOpener("data.txt.gz", "r") as fh:
    ...     lines = fh.readlines()
    >>>
    >>> # Iterate over FASTA files in a folder
    >>> for filepath in sc_iter_filepath_folder("genomes/", mask=r"\\.fa$"):
    ...     process_fasta(filepath)
    >>>
    >>> # Build custom reader
    >>> class MyReader(AbstractFileIO):
    ...     def read_from_file(self, input_file):
    ...         # Custom parsing logic
    ...         pass

See Also:
    satellome.core_functions.io.tab_file: Tab-delimited file I/O
    satellome.core_functions.io.trf_file: TRF format I/O
"""
import bz2
import gzip
import logging
import os
import pickle
import re
import shutil

from satellome.core_functions.exceptions import FileFormatError, ConfigurationError

logger = logging.getLogger(__name__)


class WiseOpener(object):
    """
    Context manager for opening regular, gzip, and bz2 files transparently.

    Automatically detects and handles compressed files (.gz, .bz2) based on
    file extension, while maintaining a consistent API for all file types.
    Ensures proper binary mode for compressed files.

    Attributes:
        file_name (str): Path to file to open
        mode (str): File opening mode ('r', 'w', 'a', 'rb', 'wb', 'ab')
        fh (file object): File handle (set during context entry)

    Example:
        >>> # Transparent gzip handling
        >>> with WiseOpener("data.txt.gz", "r") as fh:
        ...     content = fh.read()
        >>>
        >>> # Regular file (same API)
        >>> with WiseOpener("data.txt", "w") as fh:
        ...     fh.write("Hello")
        >>>
        >>> # Bzip2 files
        >>> with WiseOpener("archive.bz2", "rb") as fh:
        ...     binary_data = fh.read()

    Note:
        - Automatically adds 'b' mode for .gz and .bz2 files
        - Valid modes: 'r', 'w', 'a', 'rb', 'wb', 'ab'
        - Logs compression detection at INFO level
        - Use as context manager (with statement) for proper cleanup
    """

    def __init__(self, file_name, mode=None):
        """
        Initialize file opener with path and mode.

        Args:
            file_name (str): Path to file (can be .gz, .bz2, or regular)
            mode (str, optional): File opening mode. Defaults to "r".
                                 Valid: 'r', 'w', 'a', 'rb', 'wb', 'ab'

        Raises:
            FileFormatError: If mode is invalid

        Note:
            Binary mode ('b' suffix) is automatically added for compressed files
        """
        self.file_name = file_name
        if not mode:
            mode = "r"
        if not mode in ["w", "r", "a", "wb", "rb", "ab"]:
            logger.error("Wrong file mode: %s" % mode)
            raise FileFormatError(
                f"Invalid file mode '{mode}': expected one of ['r', 'w', 'a', 'rb', 'wb', 'ab']. "
                f"Use 'r' for reading, 'w' for writing, 'a' for appending. "
                f"Add 'b' suffix for binary mode (required for .gz and .bz2 files). "
                f"Check file opening mode specification in your code."
            )
        self.mode = mode
        self.fh = None

    def __enter__(self):
        """
        Open file with appropriate handler based on extension.

        Detects compression from extension and opens with:
        - gzip.open() for .gz files
        - bz2.BZ2File() for .bz2 files
        - built-in open() for other files

        Returns:
            file object: Opened file handle

        Note:
            Automatically adds 'b' mode for compressed files if not present
        """
        if self.file_name.endswith(".gz"):
            logger.info("Open as gz archive")
            if not "b" in self.mode:
                self.mode += "b"
            self.fh = gzip.open(self.file_name, self.mode)
        elif self.file_name.endswith(".bz2"):
            logger.info("Open as bz2 archive")
            if not "b" in self.mode:
                self.mode += "b"
            self.fh = bz2.BZ2File(self.file_name, self.mode)
        else:
            self.fh = open(self.file_name, self.mode)
        return self.fh

    def __exit__(self, *args):
        """
        Close file handle when exiting context.

        Args:
            *args: Exception info (exc_type, exc_value, traceback)
                  Passed by context manager protocol

        Note:
            Ensures file is properly closed even if exception occurs
        """
        self.fh.close()


class AbstractFileIO(object):
    """Abstract class for working with abstract data.

    Public properties:

    - data, iterable data
    - N, a number of items in data

    Public methods:

    - read_from_file(self, input_file)
    - read_online(self, input_file) ~> item
    - read_from_db(self, db_cursor) [ABSTRACT]
    - write_to_file(self, output_file)
    - write_to_db(self, db_cursor) [ABSTRACT]
    - read_as_iter(self, source)
    - iterate(self) ~> item of data
    - do(self, cf, **args) -> result
    - process(self, cf, **args)
    - clear(self)
    - do_with_iter(self, cf, **args) -> [result,]
    - process_with_iter(self, cf, **args)
    """

    def __init__(self):
        """Do nothing."""
        self._data = None

    def get_opener(self):
        return WiseOpener

    def read_from_file(self, input_file):
        """Read data from given input_file."""
        with WiseOpener(input_file) as fh:
            self._data = fh.readlines()

    def read_online(self, input_file):
        """Yield items from data online from input_file."""
        with WiseOpener(input_file) as fh:
            for item in fh:
                yield item

    def read_from_db(self, db_cursor):
        """Read data from database cursor."""
        for item in db_cursor:
            yield item

    def read_from_mongodb(self, table, query):
        """Read data online from mongodb."""
        cursor = table.find(query)
        n = cursor.count(query)
        start = 0
        limit = 2000
        end = start + limit
        while True:
            for x in cursor[start:end]:
                yield x
            start = end
            end += limit
            if start > n:
                break

    def update_mongodb(self, table, what, wherewith):
        table.update(what, wherewith, False, True)

    def write_to_file(self, output_file):
        """Write data to given output_file."""
        with WiseOpener(output_file, "w") as fh:
            fh.writelines(self._data)

    def write_to_db(self, db_cursor):
        """Write data with given database cursor."""
        raise NotImplementedError

    def write_to_mongodb(self, table, item):
        table.insert(item)

    def read_as_iter(self, source):
        """Read data from iterable source."""
        for item in source:
            self._data.append(item)

    def iterate(self, skip_empty=True):
        """Iterate over data."""
        if skip_empty:
            for item in self._data:
                if not item:
                    continue
                yield item
        else:
            for item in self._data:
                yield item

    def iterate_with_func(self, pre_func, iter_func):
        """Iterate over data with given iter_func.
        And data can be preprocessed with pre_func."""
        self._data = pre_func(self._data)
        for item in iter_func(self._data):
            yield item

    def do(self, cf, **args):
        """Do something with data with given core function and args.
        And get a result of doing.
        """
        result = cf(self._data, **args)
        return result

    def process(self, cf, **args):
        """Process data with given core function."""
        self._data = cf(self._data, **args)

    def clear(self):
        """Remove data."""
        self._data = None

    def do_with_iter(self, cf, **args):
        """Do something by iterating over data with given core function and args.
        And get a list of results of doing.
        """
        result = []
        for item in self._data:
            result.append(cf(item, **args))
        return result

    def process_with_iter(self, cf, **args):
        """Process by iterating over data with given core function."""
        for i, item in enumerate(self._data):
            self._data[i] = cf(item, **args)

    def sort(self, sort_func, reverse=False):
        """Sort data with sort_func and reversed param."""
        if not callable(sort_func):
            raise ConfigurationError(
                f"sort_func must be callable, got {type(sort_func).__name__}. "
                f"Provide a function that takes a data item and returns a sortable key. "
                f"Example: lambda x: x.some_field"
            )
        self._data.sort(key=sort_func, reverse=reverse)

    @property
    def data(self):
        return self._data

    @property
    def N(self):
        return len(self._data)


class AbstractFolderIO(object):
    """Abstract class for working with abstract data in folder.

    Public methods:

    - __init__(self, folder, mask=None)
    - iter_files(self)
    - get_files(self)
    - iter_filenames(self)
    - get_filenames(self)
    - iter_file_content(self)
    - copy_files_by_mask(self, dist_folder)

    >>> folder_reader = AbstractFolderIO(folder, mask=".")


    """

    def __init__(self, folder, mask="."):
        self.folder = folder
        self.mask = mask

    def iter_files(self):
        """iter over files in folder. Return file name."""
        for root, dirs, files in os.walk(self.folder, topdown=False):
            for name in files:
                if re.search(self.mask, name):
                    yield name

    def iter_folders(self):
        """iter over folders in folder. Return folder name."""
        for root, dirs, files in os.walk(self.folder, topdown=False):
            for folder in dirs:
                if re.search(self.mask, folder):
                    yield folder

    def get_files(self):
        """Get files in folder. Return file name."""
        result = []
        for root, dirs, files in os.walk(self.folder, topdown=False):
            for name in files:
                if re.search(self.mask, name):
                    result.append(name)
        return result

    def iter_filenames(self):
        """iter over files in folder. Return file name path."""
        for root, dirs, files in os.walk(self.folder, topdown=False):
            for name in files:
                if re.search(self.mask, name):
                    apath = os.path.join(root, name)
                    yield apath

    def get_filenames(self):
        """Get files in folder. Return path."""
        result = []
        for root, dirs, files in os.walk(self.folder, topdown=False):
            for name in files:
                if re.search(self.mask, name):
                    path = os.path.join(root, name)
                    result.append(path)
        return result

    def iter_path_names(self):
        """iter over files in folder. Return file name and path."""
        for root, dirs, files in os.walk(self.folder, topdown=False):
            for name in files:
                if re.search(self.mask, name):
                    apath = os.path.join(root, name)
                    yield name, apath

    def iter_file_content(self):
        """iter over files in folder. Return file content."""
        for root, dirs, files in os.walk(self.folder, topdown=False):
            for name in files:
                if re.search(self.mask, name):
                    path = os.path.join(root, name)
                    with WiseOpener(path, "rb") as fh:
                        yield fh.read()

    def iter_file_content_and_names(self):
        """
        Iterate over files in folder. Return file content, file_name, file_path.
        """
        for root, dirs, files in os.walk(self.folder, topdown=False):
            for name in files:
                if re.search(self.mask, name):
                    path = os.path.join(root, name)
                    with WiseOpener(path, "rb") as fh:
                        yield fh.read(), name, path

    def move_files_by_mask(self, dist_folder):
        """
        Move files matching the mask to the destination folder.

        Uses atomic operations to prevent data loss:
        - os.replace() for same-filesystem atomic replacement
        - Falls back to shutil.move() for cross-filesystem moves

        Args:
            dist_folder: Destination folder path

        Raises:
            OSError: If file move fails due to permissions or disk issues
            PermissionError: If insufficient permissions to move files
        """
        for file_path in self.iter_filenames():
            dist_file = os.path.join(dist_folder, os.path.split(file_path)[-1])
            try:
                # Use os.replace() for atomic file replacement (Python 3.3+)
                # This prevents data loss if the operation fails mid-way
                os.replace(file_path, dist_file)
            except OSError as e:
                # os.replace() may fail for cross-filesystem moves
                # Fall back to shutil.move() which handles this case
                try:
                    shutil.move(file_path, dist_file)
                    logger.warning(f"Used shutil.move() for cross-filesystem move: {file_path} -> {dist_file}")
                except (OSError, PermissionError) as move_error:
                    logger.error(f"Failed to move file {file_path} to {dist_file}: {move_error}")
                    raise

    def copy_files_by_mask(self, dist_folder):
        """
        Copy files matching the mask to the destination folder.

        Uses shutil.copy2() which preserves metadata and safely overwrites
        existing files without the risk of data loss.

        Args:
            dist_folder: Destination folder path

        Raises:
            OSError: If file copy fails due to permissions or disk issues
            PermissionError: If insufficient permissions to copy files
        """
        for file_path in self.iter_filenames():
            dist_file = os.path.join(dist_folder, os.path.split(file_path)[-1])
            try:
                # shutil.copy2() preserves metadata and safely overwrites destination if it exists
                shutil.copy2(file_path, dist_file)
            except (OSError, PermissionError) as e:
                logger.error(f"Failed to copy file {file_path} to {dist_file}: {e}")
                raise


def sc_iter_filepath_folder(folder, mask="."):
    """Shortcut for iterating file path in given folder."""
    reader = AbstractFolderIO(folder, mask=mask)
    for path in reader.iter_filenames():
        yield path


def sc_iter_filename_folder(folder, mask="."):
    """Shortcut for iterating filename in given folder."""
    reader = AbstractFolderIO(folder, mask=mask)
    for filename in reader.iter_files():
        yield filename


def sc_iter_folders(folder, mask="."):
    """Shortcut for iterating folders in given folder."""
    reader = AbstractFolderIO(folder, mask=mask)
    for folder in reader.iter_folders():
        yield folder


def sc_iter_path_name_folder(folder, mask="."):
    """Shortcut for iterating (filename, path) in given folder."""
    reader = AbstractFolderIO(folder, mask=mask)
    for filename, path in reader.iter_path_names():
        yield filename, path


def sc_iter_filedata_folder(folder, mask="."):
    """Shortcut for iterating file content in given folder."""
    reader = AbstractFolderIO(folder, mask=mask)
    for data in reader.iter_file_content():
        yield data


def sc_move_files(folder, dist_folder, mask="."):
    """Shortcut for moving file from folder to dist."""
    reader = AbstractFolderIO(folder, mask=mask)
    reader.move_files_by_mask(dist_folder)


def sc_process_file(file_name, cf, args_dict):
    """Shortcut for processing file
    with given cf funciton."""
    reader = AbstractFileIO()
    reader.read(file_name)
    args_dict["name"] = file_name
    text = cf(text, **args_dict)
    with WiseOpener(file_name, "w") as fh:
        fh.write(text)


def sc_process_folder(folder, cf, args_dict, mask="."):
    """Shortcut for processing each file in folder
    with given cf funciton."""
    reader = AbstractFolderIO(folder, mask=mask)
    for text, name, file_name in reader.iter_file_content_and_names():
        args_dict["name"] = name
        text = cf(text, **args_dict)
        with WiseOpener(file_name, "w") as fh:
            fh.write(text)


def sc_process_folder_to_other(
    folder, output_folder, cf, args_dict, mask=".", verbose=False
):
    """Shortcut for processing each file in folder
    with given cf funciton.

    To print names set *verbose* to True.
    """
    if not callable(cf):
        raise ConfigurationError(
            f"cf (core function) must be callable, got {type(cf).__name__}. "
            f"Provide a function that processes file content and returns modified content. "
            f"Example: lambda text, **args: process(text)"
        )
    reader = AbstractFolderIO(folder, mask=mask)
    for text, name, file in reader.iter_file_content_and_names():
        args_dict["name"] = name
        text = cf(text, **args_dict)
        output_file = os.path.join(output_folder, name)
        with WiseOpener(output_file, "w") as fh:
            fh.write(text)


def read_pickle_file(pickle_file):
    """Read pickle file and retrun its content."""
    with WiseOpener(pickle_file, "r") as fh:
        data = pickle.load(fh)
    return data
