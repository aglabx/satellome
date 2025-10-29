#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 04.06.2011
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com
"""
    Classes:
    
    - WiseOpener(object)
    - AbstractFileIO(object)
    - AbstractFolderIO(object)
    
    Shortcuts:
    
    - sc_iter_filepath_folder(folder, mask=".")
    - sc_iter_filename_folder(folder, mask=".")
    - sc_iter_filedata_folder(folder, mask=".")
    - sc_move_files(folder, dist_folder, mask=".")

"""
import bz2
import gzip
import logging
import os
import pickle
import re
import shutil

logger = logging.getLogger(__name__)


class WiseOpener(object):
    """Opener to open usual files and gzip or bzip archives."""

    def __init__(self, file_name, mode=None):
        self.file_name = file_name
        if not mode:
            mode = "r"
        if not mode in ["w", "r", "a", "wb", "rb", "ab"]:
            logger.error("Wrong file mode: %s" % mode)
            raise Exception("Wrong file mode: %s" % mode)
        self.mode = mode
        self.fh = None

    def __enter__(self):
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
        assert hasattr(sort_func, "__call__")
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
        for file_path in self.iter_filenames():
            dist_file = os.path.join(dist_folder, os.path.split(file_path)[-1])
            if os.path.isfile(dist_file):
                os.remove(dist_file)
                # TODO: fix me
            os.rename(file_path, dist_file)

    def copy_files_by_mask(self, dist_folder):
        for file_path in self.iter_filenames():
            dist_file = os.path.join(dist_folder, os.path.split(file_path)[-1])
            if os.path.isfile(dist_file):
                os.remove(dist_file)
            shutil.copy2(file_path, dist_file)


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
    assert hasattr(cf, "__call__")
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
