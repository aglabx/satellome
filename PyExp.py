#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.06.2011
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import os
import shutil
import bz2
import gzip
import re
import logging

try:
    import simplejson
except:
    logging.warning("Install simplejson module")


class AbstractModel(object):
    """Сlass for data wrapping.


    Private methods:

    - __str__(self) used dumpable_attributes
    - set_with_dict(self, dictionary)
    Initiation. Create attributes accordong to
    public properties:

    - dumpable_attributes
    - int_attributes
    - float_attributes
    - list_attributes
    - list_attributes_types
    - other_attributes
    """

    dumpable_attributes = []
    int_attributes = []
    float_attributes = []
    list_attributes = []
    list_attributes_types = {}
    other_attributes = {}

    def __init__(self):
        """Create attributes accordong to
        - dumpable_attributes
        - int_attributes
        - float_attributes
        - list_attributes
        - list_attributes_types
        - other_attributes
        """
        for attr in self.dumpable_attributes:
            setattr(self, attr, None)
        for attr in self.int_attributes:
            setattr(self, attr, 0)
        for attr in self.float_attributes:
            setattr(self, attr, 0.0)
        for attr in self.other_attributes:
            setattr(self, attr, self.other_attributes[attr])

    def __str__(self):
        """Get string representation with fields
        defined in dumpable_attributes."""
        self.preprocess_data()
        result = []
        for attr in self.dumpable_attributes:
            data = getattr(self, attr)
            if attr in self.list_attributes:
                if data is None:
                    data = []
                data = ",".join([str(x) for x in data])
            try:
                result.append(str(data).strip())
            except UnicodeEncodeError:
                result.append(unicode(data).strip())
        result = "%s\n" % "\t".join(result)
        return result

    def print_human_friendly(self):
        """Print human friendly string representation with fields
        defined in dumpable_attributes."""
        self.preprocess_data()
        result = []
        largest_name_length = max([len(x) for x in self.dumpable_attributes])
        print_string = "{0:%s} => {1}" % largest_name_length
        for attr in self.dumpable_attributes:
            data = getattr(self, attr)
            if attr in self.list_attributes:
                if data is None:
                    data = []
                data = ",".join([str(x) for x in data])
            try:
                print(print_string.format(attr, data.strip()))
            except UnicodeEncodeError:
                data = unicode(data).strip()
                result.append()
                print(print_string.format(attr, data.strip()))

    def get_as_string(self, dumpable_attributes):
        """Get string representation with fields
        defined in dumpable_attributes."""
        return self.__str__()

    def set_with_dict(self, dictionary):
        """Set object with dictionaty."""
        for key, value in dictionary.items():
            key, value = self.preprocess_pair(key, value)
            try:
                if value == "None" or value is None:
                    value = None
                elif key in self.int_attributes:
                    value = int(value)
                elif key in self.float_attributes:
                    value = float(value)
                elif key in self.list_attributes:
                    if not value:
                        value = []
                        continue
                    value = value.split(",")
                    value = [self.list_attributes_types[key](x) for x in value]
                setattr(self, key, value)
            except ValueError as e:
                logging.warning(self.dumpable_attributes)
                logging.warning(dictionary.items())
                raise ValueError(e)
            except TypeError as e:
                logging.warning(self.dumpable_attributes)
                logging.warning(dictionary.items())
                raise TypeError(e)

    def set_with_list(self, data):
        """Set object with list."""
        n = len(data)
        dumpable_attributes = self.dumpable_attributes
        if n != len(self.dumpable_attributes):
            if (
                hasattr(self, "alt_dumpable_attributes")
                and len(self.alt_dumpable_attributes) == n
            ):
                dumpable_attributes = self.alt_dumpable_attributes
            else:
                logging.warning(data)
                raise Exception("Wrong number of fields in data.")
        for i, value in enumerate(data):
            key = dumpable_attributes[i]
            if value == "None":
                value = None
            elif key in self.int_attributes:
                value = int(value)
            elif key in self.float_attributes:
                value = float(value)
            elif key in self.list_attributes:
                value = value.split(",")
                value = [self.list_attributes_types[key](x) for x in value]
            setattr(self, key, value)

    def as_dict(self):
        """ """
        return self.get_as_dict()

    def get_as_dict(self):
        """Get dictionary representation with fields
        defined in dumpable_attributes"""
        self.preprocess_data()
        result = {}
        for attr in self.dumpable_attributes:
            result[attr] = getattr(self, attr)
        return result

    def get_as_json(self, preprocess_func=None):
        """Return JSON representation."""
        self.preprocess_data()
        d = self.get_as_dict()
        if preprocess_func:
            d = preprocess_func(d)
        return simplejson.dumps(d)

    def preprocess_data(self):
        """Any data preprocessing before returning."""
        pass

    def preprocess_pair(self, key, value):
        """Any data preprocessing before initiation from dictionary."""
        return key, value

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        return setattr(self, key, value)


class WiseOpener(object):
    """Opener to open usual files and gzip or bzip archives."""

    def __init__(self, file_name, mode=None):
        self.file_name = file_name
        if not mode:
            mode = "r"
        if not mode in ["w", "r", "a", "wb", "rb", "ab"]:
            logging.warning("Wrong file mode: %s" % mode)
            raise Exception("Wrong file mode: %s" % mode)
        self.mode = mode
        self.fh = None

    def __enter__(self):
        if self.file_name.endswith(".gz"):
            logging.info("Open as gz archive")
            if not "b" in self.mode:
                self.mode += "b"
            self.fh = gzip.open(self.file_name, self.mode)
        elif self.file_name.endswith(".bz2"):
            logging.info("Open as bz2 archive")
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
            logging.info("Move: ", file_path, dist_file)
            if os.path.isfile(dist_file):
                os.remove(dist_file)
                # TODO: fix me
            os.rename(file_path, dist_file)

    def copy_files_by_mask(self, dist_folder):
        for file_path in self.iter_filenames():
            dist_file = os.path.join(dist_folder, os.path.split(file_path)[-1])
            logging.info("Copy: ", file_path, dist_file)
            if os.path.isfile(dist_file):
                os.remove(dist_file)
            shutil.copy2(file_path, dist_file)


def sc_iter_filepath_folder(folder, mask="."):
    """Shortcut for iterating file path in given folder."""
    reader = AbstractFolderIO(folder, mask=mask)
    for path in reader.iter_filenames():
        yield path
