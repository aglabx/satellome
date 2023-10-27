#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 05.06.2011
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com
"""
Classes:
    
- TabDelimitedFileIO(AbstractFileIO)
    
"""
import csv
import os
import tempfile

from PyExp import AbstractFileIO

csv.field_size_limit(1000000000)


class TabDelimitedFileIO(AbstractFileIO):
    """Working with tab delimited file.

    Public methods:

    - sort(self, sort_func, reversed=False)

    Inherited public properties:

    - data  - iterable data, each item is tuple (head, body)
    - N     - a number of items in data

    Overrided public methods:

    - __init__(self, skip_first=False, format_func=None, delimeter="\\t", skip_startswith=None)
    - write_to_file(self, output_file)
    - read_from_file(self, input_file)
    - read_online(self, input_file) ~> item

    Inherited public methods:

    - read_from_db(self, db_cursor)
    - write_to_db(self, db_cursor)
    - read_as_iter(self, source)
    - iterate(self) ~> item of data
    - do(self, cf, args) -> result
    - process(self, cf, args)
    - clear(self)
    - do_with_iter(self, cf, args) -> [result,]
    - process_with_iter(self, cf, args)

    """

    def __init__(
        self, skip_first=False, format_func=None, delimeter="\t", skip_startswith="#"
    ):
        """Overrided. Set token velue.

        Keyword arguments:

        - skip_first   -- skip first line in file
        - format_funcs -- list of functions that format corresponding item in tab delimited line
        - delimeter    -- line delimeter
        - skip_startswith -- skip lines starts with this value
        """
        super(TabDelimitedFileIO, self).__init__()

        self.skip_first = skip_first
        self.format_func = format_func
        self.delimeter = delimeter
        self.skip_startswith = skip_startswith

    def read_from_file(self, input_file):
        """OVerrided. Read data from tab delimeted input_file."""
        with open(input_file) as fh:
            self._data = fh.readlines()
        if self.skip_first:
            self._data.pop(0)
        if self.skip_startswith:
            self._data = [
                line for line in self._data if not line.startswith(self.skip_startswith)
            ]
        self._data = [self._process_tab_delimeited_line(line) for line in self._data]

    def read_online(self, input_file):
        """Overrided. Yield items online from data from input_file."""
        with open(input_file) as fh:
            for i, line in enumerate(fh):
                if self.skip_first and i == 0:
                    continue
                if self.skip_startswith and line.startswith(self.skip_startswith):
                    continue
                yield self._process_tab_delimeited_line(line)

    def _process_tab_delimeited_line(self, line):
        """Format line with format_func."""

        line = line.strip().split(self.delimeter)
        if self.format_func:
            assert hasattr(self.format_func, "__call__")
            line = self.format_func(line)
        return line

    def _all_str(self, line):
        """Convert to string all items in line."""
        return [str(x) for x in line]

    def write_to_file(self, output_file):
        """Overrided. Write data to tab delimited output_file."""
        self._data = ["\t".join(self._all_str(line)) for line in self._data]
        with open(output_file, "w") as fh:
            fh.writelines(self._data)


def sc_iter_tab_file(
    input_file,
    data_type,
    skip_starts_with="#",
    remove_starts_with=None,
    preprocess_function=None,
    check_function=None,
):
    """Iter over tab file, yield an object of given data_type."""

    temp_file = tempfile.NamedTemporaryFile(delete=False)
    temp_file_name = temp_file.name

    if remove_starts_with:
        with open(input_file, "r") as fh:
            data = fh.readlines()
        data = [x for x in data if not x.startswith(remove_starts_with)]
        with open(temp_file_name, "w") as fh:
            fh.writelines(data)
        input_file = temp_file_name
    if preprocess_function:
        with open(input_file, "r") as fh:
            data = fh.readlines()
        data = [preprocess_function(x) for x in data]
        with open(temp_file_name, "w") as fh:
            fh.writelines(data)
        input_file = temp_file_name
    if check_function:
        with open(input_file, "r") as fh:
            data = fh.readlines()
        data = [x for x in data if check_function(x)]
        with open(temp_file_name, "w") as fh:
            fh.writelines(data)
        input_file = temp_file_name
    with open(input_file) as fh:
        fields = data_type().dumpable_attributes
        for data in csv.DictReader(
            fh, fieldnames=fields, delimiter="\t", quoting=csv.QUOTE_NONE
        ):
            if skip_starts_with:
                if data[fields[0]].startswith(skip_starts_with):
                    continue
            obj = data_type()
            obj.set_with_dict(data)
            yield obj
    if os.path.isfile(temp_file_name):
        os.unlink(temp_file_name)


def sc_iter_simple_tab_file(input_file):
    """Iter tab file, yield a list."""
    with open(input_file) as fh:
        for data in csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
            yield data


def sc_read_dictionary(dict_file, value_func=None):
    """Read file of tab-dilimited pairs.
    key\tvalue
    """
    result = {}
    with open(dict_file) as fh:
        for data in csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
            if hasattr(value_func, "__call__"):
                data[1] = value_func(data[1])
            result[data[0]] = data[1]
    return result


def sc_write_model_to_tab_file(output_file, objs):
    """Write model obj to tab-delimited file."""
    with open(output_file, "w") as fh:
        for obj in objs:
            fh.write(str(obj))


def sc_read_simple_tab_file(input_file, skip_first=False):
    """Iter tab file, yield a list."""
    result = []
    with open(input_file) as fh:
        if skip_first:
            fh.readline()
        for data in csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
            result.append(data)
    return result
