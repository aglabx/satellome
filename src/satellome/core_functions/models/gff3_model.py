#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 28.10.2014
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com
"""
GFF3 (Generic Feature Format version 3) data models and I/O handlers.

This module provides classes for parsing, manipulating, and writing GFF3 files,
which are commonly used for genome annotations. GFF3 is a tab-delimited format
with 9 columns: seqid, source, type, start, end, score, strand, phase, attributes.

Classes:
    Gff3Model: Data model for a single GFF3 record
    Gff3FeatureDict: Dictionary for storing GFF3 feature attributes
    Gff3FileIO: File I/O handler for reading/writing GFF3 files
"""


import collections
import csv
import logging
import sys

from satellome.core_functions.models.abstract_model import AbstractModel

logger = logging.getLogger(__name__)

from satellome.core_functions.io.tab_file import TabDelimitedFileIO

if sys.version_info.major == 3 and sys.version_info.minor >= 10:
    from collections.abc import MutableMapping
else:
    from collections import MutableMapping


class Gff3Model(AbstractModel):
    """Class for gff3 data wrapping."""

    dumpable_attributes = [
        "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes",
    ]

    int_attributes = [
        "start",
        "end",
    ]

    @property
    def target(self):
        return self.seqid

    @property
    def chrm(self):
        return self.seqid

    @property
    def chromosome(self):
        return self.seqid

    @property
    def contig(self):
        return self.seqid

    @property
    def length(self):
        return abs(self.end - self.start)

    def get_coordinates(self):
        """
        Get normalized genomic coordinates as a tuple.

        Returns coordinates in ascending order regardless of strand direction.
        Handles cases where end < start (reverse strand features).

        Returns:
            tuple: (start, end) coordinates in ascending order

        Example:
            >>> record = Gff3Model()
            >>> record.start, record.end = 1000, 2000
            >>> record.get_coordinates()
            (1000, 2000)
            >>> record.start, record.end = 2000, 1000  # Reverse strand
            >>> record.get_coordinates()
            (1000, 2000)
        """
        if self.end < self.start:
            return (self.end, self.start)
        return (self.start, self.end)

    def save_original(self, line):
        """
        Store the original GFF3 line for reference.

        Args:
            line (str): Original GFF3 line from file
        """
        self.original = line

    def as_gff3(self):
        """
        Convert the GFF3 record to a properly formatted GFF3 line.

        Serializes all dumpable attributes into a tab-delimited GFF3 format string.
        Handles both simple and dict-type attributes in the attributes field.
        Dict-type attributes are formatted as key:value pairs separated by colons.

        Returns:
            str: GFF3-formatted line ending with newline character

        Example:
            >>> record = Gff3Model()
            >>> record.seqid = "chr1"
            >>> record.source = "test"
            >>> record.type = "gene"
            >>> record.start, record.end = 1000, 2000
            >>> record.score, record.strand, record.phase = ".", "+", "."
            >>> record.attributes = {"ID": "gene001", "Name": "test_gene"}
            >>> line = record.as_gff3()
            >>> "chr1\\ttest\\tgene\\t1000\\t2000" in line
            True
        """
        s = []
        for attr in self.dumpable_attributes:
            if not attr == "attributes":
                s.append(str(getattr(self, attr)))
        if not hasattr(self, "attributes"):
            self.attributes = {}
        if self.attributes:
            attr_keys = [
                key
                for key in list(self.attributes.keys())
                if not isinstance(self.attributes[key], dict)
            ]
            dict_attr_keys = [
                key
                for key in list(self.attributes.keys())
                if isinstance(self.attributes[key], dict)
            ]
            attr_keys.sort()
            attr = []
            for k in attr_keys:
                attr.append("%s=%s" % (k, self.attributes[k]))
            for k in dict_attr_keys:
                data = []
                for key in self.attributes[k]:
                    data.append(f"{key}:{self.attributes[k][key]}")
                attr.append("%s=%s" % (k, ",".join(data)))
            self.raw_features = ";".join(attr)
        s.append(self.raw_features)
        s = "\t".join(s)
        return "%s\n" % s


class Gff3FeatureDict(MutableMapping):
    """A dictionary for gff3 features"""

    def __init__(self, *args, **kwargs):
        self.store = dict()
        self.update(dict(*args, **kwargs))  # use the free update to set keys

    def __getitem__(self, key):
        return self.store[self.__keytransform__(key)]

    def __setitem__(self, key, value):
        self.store[self.__keytransform__(key)] = value

    def __delitem__(self, key):
        del self.store[self.__keytransform__(key)]

    def __iter__(self):
        return iter(self.store)

    def __len__(self):
        return len(self.store)

    def __keytransform__(self, key):
        return key


class Gff3FileIO(TabDelimitedFileIO):
    """
    File I/O handler for reading and writing GFF3 format files.

    Extends TabDelimitedFileIO with GFF3-specific parsing logic:
    - Preserves header/comment lines starting with '#'
    - Parses the 9-column GFF3 format
    - Handles complex attribute fields with key=value pairs
    - Supports filtering by feature type

    Attributes:
        headers (list): List of GFF3 header/comment lines from the file
    """

    def __init__(self, *args, **kwargs):
        """
        Initialize GFF3 file I/O handler.

        Args:
            *args: Positional arguments passed to parent TabDelimitedFileIO
            **kwargs: Keyword arguments passed to parent TabDelimitedFileIO
        """
        super(TabDelimitedFileIO, self).__init__(*args, **kwargs)
        self.headers = []

    def read_online(self, file_name, only_fields=None):
        """
        Stream GFF3 records from file one at a time (memory-efficient).

        Reads GFF3 file line by line, yielding Gff3Model objects without loading
        the entire file into memory. Preserves header lines and parses attributes
        field into structured dictionaries. Handles special cases like Dbxref
        attributes with multiple colon-separated values.

        Args:
            file_name (str): Path to GFF3 file
            only_fields (list, optional): List of feature types to include (e.g., ["gene", "CDS"]).
                                         If None, returns all features. Defaults to None.

        Yields:
            Gff3Model: Parsed GFF3 record object

        Example:
            >>> io_handler = Gff3FileIO()
            >>> for record in io_handler.read_online("annotations.gff3", only_fields=["gene"]):
            ...     print(record.seqid, record.start, record.end)
        """

        def skip_comments(iterable):
            for line in iterable:
                if not line.startswith("#"):
                    yield line

        with open(file_name) as fh:
            for i, line in enumerate(fh):
                if line.startswith("#"):
                    self.headers.append(line)
                else:
                    break

        fields = Gff3Model().dumpable_attributes

        with open(file_name) as fh:
            for data in csv.DictReader(
                skip_comments(fh),
                fieldnames=fields,
                delimiter="\t",
                quoting=csv.QUOTE_NONE,
            ):
                if only_fields and data["type"] not in only_fields:
                    continue
                _features = {}
                if data["attributes"]:
                    data["raw_features"] = data["attributes"]
                    for item in data["attributes"].split(";"):
                        if not item.strip():
                            continue
                        k, v = item.strip().split("=")
                        if k == "Dbxref":
                            try:
                                v = dict(
                                    [
                                        (
                                            ref.split(":")[0],
                                            ":".join(ref.split(":")[1:]),
                                        )
                                        for ref in v.split(",")
                                    ]
                                )
                            except (ValueError, IndexError, AttributeError) as e:
                                logger.error(f"Error parsing GFF3 attribute value: {v}, error: {e}")
                        _features[k] = v
                data["attributes"] = _features
                obj = Gff3Model()
                try:
                    obj.set_with_dict(data)
                except (AttributeError, KeyError, TypeError) as e:
                    logger.error(f"Can't parse features for {data}: {e}")
                    continue
                yield obj
