#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 05.06.2011
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com
"""
Tab-delimited file I/O utilities for genomic data.

Provides readers, writers, and iterators for tab-delimited files containing
genomic annotations, tandem repeats, and other structured data. Supports
CSV-based parsing with custom delimiters, field mapping to data models,
and various preprocessing options.

Classes:
    TabDelimitedFileIO: General-purpose tab-delimited file reader/writer
                       with sorting and formatting capabilities

Functions:
    sc_iter_tab_file: Iterate over tab file yielding model objects with
                     optional preprocessing and filtering
    sc_iter_simple_tab_file: Simple tab file iterator yielding raw lists
    sc_read_dictionary: Read two-column tab file into dictionary
    sc_write_model_to_tab_file: Write model objects to tab file
    sc_read_simple_tab_file: Read entire tab file into list of lists

Key Features:
    - Memory-efficient streaming for large files
    - Automatic field mapping to data model attributes
    - Flexible preprocessing and filtering pipelines
    - Comment line skipping (lines starting with #)
    - Custom delimiters and field formatting
    - Integration with AbstractModel-based classes

Example:
    >>> # Read TRF data as model objects
    >>> for trf_obj in sc_iter_tab_file("repeats.tab", TRModel):
    ...     print(f"{trf_obj.trf_consensus}: {trf_obj.trf_l_ind}-{trf_obj.trf_r_ind}")
    >>>
    >>> # Read simple tab file
    >>> for fields in sc_iter_simple_tab_file("data.tab"):
    ...     print(fields[0], fields[1])
    >>>
    >>> # Read dictionary from two-column file
    >>> name2id = sc_read_dictionary("names.tab", value_func=int)

See Also:
    satellome.core_functions.io.abstract_reader: Base I/O classes
    satellome.core_functions.models.abstract_model: Data model base class
"""
import csv
import os
import tempfile

from satellome.core_functions.exceptions import ConfigurationError
from satellome.core_functions.io.abstract_reader import AbstractFileIO

csv.field_size_limit(1000000000)


class TabDelimitedFileIO(AbstractFileIO):
    """
    General-purpose reader/writer for tab-delimited files.

    Provides flexible parsing and writing of tab-delimited data with support
    for custom delimiters, header skipping, comment line filtering, and
    field-level formatting. Extends AbstractFileIO with tab-specific methods.

    Attributes:
        skip_first (bool): Skip first line (header) when reading
        format_func (callable): Optional function to format each parsed line
        delimeter (str): Field delimiter character (default: tab)
        skip_startswith (str): Skip lines starting with this string (e.g., "#")
        data (list): Parsed data rows (inherited)
        N (int): Number of data rows (inherited)

    Public Methods:
        sort: Sort data using custom function
        read_from_file: Load entire file into memory
        read_online: Stream file line-by-line (memory efficient)
        write_to_file: Write data to tab-delimited file

    Private Methods:
        _process_tab_delimeited_line: Parse and format individual line
        _all_str: Convert all fields in a line to strings

    Inherited Methods:
        read_from_db, write_to_db, read_as_iter, iterate, clear,
        do, process, do_with_iter, process_with_iter
        (see AbstractFileIO for details)

    Example:
        >>> # Read with custom formatting
        >>> reader = TabDelimitedFileIO(
        ...     skip_first=True,
        ...     format_func=lambda fields: [f.strip().upper() for f in fields],
        ...     skip_startswith="#"
        ... )
        >>> for line in reader.read_online("data.tab"):
        ...     print(line[0], line[1])
        >>>
        >>> # Write data
        >>> writer = TabDelimitedFileIO()
        >>> writer.data = [["col1", "col2"], ["val1", "val2"]]
        >>> writer.write_to_file("output.tab")

    Note:
        - For model object iteration, consider using sc_iter_tab_file()
        - CSV field size limit set to 1 billion characters (module level)
        - format_func receives list of fields and must return list
    """

    def __init__(
        self, skip_first=False, format_func=None, delimeter="\t", skip_startswith="#"
    ):
        """
        Initialize tab-delimited file reader/writer with parsing options.

        Args:
            skip_first (bool, optional): Skip first line (typically header).
                                        Defaults to False.
            format_func (callable, optional): Function to format each parsed line.
                                             Receives list of fields, returns list.
                                             Defaults to None (no formatting).
            delimeter (str, optional): Field delimiter character.
                                      Defaults to "\\t" (tab).
            skip_startswith (str, optional): Skip lines starting with this string.
                                            Common value: "#" for comments.
                                            Defaults to "#".

        Raises:
            ConfigurationError: If format_func is not callable (checked during parsing)

        Example:
            >>> # Skip header and comments, uppercase all fields
            >>> reader = TabDelimitedFileIO(
            ...     skip_first=True,
            ...     format_func=lambda fields: [f.upper() for f in fields],
            ...     skip_startswith="#"
            ... )
        """
        super(TabDelimitedFileIO, self).__init__()

        self.skip_first = skip_first
        self.format_func = format_func
        self.delimeter = delimeter
        self.skip_startswith = skip_startswith

    def read_from_file(self, input_file):
        """
        Load entire tab-delimited file into memory.

        Reads all lines from file, applies skip_first and skip_startswith
        filters, parses each line into fields, and stores in self.data.

        Args:
            input_file (str): Path to input tab-delimited file

        Note:
            - Not memory-efficient for large files; consider read_online()
            - Processes all lines through _process_tab_delimeited_line()
            - Result stored in self._data (list of field lists)
        """
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
        """
        Stream tab-delimited file line-by-line (memory efficient).

        Yields parsed lines one at a time without loading entire file into
        memory. Applies skip_first and skip_startswith filters per line.

        Args:
            input_file (str): Path to input tab-delimited file

        Yields:
            list: Parsed and formatted fields from each line

        Example:
            >>> reader = TabDelimitedFileIO(skip_first=True)
            >>> for fields in reader.read_online("large_file.tab"):
            ...     process(fields)  # Process one line at a time

        Note:
            - Preferred for large files (GB+) to avoid memory issues
            - Each line processed through _process_tab_delimeited_line()
        """
        with open(input_file) as fh:
            for i, line in enumerate(fh):
                if self.skip_first and i == 0:
                    continue
                if self.skip_startswith and line.startswith(self.skip_startswith):
                    continue
                yield self._process_tab_delimeited_line(line)

    def _process_tab_delimeited_line(self, line):
        """
        Parse and optionally format a single line.

        Strips whitespace, splits by delimiter, and applies format_func if set.

        Args:
            line (str): Raw line from file

        Returns:
            list: Parsed (and optionally formatted) field list

        Raises:
            ConfigurationError: If format_func is not callable

        Note:
            - Called internally by read_from_file() and read_online()
            - format_func must accept and return a list
        """
        line = line.strip().split(self.delimeter)
        if self.format_func:
            if not callable(self.format_func):
                raise ConfigurationError(
                    f"format_func must be callable, got {type(self.format_func).__name__}. "
                    f"Provide a function that takes a list of tab-separated fields and returns processed list. "
                    f"Example: lambda fields: [field.strip() for field in fields]"
                )
            line = self.format_func(line)
        return line

    def _all_str(self, line):
        """
        Convert all fields in a line to strings.

        Args:
            line (list): List of fields (any types)

        Returns:
            list: List of stringified fields

        Note:
            - Called by write_to_file() to ensure tab-delimited output compatibility
        """
        return [str(x) for x in line]

    def write_to_file(self, output_file):
        """
        Write self.data to tab-delimited file.

        Converts all fields to strings, joins with tabs, and writes to file.

        Args:
            output_file (str): Path to output file

        Note:
            - Overwrites existing file
            - All fields converted to strings via _all_str()
            - Lines joined with tabs (hardcoded, ignores self.delimeter)
        """
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
    """
    Iterate over tab-delimited file yielding data model objects.

    Flexible streaming parser that maps tab-delimited rows to AbstractModel
    objects with optional preprocessing and filtering pipeline. Uses temporary
    files for multi-stage processing.

    Processing Pipeline:
        1. Remove lines (remove_starts_with) → temp file
        2. Preprocess lines (preprocess_function) → temp file
        3. Filter lines (check_function) → temp file
        4. Parse to model objects via CSV DictReader
        5. Skip comment lines (skip_starts_with)
        6. Yield model objects

    Args:
        input_file (str): Path to input tab-delimited file
        data_type (type): AbstractModel subclass (e.g., TRModel, Gff3Model)
                         Must have dumpable_attributes and set_with_dict()
        skip_starts_with (str, optional): Skip lines starting with this string
                                         after all preprocessing. Defaults to "#".
        remove_starts_with (str, optional): Remove lines starting with this string
                                           before preprocessing. Defaults to None.
        preprocess_function (callable, optional): Function to transform each line
                                                 (receives str, returns str).
                                                 Defaults to None.
        check_function (callable, optional): Function to filter lines (receives str,
                                            returns bool; True=keep). Defaults to None.

    Yields:
        AbstractModel: Model objects with fields populated from tab-delimited data

    Example:
        >>> # Basic usage
        >>> for trf in sc_iter_tab_file("repeats.tab", TRModel):
        ...     print(f"{trf.trf_consensus}: {trf.trf_l_ind}-{trf.trf_r_ind}")
        >>>
        >>> # With preprocessing and filtering
        >>> def preprocess(line):
        ...     return line.replace("chr", "").upper()
        >>> def check(line):
        ...     return not line.startswith("chrM")  # Skip mitochondrial
        >>> for trf in sc_iter_tab_file(
        ...     "repeats.tab",
        ...     TRModel,
        ...     preprocess_function=preprocess,
        ...     check_function=check
        ... ):
        ...     print(trf.trf_head)

    Note:
        - Uses temporary files for multi-stage processing (cleaned up automatically)
        - CSV field size limit: 1 billion characters (module level setting)
        - Field names extracted from data_type().dumpable_attributes
        - Preprocessing happens line-by-line; each stage creates new temp file
        - All stages are optional (can use just skip_starts_with for simple cases)
    """

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
            # Skip header row if present (first field value equals first field name)
            if data[fields[0]] == fields[0]:
                continue
            obj = data_type()
            obj.set_with_dict(data)
            yield obj
    if os.path.isfile(temp_file_name):
        os.unlink(temp_file_name)


def sc_iter_simple_tab_file(input_file):
    """
    Simple streaming iterator for tab-delimited files.

    Memory-efficient line-by-line reader that yields raw field lists
    without model object conversion or preprocessing.

    Args:
        input_file (str): Path to input tab-delimited file

    Yields:
        list: Raw fields from each line (no type conversion)

    Example:
        >>> for fields in sc_iter_simple_tab_file("data.tab"):
        ...     chrom, start, end = fields[0], int(fields[1]), int(fields[2])
        ...     print(f"{chrom}:{start}-{end}")

    Note:
        - No header skipping, comment filtering, or preprocessing
        - For model objects, use sc_iter_tab_file()
        - For full data loading, use sc_read_simple_tab_file()
    """
    with open(input_file) as fh:
        for data in csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
            yield data


def sc_read_dictionary(dict_file, value_func=None):
    """
    Read two-column tab-delimited file into dictionary.

    Parses file with key-value pairs (one per line) into a dictionary,
    optionally applying transformation to values.

    Args:
        dict_file (str): Path to input file (format: key<tab>value per line)
        value_func (callable, optional): Function to transform values
                                        (e.g., int, float, str.upper).
                                        Defaults to None (strings).

    Returns:
        dict: Dictionary mapping keys to (optionally transformed) values

    Example:
        >>> # Read as strings
        >>> name2desc = sc_read_dictionary("gene_names.tab")
        >>> name2desc["TP53"]
        'Tumor protein p53'
        >>>
        >>> # Read with integer values
        >>> gene2count = sc_read_dictionary("counts.tab", value_func=int)
        >>> gene2count["TP53"]
        150

    Note:
        - File format: two tab-separated columns per line
        - First column: key, second column: value
        - No header line expected
    """
    result = {}
    with open(dict_file) as fh:
        for data in csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
            if hasattr(value_func, "__call__"):
                data[1] = value_func(data[1])
            result[data[0]] = data[1]
    return result


def sc_write_model_to_tab_file(output_file, objs):
    """
    Write model objects to tab-delimited file.

    Serializes AbstractModel objects to tab-delimited format using their
    __str__() method (which formats dumpable_attributes).

    Args:
        output_file (str): Path to output file
        objs (iterable): AbstractModel objects (e.g., list of TRModel)

    Example:
        >>> trf_objects = [TRModel(), TRModel()]
        >>> # ... populate objects ...
        >>> sc_write_model_to_tab_file("output.tab", trf_objects)

    Note:
        - Uses model's __str__() which formats according to dumpable_attributes
        - No header line written
        - For reading back, use sc_iter_tab_file()
    """
    with open(output_file, "w") as fh:
        for obj in objs:
            fh.write(str(obj))


def sc_read_simple_tab_file(input_file, skip_first=False):
    """
    Load entire tab-delimited file into memory as list of lists.

    Reads all lines from file into a list, optionally skipping header.
    Not memory-efficient for large files.

    Args:
        input_file (str): Path to input tab-delimited file
        skip_first (bool, optional): Skip first line (header). Defaults to False.

    Returns:
        list: List of field lists (one per line)

    Example:
        >>> # Read entire file
        >>> data = sc_read_simple_tab_file("small_file.tab", skip_first=True)
        >>> for fields in data:
        ...     print(fields[0], fields[1])
        >>> print(f"Total lines: {len(data)}")

    Note:
        - Loads entire file into memory
        - For large files, use sc_iter_simple_tab_file() instead
        - For model objects, use sc_iter_tab_file()
    """
    result = []
    with open(input_file) as fh:
        if skip_first:
            fh.readline()
        for data in csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
            result.append(data)
    return result
