#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 05.06.2011
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com


import os
import tempfile
from sequence_model import SequenceModel
from gff_model import Gff3Model
import csv
from PyExp import AbstractFileIO


class AbstractBlockFileIO(AbstractFileIO):
    """Working with file with data organized in block, where each block starts with same token.

    Public methods:
    - get_block_sequence(self, head_start, next_head, fh)
    - get_blocks(self, token, fh)
    - gen_block_sequences(self, token, fh)

    Inherited public properties:

    - data  - iterable data, each item is tuple (head, body)
    - N     - a number of items in data

    Overrided public methods:

    - __init__(self, token)
    - read_from_file(self, input_file)
    - read_online(self, input_file) ~> item

    Inherited public methods:

    - read_from_db(self, db_cursor)
    - write_to_file(self, output_file)
    - write_to_db(self, db_cursor)
    - read_as_iter(self, source)
    - iterate(self) ~> item of data
    - do(self, cf, args) -> result
    - process(self, cf, args)
    - clear(self)
    - do_with_iter(self, cf, args) -> [result,]
    - process_with_iter(self, cf, args)

    """

    def __init__(self, token, **args):
        """Overrided. Set token velue."""
        super(AbstractBlockFileIO, self).__init__(**args)
        self.token = token
        self.wise_opener = self.get_opener()

    def read_from_file(self, input_file):
        """Overrided. Read data from given input_file."""
        with self.wise_opener(input_file, "r") as fh:
            for head, body, start, next in self.gen_block_sequences(self.token, fh):
                self.data.append((head, body, start, next))

    def read_online(self, input_file):
        """Overrided. Yield items from data online from input_file."""
        with self.wise_opener(input_file, "r") as fh:
            for head, body, start, next in self.gen_block_sequences(self.token, fh):
                yield (head, body, start, next)

    def get_block_sequence(self, head_start, next_head, fh):
        """Get a data block (head, seq, head_start, head_end).

        Arguments:

        - head_start -- a head starting position in a file
        - next_head  -- a next to head starting position in a file
        - fh         -- an open file handler

        Return format:

        - head       -- a block head
        - seq        -- a block body
        - head_start -- a file pointer to block start
        - head_end   -- a file pointer to next block start or 0
        """
        head_start = int(head_start)
        next_head = int(next_head)
        sequence = ""
        fh.seek(head_start)
        pos = head_start
        if next_head:
            head = fh.readline()
            pos = fh.tell()
            while pos != next_head:
                temp_seq = fh.readline()
                pos = fh.tell()
                if pos != next_head:
                    sequence += temp_seq
            sequence += temp_seq
        else:
            head = fh.readline()
            fasta_list = fh.readlines()
            sequence = "".join(fasta_list)
        return (head, sequence, head_start, next_head)

    def get_blocks(self, token, fh):
        """Get a list of the token positions in given file (first, next).
        For the last string the function returns (last string, 0).

        Arguments:

        - token -- the token indicating a block start
        - fh    -- an open file handler

        Return format: a list of (start, next start) tuples
        """
        fh.seek(0)
        header_start_list = []
        here = 0
        wrong = 0
        start = -1
        while fh:
            here = fh.tell()
            line = fh.readline()
            if len(line) == 0:
                header_start_list.append((start, 0))
                break
            if line.startswith(token):
                if start == -1:
                    start = here
                else:
                    header_start_list.append((start, here))
                    start = here
        return header_start_list

    def gen_block_sequences(self, token, fh):
        """Yield (head, seq, head_start, head_end) tuplefor given fh for open file.

        Arguments:

        - token -- the token indicating a block start
        - fh    -- an open file handler

        Return format:

        - head       -- a block head
        - seq        -- a block body
        - head_start -- a file pointer to block start
        - head_end   -- a file pointer to next block start or 0

        """
        header_start_list = self.get_blocks(token, fh)
        for x, y in header_start_list:
            # yeild sequence by coordinates
            yield self.get_block_sequence(x, y, fh)


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


def sc_iter_fasta_brute(file_name, inmem=False, lower=False):
    """Iter over fasta file."""

    header = None
    seq = []
    with open(file_name) as fh:
        if inmem:
            data = fh.readlines()
        else:
            data = fh
        for line in data:
            if line.startswith(">"):
                if seq or header:
                    sequence = "".join(seq)
                    if lower:
                        sequence = sequence.lower()
                    yield header, sequence
                header = line.strip()
                seq = []
                continue
            seq.append(line.strip())
        if seq or header:
            sequence = "".join(seq)
            if lower:
                sequence = sequence.lower()
            yield header, sequence


def sc_iter_tab_file(
    input_file,
    data_type,
    skip_starts_with=None,
    remove_starts_with=None,
    preprocess_function=None,
    check_function=None,
):
    """Iter over tab file, yield an object of given data_type."""

    temp_file = tempfile.NamedTemporaryFile(delete=False)
    temp_file_name = temp_file.name
    csv.field_size_limit(256 << 30)
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


def save_trs_dataset(trs_dataset, output_file, dataset_id=None):
    """Save trs dataset to file."""
    if isinstance(trs_dataset, dict):
        trs_dataset = trs_dataset.items()
        trs_dataset.sort()
        trs_dataset = [x[1] for x in trs_dataset]
    if dataset_id is None:
        with open(output_file, "w") as fh:
            for trf_obj in trs_dataset:
                data = str(trf_obj)
                fh.write(data)
    else:
        with open(output_file, "a") as fh:
            for trf_obj in trs_dataset:
                trf_obj.giid = dataset_id
                data = str(trf_obj)
                fh.write(data)


def sc_iter_fasta(file_name, lower=False, protein=False, inmem=True, skip_clean=False):
    """Iter over fasta file."""
    header = None
    seq = []
    with open(file_name) as fh:
        data = fh
        if inmem:
            data = fh.readlines()
        for line in data:
            if line.startswith(">"):
                if seq:
                    sequence = "".join(seq)
                    seq_obj = SequenceModel(lower=lower, protein=protein)
                    seq_obj.set_ncbi_sequence(header, sequence, skip_clean=skip_clean)
                    yield seq_obj
                header = line
                seq = []
                continue
            seq.append(line)
        if seq or header:
            sequence = "".join(seq)
            seq_obj = SequenceModel(lower=lower, protein=protein)
            seq_obj.set_ncbi_sequence(header, sequence)
            yield seq_obj


def sc_iter_fasta(file_name, lower=False, protein=False, inmem=True, skip_clean=False):
    """Iter over fasta file."""
    header = None
    seq = []
    with open(file_name) as fh:
        data = fh
        if inmem:
            data = fh.readlines()
        for line in data:
            if line.startswith(">"):
                if seq:
                    sequence = "".join(seq)
                    seq_obj = SequenceModel(lower=lower, protein=protein)
                    seq_obj.set_ncbi_sequence(header, sequence, skip_clean=skip_clean)
                    yield seq_obj
                header = line
                seq = []
                continue
            seq.append(line)
        if seq or header:
            sequence = "".join(seq)
            seq_obj = SequenceModel(lower=lower, protein=protein)
            seq_obj.set_ncbi_sequence(header, sequence)
            yield seq_obj


class Gff3FileIO(TabDelimitedFileIO):
    """ """

    def __init__(self, *args, **kwargs):
        """ """
        super(TabDelimitedFileIO, self).__init__(*args, **kwargs)
        self.headers = []

    def read_online(self, file_name, only_fields=None):
        """Overrided. Yield items online from data from input_file."""

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
                            except:
                                print(v)
                        _features[k] = v
                data["attributes"] = _features
                obj = Gff3Model()
                try:
                    obj.set_with_dict(data)
                except:
                    print("Can't parse features for %s" % str(data))
                    continue
                yield obj


def sc_gff3_reader(gff3_file, only_fields=None):
    """Iter over gff3 file."""
    reader = Gff3FileIO()
    for gff3_obj in reader.read_online(gff3_file, only_fields=only_fields):
        yield gff3_obj
