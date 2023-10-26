#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 05.06.2011
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

"""
Classes:
    
- AbstractBlockFileIO(AbstractFileIO)
    
"""
import gzip

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
