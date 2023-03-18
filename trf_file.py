#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 05.06.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 

"""
Classes:

- TRFFileIO(AbstractBlockFileIO)

"""
import re, os
from trf_model import TRModel
from trseeker.seqio.block_file import AbstractBlockFileIO
from trseeker.tools.sequence_tools import get_gc, remove_consensus_redundancy
from PyExp import sc_iter_filepath_folder, WiseOpener
from trseeker.settings import load_settings
from trseeker.seqio.tab_file import sc_iter_tab_file
from trf_embedings import get_cosine_distance

settings = load_settings()


class TRFFileIO(AbstractBlockFileIO):
    """  Working with raw ouput from TRF, where each block starts with '>' token.
    
    Public parameters:
    
    - self.use_mongodb -- Bool

    Public methods:
    
    - iter_parse(self, trf_file, filter=True)
    - parse_to_file(self, file_path, output_path, trf_id=0) -> trf_id
    
    Private methods:
    
    - _gen_data_line(self, data)
    - _filter_obj_set(self, obj_set)
    - _join_overlapped(self, obj1, obj2)
    
    Inherited public properties:
    
    - data  - iterable data, each item is tuple (head, body)
    - N     - a number of items in data
    
    Inherited public methods:
    
    - [OR] __init__(self)
    - read_from_file(self, input_file)
    - read_online(self, input_file) ~> item
    - get_block_sequence(self, head_start, next_head, fh)
    - get_blocks(self, token, fh)
    - gen_block_sequences(self, token, fh)
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

    def __init__(self):
        """ Overrided. Hardcoded start token."""
        token = "Sequence:"
        super(TRFFileIO, self).__init__(token)

    def iter_parse(self, trf_file, filter=True):
        """ Iterate over raw trf data and yield TRFObjs.
        """
        trf_id = 1
        for ii, (head, body, start, next) in enumerate(self.read_online(trf_file)):
            head = head.replace("\t", " ")
            obj_set = []
            # print(" processing:", head)
            n = body.count("\n")
            for i, line in enumerate(self._gen_data_line(body)):
                # print("   %s/%s" % (i,n), "\r", end=" ")
                trf_obj = TRModel()
                trf_obj.set_raw_trf(head, None, line)
                obj_set.append(trf_obj)
            # print(" filtering...")
            if filter:
                # Filter object set
                trf_obj_set = self._filter_obj_set(obj_set)
                obj_set = [x for x in trf_obj_set if x]
            ### set trf_id
            for trf_obj in obj_set:
                trf_obj.trf_id = trf_id
                trf_id += 1
            # print(" fixing monomers...")
            obj_set, variants2df = remove_consensus_redundancy(obj_set)
            yield obj_set

    def parse_to_file(self, file_path, output_path, trf_id=0, project=None, verbose=True):
        """ Parse trf file in tab delimited file."""
        if trf_id == 0:
            mode = "w"
        else:
            mode = "a"
        
        with WiseOpener(output_path, mode) as fw:
            for trf_obj_set in self.iter_parse(file_path):
                for trf_obj in trf_obj_set:
                    trf_obj.trf_id = trf_id

                    if project:
                        trf_obj.set_project_data(project)

                    fw.write(str(trf_obj))

                    trf_id += 1
        return trf_id

    def _gen_data_line(self, data):
        for line in data.split("\n"):
            line = line.strip()
            if line.startswith("Sequence"):
                continue
            if line.startswith("Parameters"):
                continue
            if not line:
                continue
            yield line

    def _join_overlapped(self, obj1, obj2, cutoff_distance=0.1):
        """ Join two overlapped objects.
        We will join two objects if they are overlapped and have the consinus distance by 5-mers vectors less than 0.1. 
        """
        # a ------ 
        # b    -----
        if obj1.trf_r_ind > obj2.trf_l_ind and obj1.trf_r_ind < obj2.trf_r_ind:
            vector1 = obj1.get_vector()
            vector2 = obj2.get_vector()
            dist = get_cosine_distance(vector1, vector2)
            if dist < cutoff_distance:
                obj1.set_form_overlap(obj2)
                return True
        return False
        
    def _filter_obj_set(self, obj_set):
        # NB: I removed the overlaping part due to suspicious results.
        # Complex filter
        is_overlapping = False
        n = len(obj_set)

        obj_set.sort(key=lambda x: (x.trf_l_ind, x.trf_r_ind))
        for a in range(0, n):
            # print("\t%s\%s" % (a,n), "\r", end=" ")

            obj1 = obj_set[a]
            if not obj1:
                continue
            for b in range(a + 1, n):
                obj2 = obj_set[b]
                if not obj2:
                    continue
                # a ------ 
                # b ------
                if obj1.trf_l_ind == obj2.trf_l_ind and obj1.trf_r_ind == obj2.trf_r_ind:
                    # Check period
                    if obj1.trf_pmatch >= obj2.trf_pmatch:
                        obj_set[b] = None
                    else:
                        obj_set[a] = None
                    continue
                # a ------ ------  ------- 
                # b ---       ---    ---
                if obj1.trf_l_ind <= obj2.trf_l_ind and obj1.trf_r_ind >= obj2.trf_r_ind:
                    obj_set[b] = None
                    continue
                # a ---       ---    ---
                # b ------ ------  ------- 
                if obj2.trf_l_ind <= obj1.trf_l_ind and obj2.trf_r_ind >= obj1.trf_r_ind:
                    obj_set[a] = None
                    continue
                # a ------ 
                # b    -----
                if obj1.trf_r_ind > obj2.trf_l_ind and obj1.trf_r_ind < obj2.trf_r_ind:
                    if self._join_overlapped(obj1, obj2, cutoff_distance=0.1):
                        obj_set[b] = None
                    continue
                # a ------ 
                # b                -----
                if obj1.trf_r_ind < obj2.trf_l_ind:
                    break
                # a               ------ 
                # b -----
                if obj2.trf_r_ind < obj1.trf_l_ind:
                    break
        obj_set = [a for a in obj_set if not a is None]
        n = len(obj_set)

        while is_overlapping:
            is_overlapping = False

            for a in range(0, n):
                obj1 = obj_set[a]
                if not obj1:
                    continue
                for b in range(a + 1, n):
                    obj2 = obj_set[b]
                    if not obj2:
                        continue
                    # a ------ 
                    # b               -----
                    if obj1.trf_r_ind < obj2.trf_l_ind:
                        break
                    # a              ------ 
                    # b -----
                    if obj2.trf_r_ind < obj1.trf_l_ind:
                        break
                    # a ------ 
                    # b    -----
                    if obj1.trf_r_ind > obj2.trf_l_ind and obj1.trf_r_ind < obj2.trf_r_ind:

                        overlap = float(abs(obj1.trf_r_ind - obj2.trf_l_ind))
                        min_length = min(obj1.trf_array_length, obj2.trf_array_length)
                        overlap_proc_diff = overlap * 1.0 / min_length
                        gc_dif = abs(obj1.trf_array_gc - obj2.trf_array_gc)

                        if overlap_proc_diff >= settings["trf_settings"]["overlapping_cutoff_proc"] \
                                    and gc_dif <= settings["trf_settings"]["overlapping_gc_diff"]:
                            is_overlapping = True
                            obj1 = self._join_overlapped(obj1, obj2)
                            obj2 = None
                            # print("overlap: ", overlap, "min_length:", min_length, "overlap_proc_diff:", overlap_proc_diff, "gc_dif:", gc_dif)
                            # print("JOINED")
                        continue
                    # a ------ 
                    # b ------
                    if obj1.trf_l_ind == obj2.trf_l_ind and obj1.trf_r_ind == obj2.trf_r_ind:
                        # Check period
                        if obj1.trf_pmatch >= obj2.trf_pmatch:
                            obj_set[b] = None
                            continue
                        else:
                            obj_set[a] = None
                            continue
                    # a ------ ------  ------- 
                    # b ---       ---     ---
                    if obj1.trf_l_ind <= obj2.trf_l_ind and obj1.trf_r_ind >= obj2.trf_r_ind:
                        obj_set[b] = None
                        continue
                    # a ---       ---            ---
                    # b ------ ------  ------- 
                    if obj2.trf_l_ind <= obj1.trf_l_ind and obj2.trf_r_ind >= obj1.trf_r_ind:
                        obj_set[a] = None
                        continue

            obj_set = [a for a in obj_set if not a is None]
        
        return obj_set

    def _join_overlapped(self, obj1, obj2):
        ''' Join overlapping sequences.'''
        obj1.trf_pmatch = int(obj1.trf_pmatch)
        obj2.trf_pmatch = int(obj2.trf_pmatch)
        obj1.trf_array_length = int(obj1.trf_array_length)
        obj2.trf_array_length = int(obj2.trf_array_length)

        obj1.trf_array = obj1.trf_array + obj2.trf_array[obj1.trf_r_ind - obj2.trf_l_ind:]

        if obj1.trf_array_length < obj2.trf_array_length:

            obj1.trf_consensus = obj2.trf_consensus
            obj1.trf_consensus_gc = obj2.trf_consensus_gc
            obj1.trf_period = obj2.trf_period
            obj1.trf_l_cons = obj2.trf_l_cons

            obj1.trf_n_a = obj2.trf_n_a
            obj1.trf_n_t = obj2.trf_n_t
            obj1.trf_n_c = obj2.trf_n_c
            obj1.trf_n_g = obj2.trf_n_g


        obj1.trf_pmatch = int((obj1.trf_pmatch * obj1.trf_array_length + obj2.trf_pmatch * obj2.trf_array_length) / (obj1.trf_array_length + obj2.trf_array_length))

        obj1.trf_l_ind = min(obj1.trf_l_ind, obj2.trf_l_ind)
        obj1.trf_r_ind = max(obj1.trf_r_ind, obj2.trf_r_ind)
        obj1.trf_array_length = obj1.trf_r_ind - obj1.trf_l_ind
        obj1.trf_n_copy = obj1.trf_array_length / obj1.trf_period

        obj1.trf_array_gc = get_gc(obj1.trf_array)

        obj1.trf_joined = 1

        return obj1

def sc_parse_raw_trf_folder(trf_raw_folder, output_trf_file, project=None):
    """ Parse raw TRF output in given folder to output_trf_file."""
    reader = TRFFileIO()
    trf_id = 1
    if os.path.isfile(output_trf_file):
        os.remove(output_trf_file)
    for file_path in sc_iter_filepath_folder(trf_raw_folder, mask="dat"):
        print("Start parse file %s..." % file_path)
        trf_id = reader.parse_to_file(file_path, output_trf_file, trf_id=trf_id, project=project)

def sc_trf_to_fasta(trf_file, fasta_file):
    """ Convert TRF file to fasta file.
    """
    with open(fasta_file, "w") as fw:
        for trf_obj in sc_iter_tab_file(trf_file, TRModel):
            fw.write(trf_obj.fasta)


