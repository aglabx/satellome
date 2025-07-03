#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 05.06.2011
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

"""
Classes:

- TRFFileIO(AbstractBlockFileIO)

"""
import os
from collections import defaultdict

from PyExp import WiseOpener

from satellome.core_functions.io.block_file import AbstractBlockFileIO
from satellome.core_functions.io.file_system import iter_filepath_folder
from satellome.core_functions.io.tab_file import sc_iter_tab_file
from satellome.core_functions.models.trf_model import TRModel
from satellome.core_functions.settings import load_settings
from satellome.core_functions.tools.parsers import refine_name
from satellome.core_functions.tools.processing import get_revcomp
from satellome.core_functions.trf_embedings import get_cosine_distance

settings = load_settings()


def join_overlapped(obj1, obj2, cutoff_distance=0.1):

    """Join overlapping sequences."""
    """ Join two overlapped objects.
    We will join two objects if they are overlapped and have the consinus distance by 5-mers vectors less than 0.1. 
    """
    # a ------
    # b    -----
    if obj1.trf_r_ind > obj2.trf_l_ind and obj1.trf_r_ind < obj2.trf_r_ind:
        vector1 = obj1.get_vector()
        vector2 = obj2.get_vector()
        dist = get_cosine_distance(vector1, vector2)
        left_part = obj2.trf_l_ind - obj1.trf_l_ind
        right_part = obj2.trf_r_ind - obj1.trf_r_ind
        middle_part = obj1.trf_r_ind - obj2.trf_l_ind

        intersect_fraction = middle_part / (left_part + right_part + middle_part)

        # print("DIST", dist, "FRAQ", intersect_fraction)

        if dist < cutoff_distance or intersect_fraction > 0.2:
            obj1.set_form_overlap(obj2)
            return True
    return False


def get_gc(sequence):
    """Count GC content."""
    length = len(sequence)
    if not length:
        return 0
    count_c = sequence.count("c") + sequence.count("C")
    count_g = sequence.count("g") + sequence.count("G")
    gc = float(count_c + count_g) / float(length)
    return float(gc)


def get_int_gc(sequence):
    """Get GC content from 0 to 100."""
    gc = get_gc(sequence)
    return int(100 * round(gc, 2))


def get_shifts_variants(sequence):
    """ """
    shifts = set()
    for i in range(len(sequence)):
        shifts.add(sequence[i:] + sequence[:i])
    return list(shifts)


def sort_dictionary_by_value(d, reverse=False):
    """Sort dictionary by value. Retrun list of (v, k) pairs."""
    result = [(v, k) for k, v in list(d.items())]
    result.sort(reverse=reverse, key=lambda x: x[0])
    return result


def remove_consensus_redundancy(trf_objs):
    """Take a minimal sequence from lexicographically sorted rotations of sequence and its reverse complement
    Example: ACT, underlined - reverse complement seqeunces
    ACT, AGT, CTA, GTA, TAC, TAG
    Find all possible multimers, e.g. replace GTAGTAGTA consensus sequence with ACT
    Return:
    1) list of sorted TRs
    2) list of (df, consensus) pairs
    """
    # sort by length
    consensuses = [x.trf_consensus for x in trf_objs]
    consensuses = list(set(consensuses))
    consensuses.sort(key=lambda x: len(x))
    length2consensuses = {}
    # print "Group monomers by length and GC"
    for i, monomer in enumerate(consensuses):
        n = len(monomer)
        length2consensuses.setdefault(n, {})
        gc = get_int_gc(monomer)
        length2consensuses[n].setdefault(gc, [])
        length2consensuses[n][gc].append(i)
    # print "Iterate over consensuses"
    N = len(consensuses)
    result_rules = {}
    for i, monomer in enumerate(consensuses):
        if not monomer:
            continue
        # print i, N, "\r",
        if monomer in result_rules:
            continue
        gc = get_int_gc(monomer)
        base = len(monomer)
        n = base
        # maximal consensus length from TRF is 2000 bp
        variants = set(
            get_shifts_variants(monomer) + get_shifts_variants(get_revcomp(monomer))
        )
        if not variants:
            raise Exception("Wrong monomer sequence for '%s'" % monomer)
        lex_consensus = min(variants)
        while n <= 2020:
            if n in length2consensuses and gc in length2consensuses[n]:
                for k in length2consensuses[n][gc]:
                    monomer_b = consensuses[k]
                    if monomer_b in result_rules:
                        continue
                    s = n // base
                    v = set()
                    for p in range(s):
                        v.add(monomer_b[p * base : (p + 1) * base])
                    if len(v) > 1:
                        continue
                    item = v.pop()
                    if item in variants:
                        # if consensuses[k] != lex_consensus:
                        # print
                        # print i, base, consensuses[k], "->", lex_consensus
                        result_rules[consensuses[k]] = lex_consensus

            n += base
    # print
    # print "Fix momomers"
    variants2df = defaultdict(int)
    for i, trf_obj in enumerate(trf_objs):
        if not trf_obj.trf_consensus:
            trf_objs[i] = None
            continue
        if trf_obj.trf_consensus in result_rules:
            variants2df[result_rules[trf_obj.trf_consensus]] += 1
        else:
            message = (
                "Error key with length",
                len(trf_obj.trf_consensus),
                trf_obj.trf_consensus,
            )
            # print str(trf_obj)
            # print message
            raise Exception(message)
        trf_obj.trf_consensus = result_rules[trf_obj.trf_consensus]
    # print "Sort families by df..."
    variants2df = sort_dictionary_by_value(variants2df, reverse=True)
    trf_objs = [x for x in trf_objs if x is not None]
    return trf_objs, variants2df


class TRFFileIO(AbstractBlockFileIO):
    """Working with raw ouput from TRF, where each block starts with '>' token.

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
        """Overrided. Hardcoded start token."""
        token = "Sequence:"
        super(TRFFileIO, self).__init__(token)

    def iter_parse(self, trf_file, filter=True):
        """Iterate over raw trf data and yield TRFObjs."""
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

    def parse_to_file(
        self, file_path, output_path, trf_id=0, project=None, verbose=True
    ):
        """Parse trf file in tab delimited file."""
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
                    refine_name(trf_id, trf_obj)

                    fw.write(str(trf_obj))

                    trf_id += 1
        return trf_id

    def refine_old_to_file(
        self, file_path, output_path, trf_id=0, project=None, verbose=True
    ):
        """Parse trf file in tab delimited file."""
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
                    refine_name(trf_id, trf_obj)

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
                if (
                    obj1.trf_l_ind == obj2.trf_l_ind
                    and obj1.trf_r_ind == obj2.trf_r_ind
                ):
                    # Check period
                    if obj1.trf_pmatch >= obj2.trf_pmatch:
                        obj_set[b] = None
                    else:
                        obj_set[a] = None
                    continue
                # a ------ ------  -------
                # b ---       ---    ---
                if (
                    obj1.trf_l_ind <= obj2.trf_l_ind
                    and obj1.trf_r_ind >= obj2.trf_r_ind
                ):
                    obj_set[b] = None
                    continue
                # a ---       ---    ---
                # b ------ ------  -------
                if (
                    obj2.trf_l_ind <= obj1.trf_l_ind
                    and obj2.trf_r_ind >= obj1.trf_r_ind
                ):
                    obj_set[a] = None
                    continue
                # a ------
                # b    -----
                if obj1.trf_r_ind > obj2.trf_l_ind and obj1.trf_r_ind < obj2.trf_r_ind:
                    # print(obj1.as_dict(), obj2.as_dict())
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
                    if (
                        obj1.trf_r_ind > obj2.trf_l_ind
                        and obj1.trf_r_ind < obj2.trf_r_ind
                    ):

                        overlap = float(abs(obj1.trf_r_ind - obj2.trf_l_ind))
                        min_length = min(obj1.trf_array_length, obj2.trf_array_length)
                        overlap_proc_diff = overlap * 1.0 / min_length
                        gc_dif = abs(obj1.trf_array_gc - obj2.trf_array_gc)

                        if (
                            overlap_proc_diff
                            >= settings["trf_settings"]["overlapping_cutoff_proc"]
                            and gc_dif
                            <= settings["trf_settings"]["overlapping_gc_diff"]
                        ):
                            is_overlapping = True
                            if self._join_overlapped(obj1, obj2):
                                obj2 = None
                            # print("overlap: ", overlap, "min_length:", min_length, "overlap_proc_diff:", overlap_proc_diff, "gc_dif:", gc_dif)
                            # print("JOINED")
                        continue
                    # a ------
                    # b ------
                    if (
                        obj1.trf_l_ind == obj2.trf_l_ind
                        and obj1.trf_r_ind == obj2.trf_r_ind
                    ):
                        # Check period
                        if obj1.trf_pmatch >= obj2.trf_pmatch:
                            obj_set[b] = None
                            continue
                        else:
                            obj_set[a] = None
                            continue
                    # a ------ ------  -------
                    # b ---       ---     ---
                    if (
                        obj1.trf_l_ind <= obj2.trf_l_ind
                        and obj1.trf_r_ind >= obj2.trf_r_ind
                    ):
                        obj_set[b] = None
                        continue
                    # a ---       ---            ---
                    # b ------ ------  -------
                    if (
                        obj2.trf_l_ind <= obj1.trf_l_ind
                        and obj2.trf_r_ind >= obj1.trf_r_ind
                    ):
                        obj_set[a] = None
                        continue

            obj_set = [a for a in obj_set if not a is None]

        return obj_set

    def _join_overlapped(self, obj1, obj2, cutoff_distance=0.1):
        return join_overlapped(obj1, obj2, cutoff_distance=cutoff_distance)


def sc_parse_raw_trf_folder(trf_raw_folder, output_trf_file, project=None):
    """Parse raw TRF output in given folder to output_trf_file."""
    reader = TRFFileIO()
    trf_id = 1
    if os.path.isfile(output_trf_file):
        os.remove(output_trf_file)
    for file_path in iter_filepath_folder(trf_raw_folder, mask="dat"):
        if not file_path.endswith(".dat"):
            continue
        print("Start parse file %s..." % file_path)
        trf_id = reader.parse_to_file(
            file_path, output_trf_file, trf_id=trf_id, project=project
        )


def sc_trf_to_fasta(trf_file, fasta_file):
    """Convert TRF file to fasta file."""
    with open(fasta_file, "w") as fw:
        for trf_obj in sc_iter_tab_file(trf_file, TRModel):
            fw.write(trf_obj.fasta)
