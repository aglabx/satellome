#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.02.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import argparse
import shutil

from satellome.core_functions.io.tab_file import sc_iter_tab_file
from satellome.core_functions.io.trf_file import join_overlapped
from satellome.core_functions.models.trf_model import TRModel
from satellome.core_functions.tools.parsers import refine_name


def refine_names(trf_file):
    data = []
    last_head = None
    last_end = None
    for i, trf_obj in enumerate(sc_iter_tab_file(trf_file, TRModel)):

        refine_name(i, trf_obj)

        # check overlapped
        if last_head is not None:
            if trf_obj.trf_l_ind < last_end and last_head == trf_obj.trf_head:
                if join_overlapped(data[-1], trf_obj, cutoff_distance=0.1):
                    last_head = trf_obj.trf_head
                    last_end = trf_obj.trf_r_ind
                    continue

        last_head = trf_obj.trf_head
        last_end = trf_obj.trf_r_ind

        data.append(trf_obj)

    with open(trf_file + ".1", "w") as fw:
        for obj in data:
            fw.write(obj.get_as_string(obj.dumpable_attributes))

    shutil.move(trf_file + ".1", trf_file)


def main(args):
    trf_file = args.input
    refine_names(trf_file)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Refine TRF names")
    parser.add_argument("-i", "--input", type=str, help="TRF file", required=True)
    args = parser.parse_args()

    main(args)
