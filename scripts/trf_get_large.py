#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 10.10.2013
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import argparse

from satellome.core_functions.io.tab_file import sc_iter_tab_file
from satellome.core_functions.io.tr_file import save_trs_dataset
from satellome.core_functions.models.trf_model import TRModel

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Get large TRs.")
    parser.add_argument("-i", "--input", help="Input file", required=True)
    parser.add_argument("-o", "--output", help="Output file", required=True)
    parser.add_argument("-c", "--cutoff", help="Cutoff", required=True)
    args = vars(parser.parse_args())

    input_file = args["input"]
    output_file = args["output"]
    cutoff = int(args["cutoff"])

    trf_objs = []
    i = 0
    for j, trf_obj in enumerate(sc_iter_tab_file(input_file, TRModel)):
        if trf_obj.trf_array_length > cutoff:
            i += 1
            trf_objs.append(trf_obj)
            print(
                i,
                j,
                trf_obj.trf_array_length,
                trf_obj.trf_consensus,
                trf_obj.trf_head,
                trf_obj.trf_l_ind,
            )
    print
    save_trs_dataset(trf_objs, output_file)
