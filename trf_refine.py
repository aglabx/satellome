#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.02.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import argparse
import shutil
import logging
from trf_model import TRModel
from seqio import sc_iter_tab_file


def refine_names(trf_file):
    data = []
    for i, obj in enumerate(sc_iter_tab_file(trf_file, TRModel)):
        name = obj.trf_head.split()
        if len(name):
            name = name[0]
        else:
            name = name
        obj.trf_id = f"{name}_{obj.trf_l_ind}_{obj.trf_r_ind}"
        obj.id = f"AGT{(i+1) * 100:013d}"
        obj.trf_consensus = obj.trf_consensus.upper()
        obj.trf_array = obj.trf_array.upper()
        data.append(obj)
    
    with open(trf_file + ".1", "w") as fw:
        for obj in data:
            fw.write(obj.get_as_string(obj.dumpable_attributes))

    shutil.move(trf_file + ".1", trf_file)
    


def main():
    args = get_args()
    trf_file = args.input
    
    logging.info("Refining names...")
    refine_names(trf_file)


def get_args():
    parser = argparse.ArgumentParser(description='Refine TRF names')
    parser.add_argument('-i', '--input', type=str, help='TRF file')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()