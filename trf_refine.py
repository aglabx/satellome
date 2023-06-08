#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.02.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import argparse
from trf_model import TRModel
from trseeker.seqio.tab_file import sc_iter_tab_file
import shutil
from satelome.parsers import refine_name


def refine_names(trf_file):
    data = []
    for i, trf_obj in enumerate(sc_iter_tab_file(trf_file, TRModel)):
        refine_name(trf_obj)
        data.append(trf_obj)
    
    with open(trf_file + ".1", "w") as fw:
        for obj in data:
            fw.write(obj.get_as_string(obj.dumpable_attributes))

    shutil.move(trf_file + ".1", trf_file)
    


def main():
    args = get_args()
    trf_file = args.input
    
    print("Refining names...")
    refine_names(trf_file)


def get_args():
    parser = argparse.ArgumentParser(description='Refine TRF names')
    parser.add_argument('-i', '--input', type=str, help='TRF file')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()