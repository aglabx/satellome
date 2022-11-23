#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 10.10.2013
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

from trseeker.models.trf_model import TRModel
from trseeker.seqio.tab_file import sc_iter_tab_file
from trseeker.seqio.tr_file import save_trs_dataset
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Get large TRs.")
    parser.add_argument("-i", "--input", help="Input file", required=True)
    parser.add_argument("-o", "--output", help="Output file", required=True)
    args = vars(parser.parse_args())

    input_file = args["input"]
    output_file = args["output"]

    i = 0
    with open(output_file, "w") as fh:
        for i, trf_obj in enumerate(sc_iter_tab_file(input_file, TRModel)):
            s = ">%s\t%s\n%s\n" % (
                trf_obj.trf_id,
                trf_obj.trf_family,
                trf_obj.trf_consensus * 3,
            )
            print(i, s)
            fh.write(s)
        print()
