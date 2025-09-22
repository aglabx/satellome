#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# @created: 10.03.2019
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import argparse
import sys
import os

# Add parent directories to path for module imports
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(os.path.dirname(current_dir))
sys.path.insert(0, parent_dir)

from satellome.core_functions.io.trf_file import TRFFileIO

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Parse TRF output.")
    parser.add_argument("-i", "--input", help="Input file", required=True)
    parser.add_argument("-o", "--output", help="Output file", required=True)
    parser.add_argument("-p", "--project", help="Project", required=True)
    args = vars(parser.parse_args())

    input_file = args["input"]
    output_file = args["output"]
    project = args["project"]

    reader = TRFFileIO()
    trf_id = 1
    trf_id = reader.parse_to_file(
        input_file, output_file, trf_id=trf_id, project=project
    )
