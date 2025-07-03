#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 26.10.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import argparse
import os
import pathlib
import sys
sys.path.append("/home/akomissarov/Dropbox/workspace/PyBioSnippets/satellome/src")

from satellome.core_functions.tools.processing import get_genome_size
from satellome.core_functions.tools.trf_tools import trf_search_by_splitting


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Parse TRF output.")
    parser.add_argument("-i", "--input", help="Input fasta file", required=True)
    parser.add_argument("-o", "--output", help="Output folder", required=True)
    parser.add_argument("-p", "--project", help="Project", required=True)
    parser.add_argument("-t", "--threads", help="Threads", required=True)
    parser.add_argument(
        "--trf", help="Path to trf [trf]", required=False, default="trf"
    )
    parser.add_argument(
        "--genome_size", help="Expected genome size", required=False, default=0
    )
    args = vars(parser.parse_args())

    fasta_file = args["input"]
    output_dir = args["output"]
    project = args["project"]
    threads = args["threads"]
    trf_path = args["trf"]
    genome_size = int(args["genome_size"])

    settings = {
        "fasta_file": fasta_file,
        "output_dir": output_dir,
        "project": project,
        "threads": threads,
        "trf_path": trf_path,
        "genome_size": genome_size,
    }

    if not output_dir.startswith("/"):
        print(f"Error: please provide the full path for output: {output_dir}")
        sys.exit(1)

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    if genome_size == 0:
        genome_size = get_genome_size(fasta_file)

    code_dir = pathlib.Path(__file__).parent.resolve()
    parser_program = os.path.join(code_dir, "trf_parse_raw.py")

    ### PART 1. Running TRF in parallel

    fasta_name = ".".join(fasta_file.split("/")[-1].split(".")[:-1])
    output_file = os.path.join(output_dir, fasta_name + ".trf")

    if os.path.isfile(output_file):
        print("TRF output file already exists. Skipping TRF.")
    else:
        # print("Running TRF...")
        output_file = trf_search_by_splitting(
            fasta_file,
            threads=threads,
            wdir=output_dir,
            project=project,
            trf_path=trf_path,
            parser_program=parser_program,
            genome_size=genome_size,
        )
