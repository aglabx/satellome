#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 26.10.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import argparse
import subprocess
import sys
import os

from satelome.core_functions.tools.processing import get_genome_size

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


    input_filename_without_extension = os.path.splitext(fasta_file)[0]

    trf_prefix = os.path.join(
        output_dir,
        input_filename_without_extension
    ) 

    settings = {
        "fasta_file": fasta_file,
        "output_dir": output_dir,
        "project": project,
        "threads": threads,
        "trf_path": trf_path,
        "genome_size": genome_size,
        "trf_prefix": trf_prefix,
    }

    if genome_size == 0:
        genome_size = get_genome_size(fasta_file)

    current_file_path = os.path.abspath(__file__)
    current_directory = os.path.dirname(current_file_path)

    trf_search_path = os.path.join(current_directory, "trf_search.py")

    command = f"time python {trf_search_path} -i {fasta_file} \
                                   -o {output_dir} \
                                   -p {project} \
                                   -t {threads} \
                                   --trf {trf_path} \
                                   --genome_size {genome_size}"
    print(command)
    completed_process = subprocess.run(command, shell=True)

    if completed_process.returncode == 0:
        print("trf_search.py executed successfully!")
    else:
        print(f"trf_search.py failed with return code {completed_process.returncode}")
        sys.exit(1)

    trf_search_path = os.path.join(current_directory, "trf_classify.py")

    command = f"time python {trf_search_path} -i {trf_prefix} -o {output_dir} -l {genome_size}"
    print(command)
    completed_process = subprocess.run(command, shell=True)
    if completed_process.returncode == 0:
        print("trf_classify.py executed successfully!")
    else:
        print(f"trf_classify.py failed with return code {completed_process.returncode}")
        sys.exit(1)
