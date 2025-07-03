#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 15.02.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import argparse
import os
import statistics

import aindex

from satellome.core_functions.io.tab_file import sc_iter_tab_file
from satellome.core_functions.models.trf_model import TRModel


def main(args):

    fastq_files = args.input
    output_prefix = args.output
    freq_cutoff = args.cutoff
    aindex_path = args.path
    threads = args.threads
    trf_file = args.trf
    result_file = args.results

    pf_file = os.path.join(output_prefix + ".23.pf")

    if not os.path.isfile(pf_file):
        command = f"python {aindex_path} -i {fastq_files} -t fastq -P {threads} -o {output_prefix} --lu {freq_cutoff} -s True --onlyindex True"
        print(command)
        os.system(command)

    print("Load aindex")
    aindex_settings = {
        "index_prefix": os.path.join(output_prefix + ".23"),
        "aindex_prefix": os.path.join(output_prefix + ".23"),
    }
    kmer2tf = aindex.load_aindex(aindex_settings, skip_aindex=True, skip_reads=True)

    ### trid, min freq, max freq, median freq, freqs
    N = 0
    for trf_obj in sc_iter_tab_file(trf_file, TRModel):
        N += 1

    results = []
    for ii, trf_obj in enumerate(sc_iter_tab_file(trf_file, TRModel)):
        if ii % 1000 == 0:
            print(f"{ii}/{N}", end=" ")
        if trf_obj.trf_array_length < 23:
            continue
        freqs = [
            kmer2tf[trf_obj.trf_array[i : i + 23]]
            for i in range(len(trf_obj.trf_array) - 23 + 1)
        ]
        results.append(
            (trf_obj.trf_id, min(freqs), max(freqs), statistics.median(freqs), freqs)
        )

    print("Save results")
    with open(result_file, "w") as fw:
        for d in results:
            d[-1] = ",".join(map(str, d[-1]))
            fw.write("\t".join(map(str, d)))
            fw.write("\n")

    return result_file


def get_args():
    parser = argparse.ArgumentParser(description="Compute aindex for raw reads")
    parser.add_argument(
        "-i", "--input", type=str, help="Fastq files separated by comma", required=True
    )
    parser.add_argument("-o", "--output", type=str, help="Output prefix", required=True)
    parser.add_argument(
        "-c",
        "--cutoff",
        type=int,
        default=1000,
        help="Minimal frequency of kmer for filtering",
        required=False,
    )
    parser.add_argument(
        "--path",
        type=str,
        help="Theh full path to compute_aindex.py",
        default="./compute_aindex.py",
        required=False,
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=1, help="Number of threads", required=True
    )
    parser.add_argument(
        "-j", "--trf", type=str, default=None, help="TRF file", required=True
    )
    parser.add_argument(
        "-r", "--results", type=str, default=None, help="Results file", required=True
    )
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    # parser = argparse.ArgumentParser(description='Compute index.')
    # parser.add_argument('-i', help='Fasta, comma separated fastqs or reads', required=True)
    # parser.add_argument('-j', help='JF2 file if exists (None)', required=False, default=None)
    # parser.add_argument('-t', help='Reads type reads, fasta, fastq, se (reads)', required=False, default="fastq")
    # parser.add_argument('-o', help='Output prefix', required=True)
    # parser.add_argument('-H', help='Build header file for fasta', required=False, default=True)
    # parser.add_argument('--lu', help='-L for jellyfish', required=False, default=0)
    # parser.add_argument('-s', help='Sort dat file (None)', required=False, default=None)
    # parser.add_argument('--interactive', help='Interactive (False)', required=False, default=None)
    # parser.add_argument('-P', help='Threads (12)', required=False, default=12)
    # parser.add_argument('-M', help='JF2 memory in Gb (5)', required=False, default=5)
    # parser.add_argument('--onlyindex', help='Compute only index', required=False, default=False)
    # parser.add_argument('--unzip', help='Unzip files', required=False, default=False)
    # parser.add_argument('--kmers', help='Make kmers file', required=False, default=False)

    # args = vars(parser.parse_args())

    args = get_args()
    main(args)
