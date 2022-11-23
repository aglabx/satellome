#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 10.03.2019
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import argparse
import os
import sys

import yaml
from trseeker.seqio.fasta_file import sc_iter_fasta_brute
from trseeker.tools.trf_tools import trf_search_by_splitting

from core_functions.classification_micro import scf_basic_trs_classification

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Parse TRF output.")
    parser.add_argument("-i", "--input", help="Input fasta file", required=True)
    parser.add_argument("-o", "--output", help="Output folder", required=True)
    parser.add_argument("-p", "--project", help="Project", required=True)
    parser.add_argument("-r", "--results", help="Results yaml file", required=True)
    parser.add_argument("-t", "--threads", help="Threads", required=True)
    parser.add_argument(
        "--genome_size", help="Expected genome size", required=False, default=0
    )
    args = vars(parser.parse_args())

    fasta_file = args["input"]
    output_dir = args["output"]
    project = args["project"]
    threads = args["threads"]
    results_file = args["results"]
    genome_size = int(args["genome_size"])
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    if not output_dir.startswith("/"):
        print(f"Error: please provide the full path for output: {output_dir}")
        sys.exit(1)

    if genome_size == 0:
        print("Computing genome size...", end=" ")
        for header, seq in sc_iter_fasta_brute(fasta_file):
            genome_size += len(seq)
        print(f"{genome_size} bp.")

    output_file = trf_search_by_splitting(
        fasta_file, threads=threads, wdir=output_dir, project=project
    )

    base_prefix = os.path.splitext(output_file)[0]

    settings = {
        "folders": {
            "data_gff3": os.path.join(output_dir, "data_gff3"),
            "reports": os.path.join(output_dir, "reports"),
            "trf_parsed_folder": output_dir,
            "mathematica": os.path.join(output_dir, "mathematica"),
        },
        "files": {
            "trf_all_file": base_prefix + ".trf",
            "trf_micro_file": base_prefix + ".micro.trf",
            "trf_micro_kmers_file": base_prefix + ".micro.kmers",
            "trf_pmicro_kmers_file": base_prefix + ".pmicro.kmers",
            "gff_micro_file": base_prefix + ".micro.gff",
            "report_micro_file": base_prefix + ".micro.report",
            "trf_perfect_micro_file": base_prefix + ".pmicro.trf",
            "gff_pmicro_file": base_prefix + ".pmicro.gff",
            "report_pmicro_file": base_prefix + ".pmicro.report",
            "gff_tssr_file": base_prefix + ".tssr.gff",
            "trf_tssr_file": base_prefix + ".tssr.trf",
            "report_tssr_file": base_prefix + ".tssr.report",
            "trf_fssr_file": base_prefix + ".fssr.gff",
            "report_fssr_file": base_prefix + ".fssr.report",
            "trf_complex_file": base_prefix + ".complex.%s.trf",
            "gff_complex_file": base_prefix + ".complex.%s.gff",
            "clouds_complex_file": base_prefix + ".complex.%s.clouds",
            "clouds_png_complex_file": base_prefix + ".complex.%s.clouds.png",
            "trf_1k_file": base_prefix + ".1kb.trf",
            "gff_1k_file": base_prefix + ".1kb.gff",
            "trf_1k_fasta_file": base_prefix + ".1kb.fasta",
            "trf_3k_file": base_prefix + ".3kb.trf",
            "gff_3k_file": base_prefix + ".3kb.gff",
            "trf_3k_fasta_file": base_prefix + ".3kb.fasta",
            "trf_10k_file": base_prefix + ".10kb.trf",
            "gff_10k_file": base_prefix + ".10kb.gff",
            "trf_10k_fasta_file": base_prefix + ".10kb.fasta",
        },
    }

    project = {
        "pid": project,
        "work_files": {
            "ref_assembly_name_for_trf": "dataset",
            "assembly_stats": {
                "dataset": {
                    "total_length": genome_size,
                },
            },
        },
    }

    settings, project = scf_basic_trs_classification(settings, project)

    with open(results_file, "w") as fh:
        yaml.dump(project, fh, default_flow_style=False)
