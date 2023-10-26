#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.02.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import argparse
import os

import yaml

from satelome.core_functions.classification_micro import \
    scf_basic_trs_classification

from satelome.core_functions.tools.processing import get_genome_size

def classify_trf_data(trf_prefix, output_dir, genome_size):

    base_prefix = trf_prefix

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
        "pid": "project",
        "work_files": {
            "ref_assembly_name_for_trf": "dataset",
            "assembly_stats": {
                "dataset": {
                    "genome_size": genome_size,
                },
            },
        },
    }

    ### PART 2. Classify according to monomer/array features

    results_file = os.path.join(output_dir, "results.yaml")

    print("Classifying TRF results...")
    settings, project = scf_basic_trs_classification(settings, project)

    print("Saving results...")
    with open(results_file, "w") as fh:
        yaml.dump(project, fh, default_flow_style=False)


def main():
    args = get_args()
    trf_prefix = args.prefix
    output_dir = args.output
    genome_size = args.genome_size

    if genome_size == 0:
        genome_size = get_genome_size(fasta_file)

    print("Refining names...")
    classify_trf_data(trf_prefix, output_dir, genome_size)


def get_args():
    parser = argparse.ArgumentParser(
        description="Classify TRF and write basic statistics"
    )
    parser.add_argument(
        "-i",
        "--prefix",
        type=str,
        help="TRF prefix (trf file without extension))",
        required=True,
    )
    parser.add_argument(
        "-o", "--output", type=str, help="Output directory", required=True
    )
    parser.add_argument(
        "-l",
        "--genome_size",
        type=int,
        help="Total length of the assembly",
        required=True,
    )
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main()
