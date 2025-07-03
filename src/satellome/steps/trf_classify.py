#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.02.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import argparse
import os
import sys
sys.path.append("/home/akomissarov/Dropbox/workspace/PyBioSnippets/satellome/src")
import yaml

from satellome.core_functions.classification_micro import \
    scf_basic_trs_classification

def classify_trf_data(trf_prefix, output_dir, genome_size):

    base_prefix = trf_prefix
    base_file = os.path.basename(trf_prefix)

    settings = {
        "folders": {
            "data_gff3": os.path.join(output_dir, "gff3"),
            "reports": os.path.join(output_dir, "reports"),
            "trf_parsed_folder": output_dir,
            "fasta": os.path.join(output_dir, "fasta"),
        },
        "files": {},
    }

    for folder_path in settings["folders"].values():
        if not os.path.isdir(folder_path):
            os.makedirs(folder_path)
                 
    settings["files"] = {
        "trf_all_file": base_prefix + ".trf",
        
        "trf_micro_file": base_prefix + ".micro.trf",
        "trf_perfect_micro_file": base_prefix + ".pmicro.trf",
        "trf_tssr_file": base_prefix + ".tssr.trf",
        "trf_complex_file": base_prefix + ".complex.trf",
    
        "trf_1k_file": base_prefix + ".1kb.trf",
        "trf_10k_file": base_prefix + ".10kb.trf",
        "trf_1k_fasta_file": os.path.join(settings["folders"]["fasta"], base_file + ".1kb.fasta"),
        "trf_3k_file": base_prefix + ".3kb.trf",
        "trf_3k_fasta_file": os.path.join(settings["folders"]["fasta"], base_file + ".3kb.fasta"),
        "trf_10k_fasta_file": os.path.join(settings["folders"]["fasta"], base_file + ".10kb.fasta"),
        "gff_micro_file": os.path.join(settings["folders"]["data_gff3"], base_file + ".micro.gff"),
        "report_micro_file": os.path.join(settings["folders"]["reports"], base_file + ".micro.report"),
        "gff_pmicro_file": os.path.join(settings["folders"]["data_gff3"], base_file + ".pmicro.gff"),
        "report_pmicro_file":os.path.join(settings["folders"]["reports"], base_file + ".pmicro.report"),
        "gff_tssr_file": os.path.join(settings["folders"]["data_gff3"], base_file + ".tssr.gff"),
        "report_tssr_file": os.path.join(settings["folders"]["reports"], base_file + ".tssr.report"),
        "trf_fssr_file": os.path.join(settings["folders"]["data_gff3"], base_file + ".fssr.gff"),
        "report_fssr_file": os.path.join(settings["folders"]["reports"], base_file + ".fssr.report"),
        "gff_complex_file": os.path.join(settings["folders"]["data_gff3"], base_file + ".complex.gff"),
        "gff_1k_file": os.path.join(settings["folders"]["data_gff3"], base_file + ".1kb.gff"),
        "gff_3k_file": os.path.join(settings["folders"]["data_gff3"], base_file + ".3kb.gff"),
        "gff_10k_file": os.path.join(settings["folders"]["data_gff3"], base_file + ".10kb.gff"),
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
