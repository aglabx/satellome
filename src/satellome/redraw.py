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
from satellome.core_functions.io.tab_file import sc_iter_tab_file
from satellome.core_functions.models.trf_model import TRModel
from satellome.core_functions.tools.gene_intersect import add_annotation_from_gff
from satellome.core_functions.tools.reports import create_html_report

from satellome.core_functions.tools.processing import get_genome_size
from satellome.core_functions.tools.ncbi import get_taxon_name

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
        "--genome_size", help="Expected genome size [will be computed]", required=False, default=0
    )
    parser.add_argument("--taxid", help="NCBI taxid, look here https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi", required=True)
    parser.add_argument("--gff", help="Input gff file [None]", required=False, default=None)
    parser.add_argument("--rm", help="Input RM *.ori.out file [None]", required=False, default=None)
    parser.add_argument("--srr", help="SRR index for raw reads [None]", required=False, default=None)
    parser.add_argument("-c", "--cutoff", help="Cutoff for large TRs", required=False, default=1000)
    ### add minimal_scaffold_length
    parser.add_argument("-l", "--minimal_scaffold_length", help="Minimal scaffold length", required=False, default=10000)
    parser.add_argument("-e", "--drawing_enhancing", help="Drawing enhancing", required=False, default=100000)
    parser.add_argument("--large_file", help="Suffix for TR file for analysis, it can be '', 1kb, 3kb, 10kb [1kb]", required=False, default="1kb") 
    parser.add_argument("--taxon", help="Taxon name", required=False, default=None)

    args = vars(parser.parse_args())

    fasta_file = args["input"]
    output_dir = args["output"]
    project = args["project"]
    threads = args["threads"]
    trf_path = args["trf"]
    large_cutoff = int(args["cutoff"])
    genome_size = int(args["genome_size"])
    gff_file = args["gff"]
    minimal_scaffold_length = int(args["minimal_scaffold_length"])
    drawing_enhancing = int(args["drawing_enhancing"])
    taxid = args["taxid"]
    large_file_suffix = args["large_file"]
    repeatmasker_file = args["rm"]
    taxon_name = args["taxon"]

    taxon_name = taxon_name.replace(" ", "_")

    input_filename_without_extension = os.path.basename(os.path.splitext(fasta_file)[0])

    trf_prefix = os.path.join(
        output_dir,
        input_filename_without_extension
    )
    if large_file_suffix:
        trf_file = f"{trf_prefix}.{large_file_suffix}.trf"
    else:
        trf_file = f"{trf_prefix}.trf"

    if not genome_size:
        genome_size = get_genome_size(fasta_file)

    current_file_path = os.path.abspath(__file__)
    current_directory = os.path.dirname(current_file_path)

    trf_search_path = os.path.join(current_directory, "steps", "trf_search.py")
    trf_classify_path = os.path.join(current_directory, "steps", "trf_classify.py")
    trf_draw_path = os.path.join(current_directory, "steps", "trf_draw.py")

    distance_file = os.path.join(output_dir, "distances.tsv")

    html_report_file = os.path.join(output_dir, "reports", "satellome_report.html")
    if not os.path.exists(os.path.dirname(html_report_file)):
        os.makedirs(os.path.dirname(html_report_file))

    output_image_dir = os.path.join(output_dir, "images")
    if not os.path.exists(output_image_dir):
        os.makedirs(output_image_dir)
    
    settings = {
        "fasta_file": fasta_file,
        "output_dir": output_dir,
        "project": project,
        "threads": threads,
        "trf_path": trf_path,
        "genome_size": genome_size,
        "trf_prefix": trf_prefix,
        "large_cutoff": large_cutoff,
        "trf_search_path": trf_search_path,
        "trf_classify_path": trf_classify_path,
        "gff_file": gff_file,
        "trf_file": f"{trf_prefix}.trf",
        "minimal_scaffold_length": minimal_scaffold_length,
        "drawing_enhancing": drawing_enhancing,
        "taxon_name": taxon_name,
        "srr": args["srr"],
        "taxid": taxid,
        "distance_file": distance_file,
        "output_image_dir": output_image_dir,
        "large_file_suffix": large_file_suffix,
        "repeatmasker_file": repeatmasker_file,
        "html_report_file": html_report_file,
    }

    command = f"python {trf_draw_path} -f {fasta_file} -i {trf_file} -o {output_image_dir} -c {minimal_scaffold_length} -e {drawing_enhancing} -t '{taxon_name}' -d {distance_file} -s {genome_size} "
    # print(command)
    completed_process = subprocess.run(command, shell=True)
    if completed_process.returncode == 0:
        print("trf_draw.py executed successfully!")
    else:
        print(f"trf_draw.py failed with return code {completed_process.returncode}")
        sys.exit(1)

    create_html_report(output_image_dir, html_report_file)
