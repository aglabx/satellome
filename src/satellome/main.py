#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 26.10.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))
import argparse
import subprocess

import os
from satellome.core_functions.io.tab_file import sc_iter_tab_file
from satellome.core_functions.models.trf_model import TRModel
from satellome.core_functions.tools.gene_intersect import add_annotation_from_gff
from satellome.core_functions.tools.reports import create_html_report

from satellome.core_functions.tools.processing import get_genome_size_with_progress
from satellome.core_functions.tools.ncbi import get_taxon_name

def print_logo():
    '''https://patorjk.com/software/taag/#p=display&f=Ghost&t=AGLABX%0Asatellome
    '''
    print('''
   ('-.                             ('-.    .-. .-') ) (`-.                      
  ( OO ).-.                        ( OO ).-.\  ( OO ) ( OO ).                    
  / . --. /  ,----.     ,--.       / . --. / ;-----.\(_/.  \_)-.                 
  | \-.  \  '  .-./-')  |  |.-')   | \-.  \  | .-.  | \  `.'  /                  
.-'-'  |  | |  |_( O- ) |  | OO ).-'-'  |  | | '-' /_) \     /\                  
 \| |_.'  | |  | .--, \ |  |`-' | \| |_.'  | | .-. `.   \   \ |                  
  |  .-.  |(|  | '. (_/(|  '---.'  |  .-.  | | |  \  | .'    \_)                 
  |  | |  | |  '--'  |  |      |   |  | |  | | '--'  //  .'.  \                  
  `--' `--'  `------'   `------'   `--' `--' `------''--'   '--'                 
  .-')     ('-.     .-') _     ('-.                        _   .-')       ('-.   
 ( OO ).  ( OO ).-.(  OO) )  _(  OO)                      ( '.( OO )_   _(  OO)  
(_)---\_) / . --. //     '._(,------.,--.      .-'),-----. ,--.   ,--.)(,------. 
/    _ |  | \-.  \ |'--...__)|  .---'|  |.-') ( OO'  .-.  '|   `.'   |  |  .---' 
\  :` `..-'-'  |  |'--.  .--'|  |    |  | OO )/   |  | |  ||         |  |  |     
 '..`''.)\| |_.'  |   |  |  (|  '--. |  |`-' |\_) |  |\|  ||  |'.'|  | (|  '--.  
.-._)   \ |  .-.  |   |  |   |  .--'(|  '---.'  \ |  | |  ||  |   |  |  |  .--'  
\       / |  | |  |   |  |   |  `---.|      |    `'  '-'  '|  |   |  |  |  `---. 
 `-----'  `--' `--'   `--'   `------'`------'      `-----' `--'   `--'  `------' 
''')
    


def main():
    parser = argparse.ArgumentParser(description="Parse TRF output.")
    parser.add_argument("-i", "--input", help="Input fasta file", required=True)
    parser.add_argument("-o", "--output", help="Output folder", required=True)
    parser.add_argument("-p", "--project", help="Project", required=True)
    parser.add_argument("-t", "--threads", help="Threads", required=True)
    parser.add_argument(
        "--trf", help="Path to trf [trf]", required=False, default="trf"
    )
    parser.add_argument(
        "--genome_size", help="Expected genome size [will be computed from fasta]", required=False, default=0
    )
    parser.add_argument("--taxid", help="NCBI taxid, look here https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi [None]", required=False, default=None)
    parser.add_argument("--gff", help="Input gff file [None]", required=False, default=None)
    parser.add_argument("--rm", help="Input RM *.ori.out file [None]", required=False, default=None)
    parser.add_argument("--srr", help="SRR index for raw reads [None]", required=False, default=None)
    parser.add_argument("-c", "--cutoff", help="Cutoff for large TRs [1000]", required=False, default=1000)
    ### add minimal_scaffold_length
    parser.add_argument("-l", "--minimal_scaffold_length", help="Minimal scaffold length [10000]", required=False, default=10000)
    parser.add_argument("-e", "--drawing_enhancing", help="Drawing enhancing [100000]", required=False, default=100000)
    parser.add_argument("--large_file", help="Suffix for TR file for analysis, it can be '', 1kb, 3kb, 10kb [1kb]", required=False, default="1kb") 
    parser.add_argument("--taxon", help="Taxon name [Unknown]", required=False, default=None)
    parser.add_argument("--force", help="Force rerun all steps even if output files exist", action='store_true', default=False)
    parser.add_argument("--use_kmer_filter", help="Use k-mer profiling to filter repeat-poor regions", action='store_true', default=False)
    parser.add_argument("--kmer_threshold", help="Unique k-mer threshold for repeat detection [90000]", required=False, default=90000, type=int)
    parser.add_argument("--kmer_bed", help="Pre-computed k-mer profile BED file from varprofiler", required=False, default=None)

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
    force_rerun = args["force"]
    use_kmer_filter = args["use_kmer_filter"]
    kmer_threshold = args["kmer_threshold"]
    kmer_bed_file = args["kmer_bed"]

    print_logo()

    print(f"Starting Satellome analysis...")
    print(f"Project: {project}")
    print(f"Input: {fasta_file}")
    print(f"Output: {output_dir}")
    if force_rerun:
        print("⚠️  Force rerun mode: All steps will be executed even if outputs exist")
    else:
        print("✅ Smart mode: Steps with existing outputs will be skipped")
    print("-" * 50)

    if taxon_name is None:
        if taxid is not None:
            taxon_name = get_taxon_name(taxid)
        if taxon_name is None:
            print(f"Invalid taxid or NCBI connection problem: {taxid}")
            print(f"Taxon set to 'Unknown'")
            taxon_name = "Unknown"
        else:
            print(f"Taxon name: {taxon_name}")
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
        genome_size = get_genome_size_with_progress(fasta_file)

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

    #TODO: use large_cutoff in code

    # Step 1: TRF Search - check if main TRF file exists and is not empty
    main_trf_file = f"{trf_prefix}.trf"
    if os.path.exists(main_trf_file) and os.path.getsize(main_trf_file) > 0 and not force_rerun:
        print(f"TRF search already completed! Found {main_trf_file} ({os.path.getsize(main_trf_file):,} bytes)")
        print("Use --force to rerun this step")
    else:
        if force_rerun and os.path.exists(main_trf_file):
            print("Force rerun: Running TRF search...")
        else:
            print("Running TRF search...")
        command = f"python {trf_search_path} -i {fasta_file} \
                                       -o {output_dir} \
                                       -p {project} \
                                       -t {threads} \
                                       --trf {trf_path} \
                                       --genome_size {genome_size}"
        
        # Add k-mer filtering options if enabled
        if use_kmer_filter or kmer_bed_file:
            command += " --use_kmer_filter"
            command += f" --kmer_threshold {kmer_threshold}"
            if kmer_bed_file:
                command += f" --kmer_bed {kmer_bed_file}"
        # print(command)
        completed_process = subprocess.run(command, shell=True)

        if completed_process.returncode == 0:
            print("trf_search.py executed successfully!")
        else:
            print(f"trf_search.py failed with return code {completed_process.returncode}")
            sys.exit(1)

    ### check for annotation
    was_annotated = False
    if os.path.exists(trf_file):
        for trf_obj in sc_iter_tab_file(trf_file, TRModel):
            if trf_obj.trf_ref_annotation is not None:
                was_annotated = True
            break

    if gff_file and not was_annotated:
        print("Adding annotation from GFF file...")
        reports_folder = os.path.join(
            output_dir,
            "reports")
        if not os.path.exists(reports_folder):
            os.makedirs(reports_folder)
        annotation_report_file = os.path.join(
            reports_folder,
            "annotation_report.txt"
        )
        add_annotation_from_gff(settings["trf_file"], gff_file, annotation_report_file, rm_file=repeatmasker_file)
        print("Annotation added!")
    else:
        if was_annotated:
            print("Annotation was added before!")
        else:
            print("Please provide GFF file and optionally RM file for annotation!")

    # Step 2: TRF Classification - check if classified files exist
    micro_trf_file = f"{trf_prefix}.micro.trf"
    complex_trf_file = f"{trf_prefix}.complex.trf"
    pmicro_trf_file = f"{trf_prefix}.pmicro.trf"
    tssr_trf_file = f"{trf_prefix}.tssr.trf"
    
    # Check if main classification files exist and are not empty
    # Note: some classification files can be empty if no repeats of that type found
    classification_complete = (
        os.path.exists(micro_trf_file) and 
        os.path.exists(complex_trf_file) and
        os.path.exists(pmicro_trf_file) and
        os.path.exists(tssr_trf_file)
    )
    
    if classification_complete and not force_rerun:
        print(f"TRF classification already completed! Found all classified files.")
        print("Use --force to rerun this step")
    else:
        if force_rerun and classification_complete:
            print("Force rerun: Running TRF classification...")
        else:
            print("Running TRF classification...")
        command = f"python {trf_classify_path} -i {trf_prefix} -o {output_dir} -l {genome_size}"
        # print(command)
        completed_process = subprocess.run(command, shell=True)
        if completed_process.returncode == 0:
            print("trf_classify.py executed successfully!")
        else:
            print(f"trf_classify.py failed with return code {completed_process.returncode}")
            sys.exit(1)


    # Step 3: TRF Drawing - check if distance files or images exist
    # Check for distance file with any extension (like .429, .429.vector)
    distance_files_exist = any(
        f.startswith("distances.tsv") for f in os.listdir(output_dir) 
        if os.path.isfile(os.path.join(output_dir, f))
    ) if os.path.exists(output_dir) else False
    
    html_report_exists = os.path.exists(html_report_file)
    
    if distance_files_exist and html_report_exists and not force_rerun:
        print(f"TRF drawing and HTML report already completed!")
        print("Use --force to rerun this step")
    else:
        if force_rerun and distance_files_exist:
            print("Force rerun: Running TRF drawing...")
        else:
            print("Running TRF drawing...")
        
        # Add --force flag if force_rerun is True
        force_flag = " --force" if force_rerun else ""
        command = f"python {trf_draw_path} -f {fasta_file} -i {trf_file} -o {output_image_dir} -c {minimal_scaffold_length} -e {drawing_enhancing} -t '{taxon_name}' -d {distance_file} -s {genome_size}{force_flag}"
        # print(command)
        completed_process = subprocess.run(command, shell=True)
        if completed_process.returncode == 0:
            print("trf_draw.py executed successfully!")
        else:
            print(f"trf_draw.py failed with return code {completed_process.returncode}")
            sys.exit(1)

        # Create HTML report only if drawing was successful
        create_html_report(output_image_dir, html_report_file)

    print("\n" + "="*50)
    print("SATELLOME ANALYSIS COMPLETED SUCCESSFULLY!")
    print("="*50)
    print(f"Project: {project}")
    print(f"Taxon: {taxon_name}")
    print(f"Output directory: {output_dir}")
    print(f"HTML report: {html_report_file}")
    print("="*50)


if __name__ == "__main__":
    main()
