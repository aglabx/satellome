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
import logging

try:
    from satellome import __version__
    from satellome.core_functions.io.tab_file import sc_iter_tab_file
    from satellome.core_functions.models.trf_model import TRModel
    from satellome.core_functions.tools.gene_intersect import add_annotation_from_gff
    from satellome.core_functions.tools.reports import create_html_report
    from satellome.core_functions.tools.processing import get_genome_size_with_progress
    from satellome.core_functions.tools.ncbi import get_taxon_name
    from satellome.installers import install_fastan, install_tanbed, install_trf_large
    from satellome.constants import (
        MIN_SCAFFOLD_LENGTH_DEFAULT, TR_CUTOFF_DEFAULT,
        KMER_THRESHOLD_DEFAULT, DRAWING_ENHANCING_DEFAULT,
        SEPARATOR_LINE, SEPARATOR_LINE_DOUBLE,
        DEFAULT_TAXON_NAME
    )
except ImportError:
    from src.satellome import __version__
    from src.satellome.core_functions.io.tab_file import sc_iter_tab_file
    from src.satellome.core_functions.models.trf_model import TRModel
    from src.satellome.core_functions.tools.gene_intersect import add_annotation_from_gff
    from src.satellome.core_functions.tools.reports import create_html_report
    from src.satellome.core_functions.tools.processing import get_genome_size_with_progress
    from src.satellome.core_functions.tools.ncbi import get_taxon_name
    from src.satellome.installers import install_fastan, install_tanbed, install_trf_large
    from src.satellome.constants import (
        MIN_SCAFFOLD_LENGTH_DEFAULT, TR_CUTOFF_DEFAULT,
        KMER_THRESHOLD_DEFAULT, DRAWING_ENHANCING_DEFAULT,
        SEPARATOR_LINE, SEPARATOR_LINE_DOUBLE,
        DEFAULT_TAXON_NAME
    )

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def print_logo():
    '''https://patorjk.com/software/taag/#p=display&f=Ghost&t=AGLABX%0Asatellome
    '''
    logo = '''
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
'''
    logger.info(logo)
    


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Satellome - Tandem Repeat Analysis Pipeline")

    # Version argument (handle it first)
    parser.add_argument(
        "-v", "--version",
        action="version",
        version=f"satellome v{__version__}",
        help="Show version information and exit"
    )

    parser.add_argument("-i", "--input", help="Input fasta file", required=False)
    parser.add_argument("-o", "--output", help="Output folder (must be an absolute path, e.g., /home/user/output)", required=False)
    parser.add_argument("-p", "--project", help="Project", required=False)
    parser.add_argument("-t", "--threads", help="Threads", required=False)
    parser.add_argument(
        "--trf", help="Path to TRF binary (default: trf in PATH)", required=False, default="trf"
    )
    parser.add_argument(
        "--genome_size", help="Expected genome size [will be computed from fasta]", required=False, default=0
    )
    parser.add_argument("--taxid", help="NCBI taxid, look here https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi [None]", required=False, default=None)
    parser.add_argument("--gff", help="Input gff file [None]", required=False, default=None)
    parser.add_argument("--rm", help="Input RM *.ori.out file [None]", required=False, default=None)
    parser.add_argument("--srr", help="SRR index for raw reads [None]", required=False, default=None)
    parser.add_argument("-c", "--cutoff", help=f"Cutoff for large TRs [{TR_CUTOFF_DEFAULT}]", required=False, default=TR_CUTOFF_DEFAULT)
    ### add minimal_scaffold_length
    parser.add_argument("-l", "--minimal_scaffold_length", help=f"Minimal scaffold length [{MIN_SCAFFOLD_LENGTH_DEFAULT}]", required=False, default=MIN_SCAFFOLD_LENGTH_DEFAULT)
    parser.add_argument("-e", "--drawing_enhancing", help=f"Drawing enhancing [{DRAWING_ENHANCING_DEFAULT}]", required=False, default=DRAWING_ENHANCING_DEFAULT)
    parser.add_argument("--large_file", help="Suffix for TR file for analysis, it can be '', 1kb, 3kb, 10kb [1kb]", required=False, default="1kb")
    parser.add_argument("--taxon", help="Taxon name [Unknown]", required=False, default=None)
    parser.add_argument("--force", help="Force rerun all steps even if output files exist", action='store_true', default=False)
    parser.add_argument("--recompute-failed", help="Recompute only chromosomes/contigs that failed TRF analysis (missing from TRF results)", action='store_true', default=False)
    parser.add_argument("--use_kmer_filter", help="Use k-mer profiling to filter repeat-poor regions", action='store_true', default=False)
    parser.add_argument("--kmer_threshold", help=f"Unique k-mer threshold for repeat detection [{KMER_THRESHOLD_DEFAULT}]", required=False, default=KMER_THRESHOLD_DEFAULT, type=int)
    parser.add_argument("--kmer_bed", help="Pre-computed k-mer profile BED file from varprofiler", required=False, default=None)
    parser.add_argument("--continue-on-error", help="Continue pipeline even if some TRF runs fail (results may be incomplete)", action='store_true', default=False)
    parser.add_argument("--keep-trf", help="Keep original TRF files before filtering (saved with .original suffix)", action='store_true', default=False)

    # Installation commands
    parser.add_argument("--install-fastan", help="Install FasTAN binary to ~/.satellome/bin/", action='store_true', default=False)
    parser.add_argument("--install-tanbed", help="Install tanbed binary to ~/.satellome/bin/", action='store_true', default=False)
    parser.add_argument("--install-trf-large", help="Install modified TRF (for large genomes) to ~/.satellome/bin/", action='store_true', default=False)
    parser.add_argument("--install-all", help="Install all external dependencies (FasTAN, tanbed, and modified TRF)", action='store_true', default=False)

    return vars(parser.parse_args())


def validate_and_prepare_environment(args):
    """Validate arguments and prepare the environment."""
    output_dir = args["output"]
    trf_path = args["trf"]

    # Convert relative path to absolute path if needed
    if not os.path.isabs(output_dir):
        original_path = output_dir
        output_dir = os.path.abspath(output_dir)
        logger.info(f"Converted relative path '{original_path}' to absolute path: {output_dir}")
        args["output"] = output_dir  # Update the args with the absolute path

    # Check if TRF is available
    import shutil
    if trf_path == "trf":
        trf_found = shutil.which(trf_path)
        if trf_found:
            logger.info(f"TRF binary: {trf_found}")
        else:
            logger.warning(f"TRF not found in PATH. Please install TRF or provide path with --trf")
            logger.warning("Download TRF from: https://tandem.bu.edu/trf/trf.html")
    else:
        # Check if the provided path exists
        if os.path.exists(trf_path) and os.access(trf_path, os.X_OK):
            logger.info(f"TRF binary: {trf_path}")
        else:
            logger.error(f"TRF not found or not executable at: {trf_path}")
            sys.exit(1)

    # Create necessary directories
    html_report_file = os.path.join(output_dir, "reports", "satellome_report.html")
    if not os.path.exists(os.path.dirname(html_report_file)):
        os.makedirs(os.path.dirname(html_report_file))

    output_image_dir = os.path.join(output_dir, "images")
    if not os.path.exists(output_image_dir):
        os.makedirs(output_image_dir)

    return html_report_file, output_image_dir


def build_settings(args, fasta_file, output_dir, project, threads, trf_path, genome_size, taxon_name, taxid, html_report_file, output_image_dir):
    """Build settings dictionary for the pipeline."""
    input_filename_without_extension = os.path.basename(os.path.splitext(fasta_file)[0])
    trf_prefix = os.path.join(output_dir, input_filename_without_extension)

    current_file_path = os.path.abspath(__file__)
    current_directory = os.path.dirname(current_file_path)

    distance_file = os.path.join(output_dir, "distances.tsv")

    return {
        "fasta_file": fasta_file,
        "output_dir": output_dir,
        "project": project,
        "threads": threads,
        "trf_path": trf_path,
        "genome_size": genome_size,
        "trf_prefix": trf_prefix,
        "large_cutoff": int(args["cutoff"]),
        "trf_search_path": os.path.join(current_directory, "steps", "trf_search.py"),
        "trf_classify_path": os.path.join(current_directory, "steps", "trf_classify.py"),
        "trf_draw_path": os.path.join(current_directory, "steps", "trf_draw.py"),
        "trf_parse_raw_path": os.path.join(current_directory, "steps", "trf_parse_raw.py"),
        "gff_file": args["gff"],
        "trf_file": f"{trf_prefix}.trf",
        "minimal_scaffold_length": int(args["minimal_scaffold_length"]),
        "drawing_enhancing": int(args["drawing_enhancing"]),
        "taxon_name": taxon_name,
        "srr": args["srr"],
        "taxid": taxid,
        "distance_file": distance_file,
        "output_image_dir": output_image_dir,
        "large_file_suffix": args["large_file"],
        "repeatmasker_file": args["rm"],
        "html_report_file": html_report_file,
    }


def run_trf_search(settings, args, force_rerun):
    """Run TRF search step."""
    from satellome.core_functions.tools.trf_tools import recompute_failed_chromosomes

    trf_prefix = settings["trf_prefix"]
    main_trf_file = f"{trf_prefix}.trf"

    # Check if recompute-failed mode is enabled
    if args.get("recompute_failed", False) and os.path.exists(main_trf_file) and os.path.getsize(main_trf_file) > 0:
        logger.info("RECOMPUTE-FAILED MODE: Will check and recompute only missing chromosomes...")

        try:
            success = recompute_failed_chromosomes(
                fasta_file=settings['fasta_file'],
                existing_trf_file=main_trf_file,
                output_dir=settings['output_dir'],
                project=settings['project'],
                threads=settings['threads'],
                trf_path=settings['trf_path'],
                parser_program=settings['trf_parse_raw_path'],
                min_scaffold_size=1000000,
                match_first_word=True,
            )

            if success:
                logger.info("✅ Recompute completed successfully!")
                logger.info("✅ TRF file updated. Will proceed to regenerate downstream files (1kb, 3kb, 10kb, images, report)...")
                return "recomputed"
            else:
                logger.error("❌ Recompute failed")
                sys.exit(1)

        except Exception as e:
            logger.error(f"❌ Recompute failed with error: {e}")
            sys.exit(1)

    if os.path.exists(main_trf_file) and os.path.getsize(main_trf_file) > 0 and not force_rerun:
        logger.info(f"TRF search already completed! Found {main_trf_file} ({os.path.getsize(main_trf_file):,} bytes)")
        logger.info("Use --force to rerun this step or --recompute-failed to recompute only missing chromosomes")
        return True

    if force_rerun and os.path.exists(main_trf_file):
        logger.info("Force rerun: Running TRF search...")
    else:
        logger.info("Running TRF search...")

    command = f"python {settings['trf_search_path']} -i {settings['fasta_file']} \
                                   -o {settings['output_dir']} \
                                   -p {settings['project']} \
                                   -t {settings['threads']} \
                                   --trf {settings['trf_path']} \
                                   --genome_size {settings['genome_size']}"

    # Add k-mer filtering options if enabled
    if args["use_kmer_filter"] or args["kmer_bed"]:
        command += " --use_kmer_filter"
        command += f" --kmer_threshold {args['kmer_threshold']}"
        if args["kmer_bed"]:
            command += f" --kmer_bed {args['kmer_bed']}"

    # Add continue-on-error option if enabled
    if args["continue_on_error"]:
        command += " --continue-on-error"

    logger.debug(f"Command: {command}")
    completed_process = subprocess.run(command, shell=True)

    if completed_process.returncode == 0:
        logger.info("trf_search.py executed successfully!")
        return True
    else:
        logger.error(f"trf_search.py failed with return code {completed_process.returncode}")
        sys.exit(1)


def add_annotations(settings, force_rerun):
    """Add annotations from GFF and RepeatMasker files."""
    trf_file = settings["trf_file"]
    if settings["large_file_suffix"]:
        trf_file = f"{settings['trf_prefix']}.{settings['large_file_suffix']}.trf"

    # Check if already annotated
    was_annotated = False
    if os.path.exists(trf_file):
        for trf_obj in sc_iter_tab_file(trf_file, TRModel):
            if trf_obj.trf_ref_annotation is not None:
                was_annotated = True
            break

    if settings["gff_file"] and not was_annotated:
        logger.info("Adding annotation from GFF file...")
        reports_folder = os.path.join(settings["output_dir"], "reports")
        if not os.path.exists(reports_folder):
            os.makedirs(reports_folder)
        annotation_report_file = os.path.join(reports_folder, "annotation_report.txt")
        add_annotation_from_gff(
            settings["trf_file"],
            settings["gff_file"],
            annotation_report_file,
            rm_file=settings["repeatmasker_file"]
        )
        logger.info("Annotation added!")
    else:
        if was_annotated:
            logger.info("Annotation was added before!")
        else:
            logger.info("Please provide GFF file and optionally RM file for annotation!")


def run_trf_classification(settings, args, force_rerun):
    """Run TRF classification step."""
    trf_prefix = settings["trf_prefix"]

    # Check if main classification files exist
    classification_files = [
        f"{trf_prefix}.micro.trf",
        f"{trf_prefix}.complex.trf",
        f"{trf_prefix}.pmicro.trf",
        f"{trf_prefix}.tssr.trf"
    ]

    classification_complete = all(os.path.exists(f) for f in classification_files)

    if classification_complete and not force_rerun:
        logger.info(f"TRF classification already completed! Found all classified files.")
        logger.info("Use --force to rerun this step")
        return True

    if force_rerun and classification_complete:
        logger.info("Force rerun: Running TRF classification...")
    else:
        logger.info("Running TRF classification...")

    command = f"python {settings['trf_classify_path']} -i {trf_prefix} -o {settings['output_dir']} -l {settings['genome_size']}"
    if args["keep_trf"]:
        command += " --keep-trf"

    logger.debug(f"Command: {command}")
    completed_process = subprocess.run(command, shell=True)

    if completed_process.returncode == 0:
        logger.info("trf_classify.py executed successfully!")
        return True
    else:
        logger.error(f"trf_classify.py failed with return code {completed_process.returncode}")
        sys.exit(1)


def run_trf_drawing(settings, force_rerun):
    """Run TRF drawing and report generation step."""
    output_dir = settings["output_dir"]
    html_report_file = settings["html_report_file"]

    # Check for distance file with any extension
    distance_files_exist = any(
        f.startswith("distances.tsv") for f in os.listdir(output_dir)
        if os.path.isfile(os.path.join(output_dir, f))
    ) if os.path.exists(output_dir) else False

    html_report_exists = os.path.exists(html_report_file)

    if distance_files_exist and html_report_exists and not force_rerun:
        logger.info(f"TRF drawing and HTML report already completed!")
        logger.info("Use --force to rerun this step")
        return True

    if force_rerun and distance_files_exist:
        logger.info("Force rerun: Running TRF drawing...")
    else:
        logger.info("Running TRF drawing...")

    # Build TRF file path with suffix
    trf_file = settings["trf_file"]
    if settings["large_file_suffix"]:
        trf_file = f"{settings['trf_prefix']}.{settings['large_file_suffix']}.trf"

    # Add --force flag if force_rerun is True
    force_flag = " --force" if force_rerun else ""
    command = f"python {settings['trf_draw_path']} -f {settings['fasta_file']} -i {trf_file} -o {settings['output_image_dir']} -c {settings['minimal_scaffold_length']} -e {settings['drawing_enhancing']} -t '{settings['taxon_name']}' -s {settings['genome_size']}{force_flag}"

    logger.debug(f"Command: {command}")
    completed_process = subprocess.run(command, shell=True)

    if completed_process.returncode == 0:
        logger.info("trf_draw.py executed successfully!")
        # Create HTML report only if drawing was successful
        create_html_report(settings["output_image_dir"], html_report_file)
        return True
    else:
        logger.error(f"trf_draw.py failed with return code {completed_process.returncode}")
        sys.exit(1)


def handle_installation_commands(args):
    """
    Handle installation commands (--install-fastan, --install-tanbed, --install-trf-large, --install-all).

    Args:
        args: Parsed command line arguments

    Returns:
        bool: True if installation commands were processed (and program should exit), False otherwise
    """
    install_fastan_flag = args.get("install_fastan", False)
    install_tanbed_flag = args.get("install_tanbed", False)
    install_trf_large_flag = args.get("install_trf_large", False)
    install_all_flag = args.get("install_all", False)

    # If no installation commands, return False to continue with main pipeline
    if not (install_fastan_flag or install_tanbed_flag or install_trf_large_flag or install_all_flag):
        return False

    logger.info(SEPARATOR_LINE_DOUBLE)
    logger.info("Installation mode activated")
    logger.info(SEPARATOR_LINE_DOUBLE)

    success = True

    # Install FasTAN
    if install_fastan_flag or install_all_flag:
        logger.info("Installing FasTAN...")
        if install_fastan(force=True):
            logger.info("✓ FasTAN installed successfully")
        else:
            logger.error("✗ FasTAN installation failed")
            success = False
        logger.info(SEPARATOR_LINE)

    # Install tanbed
    if install_tanbed_flag or install_all_flag:
        logger.info("Installing tanbed...")
        if install_tanbed(force=True):
            logger.info("✓ tanbed installed successfully")
        else:
            logger.error("✗ tanbed installation failed")
            success = False
        logger.info(SEPARATOR_LINE)

    # Install modified TRF (for large genomes)
    if install_trf_large_flag or install_all_flag:
        logger.info("Installing modified TRF (for large genomes)...")
        if install_trf_large(force=True):
            logger.info("✓ Modified TRF installed successfully")
        else:
            logger.error("✗ Modified TRF installation failed")
            success = False
        logger.info(SEPARATOR_LINE)

    # Print summary
    logger.info(SEPARATOR_LINE_DOUBLE)
    if success:
        logger.info("All installations completed successfully!")
        logger.info("Binaries installed to: ~/.satellome/bin/")
        logger.info("You can now use these tools with Satellome.")
    else:
        logger.error("Some installations failed. Please check the error messages above.")
        sys.exit(1)

    logger.info(SEPARATOR_LINE_DOUBLE)

    # Return True to indicate program should exit after installation
    return True


def print_summary(project, taxon_name, output_dir, html_report_file):
    """Print final summary of the analysis."""
    logger.info("\n" + SEPARATOR_LINE_DOUBLE)
    logger.info("SATELLOME ANALYSIS COMPLETED SUCCESSFULLY!")
    logger.info(SEPARATOR_LINE_DOUBLE)
    logger.info(f"Project: {project}")
    logger.info(f"Taxon: {taxon_name}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"HTML report: {html_report_file}")
    logger.info(SEPARATOR_LINE_DOUBLE)


def main():
    args = parse_arguments()

    # Handle installation commands first (exits if installation was performed)
    if handle_installation_commands(args):
        sys.exit(0)

    # Validate required arguments for pipeline mode
    required_args = ["input", "output", "project", "threads"]
    missing_args = [arg for arg in required_args if not args.get(arg)]
    if missing_args:
        logger.error(f"Missing required arguments: {', '.join(missing_args)}")
        logger.error("Use --help to see all required arguments")
        sys.exit(1)

    print_logo()

    # Validate and prepare environment
    html_report_file, output_image_dir = validate_and_prepare_environment(args)

    # Extract main parameters
    fasta_file = args["input"]
    output_dir = args["output"]
    project = args["project"]
    threads = args["threads"]
    trf_path = args["trf"]
    genome_size = int(args["genome_size"])
    taxid = args["taxid"]
    taxon_name = args["taxon"]
    force_rerun = args["force"]

    logger.info(f"Starting Satellome analysis...")
    logger.info(f"Project: {project}")
    logger.info(f"Input: {fasta_file}")
    logger.info(f"Output: {output_dir}")

    if force_rerun:
        logger.warning("Force rerun mode: All steps will be executed even if outputs exist")
    else:
        logger.info("Smart mode: Steps with existing outputs will be skipped")
    logger.info(SEPARATOR_LINE)

    # Resolve taxon name
    if taxon_name is None:
        if taxid is not None:
            taxon_name = get_taxon_name(taxid)
        if taxon_name is None:
            logger.warning(f"Invalid taxid or NCBI connection problem: {taxid}")
            logger.warning(f"Taxon set to 'Unknown'")
            taxon_name = DEFAULT_TAXON_NAME
        else:
            logger.info(f"Taxon name: {taxon_name}")
    taxon_name = taxon_name.replace(" ", "_")

    # Calculate genome size if needed
    if not genome_size:
        genome_size = get_genome_size_with_progress(fasta_file)

    # Build settings
    settings = build_settings(
        args, fasta_file, output_dir, project, threads, trf_path,
        genome_size, taxon_name, taxid, html_report_file, output_image_dir
    )

    #TODO: use large_cutoff in code

    # Step 1: TRF Search
    trf_search_result = run_trf_search(settings, args, force_rerun)

    # If recompute-failed mode was used and TRF was updated, force regeneration of downstream files
    force_downstream = force_rerun or (trf_search_result == "recomputed")

    # Step 2: Add annotations
    add_annotations(settings, force_downstream)

    # Step 3: TRF Classification
    run_trf_classification(settings, args, force_downstream)

    # Step 4: TRF Drawing and HTML report
    run_trf_drawing(settings, force_downstream)

    # Print summary
    print_summary(project, taxon_name, output_dir, html_report_file)


if __name__ == "__main__":
    main()
