#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 26.10.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import sys
import os
import argparse
import subprocess
import logging

from satellome import __version__
from satellome.core_functions.io.tab_file import sc_iter_tab_file
from satellome.core_functions.models.trf_model import TRModel
from satellome.core_functions.tools.gene_intersect import add_annotation_from_gff
from satellome.core_functions.tools.reports import create_html_report
from satellome.core_functions.tools.processing import get_genome_size_with_progress
from satellome.core_functions.tools.ncbi import get_taxon_name
from satellome.core_functions.tools.bed_tools import extract_sequences_from_bed
from satellome.core_functions.tools.validation import (
    validate_input_files, validate_fasta_file, validate_gff_file,
    validate_repeatmasker_file, validate_trf_binary, validate_output_directory,
    ValidationError, FastaValidationError, GFFValidationError,
    BinaryValidationError, OutputDirValidationError
)
from satellome.installers import install_fastan, install_tanbed, install_trf_large, install_trf_standard
from satellome.constants import (
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
    logo = r'''
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
    parser.add_argument("--nofastan", help="Skip FasTAN analysis", action='store_true', default=False)
    parser.add_argument("--run-trf", help="Run TRF analysis (disabled by default, FasTAN is the default tool)", action='store_true', default=False)
    parser.add_argument("--notrf", help="[DEPRECATED] TRF is now disabled by default. Use --run-trf to enable.", action='store_true', default=False)

    # Installation commands
    parser.add_argument("--install-fastan", help="Install FasTAN binary to ~/.satellome/bin/", action='store_true', default=False)
    parser.add_argument("--install-tanbed", help="Install tanbed binary to ~/.satellome/bin/", action='store_true', default=False)
    parser.add_argument("--install-trf", help="Install standard TRF (download pre-compiled binary)", action='store_true', default=False)
    parser.add_argument("--install-trf-large", help="Install modified TRF (for large genomes) to ~/.satellome/bin/", action='store_true', default=False)
    parser.add_argument("--install-all", help="Install all external dependencies (FasTAN, tanbed, and modified TRF)", action='store_true', default=False)

    return vars(parser.parse_args())


def validate_and_prepare_environment(args):
    """Validate arguments and prepare the environment."""
    output_dir = args["output"]
    trf_path = args["trf"]
    fasta_file = args["input"]
    gff_file = args.get("gff")
    rm_file = args.get("rm")

    # Convert relative paths to absolute paths if needed
    if not os.path.isabs(output_dir):
        output_dir = os.path.abspath(output_dir)
        args["output"] = output_dir

    if not os.path.isabs(fasta_file):
        fasta_file = os.path.abspath(fasta_file)
        args["input"] = fasta_file

    if gff_file and not os.path.isabs(gff_file):
        gff_file = os.path.abspath(gff_file)
        args["gff"] = gff_file

    if rm_file and not os.path.isabs(rm_file):
        rm_file = os.path.abspath(rm_file)
        args["rm"] = rm_file

    logger.info(SEPARATOR_LINE)
    logger.info("INPUT VALIDATION")
    logger.info(SEPARATOR_LINE)

    # Validate FASTA file first (required)
    try:
        fasta_stats = validate_fasta_file(fasta_file, check_sequences=True)
        logger.info(
            f"✓ FASTA: {fasta_stats['num_sequences']} sequences, "
            f"{fasta_stats['total_length']:,} bp total"
        )
        if fasta_stats['warnings']:
            logger.warning(f"FASTA has {len(fasta_stats['warnings'])} warnings (see details below)")
            for warning in fasta_stats['warnings'][:3]:
                logger.warning(f"  - {warning}")
            if len(fasta_stats['warnings']) > 3:
                logger.warning(f"  ... and {len(fasta_stats['warnings']) - 3} more warnings")
    except FastaValidationError as e:
        logger.error(f"✗ FASTA validation failed: {e}")
        sys.exit(1)

    # Validate GFF file if provided
    if gff_file:
        try:
            gff_stats = validate_gff_file(gff_file)
            logger.info(f"✓ GFF: {gff_stats['num_features']} features")
            if gff_stats['num_malformed'] > 0:
                logger.warning(
                    f"GFF has {gff_stats['num_malformed']} malformed lines (will be skipped)"
                )
        except GFFValidationError as e:
            logger.error(f"✗ GFF validation failed: {e}")
            sys.exit(1)

    # Validate RepeatMasker file if provided
    if rm_file:
        try:
            rm_stats = validate_repeatmasker_file(rm_file)
            logger.info(f"✓ RepeatMasker: {rm_stats['num_features']} features")
            if rm_stats['num_malformed'] > 0:
                logger.warning(
                    f"RepeatMasker has {rm_stats['num_malformed']} malformed lines (will be skipped)"
                )
        except ValidationError as e:
            logger.error(f"✗ RepeatMasker validation failed: {e}")
            sys.exit(1)

    # Check if TRF is available (only if TRF is enabled)
    run_trf = args.get("run_trf", False)
    if run_trf:
        import shutil
        try:
            trf_found = validate_trf_binary(trf_path)
            logger.info(f"✓ TRF binary: {trf_found}")
            args["trf"] = trf_found  # Update with full path
        except BinaryValidationError:
            logger.warning(f"TRF not found: {trf_path}")
            logger.info("Attempting to install TRF automatically...")

            # Try to auto-install TRF
            try:
                from satellome.installers import install_trf_large, install_trf_standard
                from satellome.installers.base import get_satellome_bin_dir

                # Try modified TRF first (for large genomes)
                logger.info("Trying modified TRF (for large genomes >2GB chromosomes)...")
                if install_trf_large(force=False):
                    logger.info("✓ Modified TRF installed successfully!")
                    trf_bin = get_satellome_bin_dir() / "trf"
                    if trf_bin.exists():
                        trf_path = str(trf_bin)
                        args["trf"] = trf_path
                        logger.info(f"Using installed TRF: {trf_path}")
                    else:
                        logger.error("TRF installation succeeded but binary not found")
                        sys.exit(1)
                else:
                    # Fallback to standard TRF (download pre-compiled binary)
                    logger.warning("Modified TRF installation failed (missing build tools?)")
                    logger.info("Falling back to standard TRF (download pre-compiled binary)...")

                    if install_trf_standard(force=False):
                        logger.info("✓ Standard TRF installed successfully!")
                        trf_bin = get_satellome_bin_dir() / "trf"
                        if trf_bin.exists():
                            trf_path = str(trf_bin)
                            args["trf"] = trf_path
                            logger.info(f"Using installed TRF: {trf_path}")
                            logger.info("Note: Using standard TRF. For genomes with chromosomes >2GB,")
                            logger.info("      install build tools and run: satellome --install-trf-large")
                        else:
                            logger.error("TRF installation succeeded but binary not found")
                            sys.exit(1)
                    else:
                        logger.error("Both TRF installers failed")
                        logger.warning("Please install TRF manually:")
                        logger.warning("  Option 1: satellome --install-trf-large (requires build tools)")
                        logger.warning("  Option 2: satellome --install-trf (pre-compiled binary)")
                        logger.warning("  Option 3: Download from https://tandem.bu.edu/trf/trf.html")
                        sys.exit(1)
            except (OSError, IOError, PermissionError) as e:
                logger.error(f"TRF auto-installation failed (I/O error): {e}")
                logger.warning("Please install TRF manually:")
                logger.warning("  Option 1: satellome --install-trf-large (requires build tools)")
                logger.warning("  Option 2: satellome --install-trf (pre-compiled binary)")
                logger.warning("  Option 3: Download from https://tandem.bu.edu/trf/trf.html")
                sys.exit(1)
            except Exception as e:
                # Catch unexpected errors but log them distinctly
                logger.error(f"TRF auto-installation failed (unexpected error): {type(e).__name__}: {e}")
                logger.warning("Please install TRF manually:")
                logger.warning("  Option 1: satellome --install-trf-large (requires build tools)")
                logger.warning("  Option 2: satellome --install-trf (pre-compiled binary)")
                logger.warning("  Option 3: Download from https://tandem.bu.edu/trf/trf.html")
                sys.exit(1)
    else:
        logger.info("✓ TRF validation skipped (TRF disabled by default, use --run-trf to enable)")

    # Validate output directory
    try:
        output_path = validate_output_directory(output_dir, create_if_missing=True)
        logger.info(f"✓ Output directory: {output_path}")
        args["output"] = output_path  # Update with absolute path
    except OutputDirValidationError as e:
        logger.error(f"✗ Output directory validation failed: {e}")
        sys.exit(1)

    logger.info("All input validations passed!")
    logger.info(SEPARATOR_LINE)

    # Create necessary subdirectories
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
        "trf_file": f"{trf_prefix}.sat",
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
    main_trf_file = f"{trf_prefix}.sat"

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

    command = f"{sys.executable} {settings['trf_search_path']} -i {settings['fasta_file']} \
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
        trf_file = f"{settings['trf_prefix']}.{settings['large_file_suffix']}.sat"

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
    # Use directory containing trf_prefix for output
    classify_output_dir = os.path.dirname(trf_prefix)

    # Check if main classification files exist
    classification_files = [
        f"{trf_prefix}.micro.sat",
        f"{trf_prefix}.complex.sat",
        f"{trf_prefix}.pmicro.sat",
        f"{trf_prefix}.tssr.sat"
    ]

    classification_complete = all(os.path.exists(f) for f in classification_files)

    if classification_complete and not force_rerun:
        logger.info(f"Classification already completed! Found all classified files.")
        logger.info("Use --force to rerun this step")
        return True

    if force_rerun and classification_complete:
        logger.info("Force rerun: Running classification...")
    else:
        logger.info("Running classification...")

    command = f"{sys.executable} {settings['trf_classify_path']} -i {trf_prefix} -o {classify_output_dir} -l {settings['genome_size']}"
    if args["keep_trf"]:
        command += " --keep-trf"

    logger.debug(f"Command: {command}")
    completed_process = subprocess.run(command, shell=True)

    if completed_process.returncode == 0:
        logger.info("Classification completed successfully!")
        return True
    else:
        logger.error(f"Classification failed with return code {completed_process.returncode}")
        sys.exit(1)


def run_trf_drawing(settings, force_rerun):
    """Run TRF drawing and report generation step."""
    # Use directory containing trf_prefix for checking outputs
    drawing_output_dir = os.path.dirname(settings["trf_prefix"])
    html_report_file = settings["html_report_file"]

    # Check for distance file with any extension
    distance_files_exist = any(
        f.startswith("distances.tsv") for f in os.listdir(drawing_output_dir)
        if os.path.isfile(os.path.join(drawing_output_dir, f))
    ) if os.path.exists(drawing_output_dir) else False

    html_report_exists = os.path.exists(html_report_file)

    if distance_files_exist and html_report_exists and not force_rerun:
        logger.info(f"Drawing and HTML report already completed!")
        logger.info("Use --force to rerun this step")
        return True

    if force_rerun and distance_files_exist:
        logger.info("Force rerun: Running drawing...")
    else:
        logger.info("Running drawing...")

    # Build TRF file path with suffix
    trf_file = settings["trf_file"]
    if settings["large_file_suffix"]:
        trf_file = f"{settings['trf_prefix']}.{settings['large_file_suffix']}.sat"

    # Add --force flag if force_rerun is True
    force_flag = " --force" if force_rerun else ""
    command = f"{sys.executable} {settings['trf_draw_path']} -f {settings['fasta_file']} -i {trf_file} -o {settings['output_image_dir']} -c {settings['minimal_scaffold_length']} -e {settings['drawing_enhancing']} -t '{settings['taxon_name']}' -s {settings['genome_size']}{force_flag}"

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


def run_fastan(settings, force_rerun):
    """Run FasTAN analysis step."""
    import shutil
    from pathlib import Path
    from satellome.installers.base import get_satellome_bin_dir

    fasta_file = settings["fasta_file"]
    output_dir = settings["output_dir"]
    project = settings["project"]

    # Create output directories
    # fastan/ - intermediate files only (.1aln, .bed, ArraySplitter outputs)
    # fasta/ - output FASTA files
    # gff3/ - output GFF3 files
    # images/ - visualizations
    # reports/ - HTML reports
    # *.sat files go at output_dir level
    fastan_dir = os.path.join(output_dir, "fastan")
    fasta_dir = os.path.join(output_dir, "fasta")
    gff3_dir = os.path.join(output_dir, "gff3")

    for dir_path in [fastan_dir, fasta_dir, gff3_dir]:
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
            logger.info(f"Created directory: {dir_path}")

    # Output files - use genome filename instead of project name
    genome_basename = os.path.splitext(os.path.basename(fasta_file))[0]
    # Remove .gz extension if present
    if genome_basename.endswith('.gz'):
        genome_basename = os.path.splitext(genome_basename)[0]

    # Intermediate files in fastan/
    aln_file = os.path.join(fastan_dir, f"{genome_basename}.1aln")
    bed_file = os.path.join(fastan_dir, f"{genome_basename}.bed")

    # Output files at output_dir level or in dedicated subdirs
    trf_file = os.path.join(output_dir, f"{genome_basename}.sat")
    fasta_output = os.path.join(fasta_dir, f"{genome_basename}.arrays.fasta")

    # Check if already completed (all main output files exist)
    if not force_rerun:
        existing_files = []
        missing_files = []
        for f in [aln_file, bed_file, trf_file, fasta_output]:
            if os.path.exists(f):
                existing_files.append(os.path.basename(f))
            else:
                missing_files.append(os.path.basename(f))

        if not missing_files:
            logger.info("FasTAN analysis already completed!")
            logger.info(f"  Output directory: {output_dir}")
            logger.info(f"  Found files: {', '.join(existing_files)}")
            # Check for size-filtered files (at output_dir level)
            for suffix in ["1kb", "3kb", "10kb"]:
                filtered_trf = os.path.join(output_dir, f"{genome_basename}.{suffix}.sat")
                if os.path.exists(filtered_trf):
                    logger.info(f"  Found filtered: {genome_basename}.{suffix}.sat")
            # Check for ArraySplitter output (in fastan/ intermediate dir)
            hors_file = os.path.join(fastan_dir, f"{genome_basename}.hors.tsv")
            if os.path.exists(hors_file):
                logger.info(f"  Found ArraySplitter: {genome_basename}.hors.tsv")
            logger.info("Use --force to rerun this step")
            return True
        elif existing_files:
            logger.info(f"Partial FasTAN results found. Missing: {', '.join(missing_files)}")
            logger.info("Continuing analysis...")

    # Find fastan binary (auto-install if not found)
    fastan_bin = shutil.which("fastan")
    if not fastan_bin:
        # Try satellome bin directory
        satellome_bin = get_satellome_bin_dir() / "fastan"
        if satellome_bin.exists():
            fastan_bin = str(satellome_bin)
        else:
            # Auto-install FasTAN
            logger.warning("FasTAN binary not found.")
            logger.info("Attempting to install FasTAN automatically...")
            try:
                if install_fastan(force=False):
                    logger.info("✓ FasTAN installed successfully!")
                    fastan_bin = str(get_satellome_bin_dir() / "fastan")
                    if not os.path.exists(fastan_bin):
                        logger.error("FasTAN installation succeeded but binary not found")
                        return False
                else:
                    logger.warning("FasTAN installation failed (missing build tools?)")
                    logger.warning("Install manually with: satellome --install-fastan")
                    logger.warning("Skipping FasTAN analysis...")
                    return False
            except Exception as e:
                logger.warning(f"FasTAN auto-installation failed: {e}")
                logger.warning("Install manually with: satellome --install-fastan")
                logger.warning("Skipping FasTAN analysis...")
                return False

    logger.info(f"FasTAN binary: {fastan_bin}")

    # Find tanbed binary (auto-install if not found)
    tanbed_bin = shutil.which("tanbed")
    if not tanbed_bin:
        # Try satellome bin directory
        satellome_bin = get_satellome_bin_dir() / "tanbed"
        if satellome_bin.exists():
            tanbed_bin = str(satellome_bin)
        else:
            # Auto-install tanbed
            logger.warning("tanbed binary not found.")
            logger.info("Attempting to install tanbed automatically...")
            try:
                if install_tanbed(force=False):
                    logger.info("✓ tanbed installed successfully!")
                    tanbed_bin = str(get_satellome_bin_dir() / "tanbed")
                    if not os.path.exists(tanbed_bin):
                        logger.error("tanbed installation succeeded but binary not found")
                        return False
                else:
                    logger.warning("tanbed installation failed (missing build tools?)")
                    logger.warning("Install manually with: satellome --install-tanbed")
                    logger.warning("Skipping FasTAN analysis...")
                    return False
            except Exception as e:
                logger.warning(f"tanbed auto-installation failed: {e}")
                logger.warning("Install manually with: satellome --install-tanbed")
                logger.warning("Skipping FasTAN analysis...")
                return False

    logger.info(f"tanbed binary: {tanbed_bin}")

    # Run FasTAN
    if force_rerun and os.path.exists(aln_file):
        logger.info("Force rerun: Running FasTAN...")
    else:
        logger.info("Running FasTAN...")

    fastan_command = f"{fastan_bin} {fasta_file} {aln_file}"
    logger.debug(f"Command: {fastan_command}")

    try:
        # Don't capture output so progress is visible
        fastan_process = subprocess.run(fastan_command, shell=True)

        if fastan_process.returncode == 0:
            logger.info("FasTAN executed successfully!")
        else:
            logger.error(f"FasTAN failed with return code {fastan_process.returncode}")
            return False
    except Exception as e:
        logger.error(f"FasTAN execution failed: {e}")
        return False

    # Run tanbed to convert to BED format
    logger.info("Converting FasTAN output to BED format...")
    tanbed_command = f"{tanbed_bin} {aln_file} > {bed_file}"
    logger.debug(f"Command: {tanbed_command}")

    try:
        tanbed_process = subprocess.run(tanbed_command, shell=True, capture_output=True, text=True)

        if tanbed_process.returncode == 0:
            logger.info(f"✓ BED file created: {bed_file}")

            # Extract sequences from FASTA based on BED coordinates
            logger.info("Extracting sequences from FASTA based on BED coordinates...")
            try:
                extracted_count = extract_sequences_from_bed(
                    fasta_file, bed_file, trf_file,
                    fasta_output_file=fasta_output,
                    project=project
                )
                logger.info(f"✓ Sequence extraction completed: {trf_file}")
                logger.info(f"✓ FASTA output created: {fasta_output}")
                logger.info(f"✓ Extracted {extracted_count} sequences from {os.path.basename(fasta_file)}")

                # Create size-filtered TRF files (1kb, 3kb, 10kb)
                from satellome.core_functions.tools.bed_tools import filter_trf_by_size

                size_cutoffs = [
                    (1000, "1kb"),
                    (3000, "3kb"),
                    (10000, "10kb"),
                ]

                logger.info("Creating size-filtered TRF files...")
                for cutoff, suffix in size_cutoffs:
                    # .sat files at output_dir level, .fasta files in fasta/ subdir
                    filtered_trf = os.path.join(output_dir, f"{genome_basename}.{suffix}.sat")
                    filtered_fasta = os.path.join(fasta_dir, f"{genome_basename}.{suffix}.arrays.fasta")
                    stats = filter_trf_by_size(trf_file, filtered_trf, cutoff, fasta_output_file=filtered_fasta)
                    logger.info(f"✓ {suffix}: {stats['filtered']} arrays > {cutoff} bp")

            except Exception as e:
                logger.error(f"Sequence extraction failed: {e}")
                logger.warning("Continuing without sequence extraction...")
                # Don't fail the whole pipeline - BED file is still useful

            # Run ArraySplitter on main FASTA to calculate consensus and HORs
            if os.path.exists(fasta_output) and os.path.getsize(fasta_output) > 0:
                logger.info("Running ArraySplitter for consensus calculation...")
                threads = settings.get("threads", 4)
                arraysplitter_prefix = os.path.join(fastan_dir, genome_basename)
                run_arraysplitter(fasta_output, arraysplitter_prefix, threads, force_rerun)

            logger.info(f"✓ FasTAN analysis completed!")
            return True
        else:
            logger.error(f"tanbed failed with return code {tanbed_process.returncode}")
            logger.error(f"Error: {tanbed_process.stderr}")
            return False
    except Exception as e:
        logger.error(f"tanbed execution failed: {e}")
        return False


def run_arraysplitter(fasta_file, output_prefix, threads, force_rerun=False):
    """
    Run ArraySplitter to calculate consensus sequences and HORs.

    Args:
        fasta_file: Input FASTA file with arrays
        output_prefix: Output prefix for ArraySplitter results
        threads: Number of threads to use
        force_rerun: Force rerun even if output exists

    Returns:
        bool: True if successful, False otherwise
    """
    import shutil

    # Check if output already exists
    hors_file = f"{output_prefix}.hors.tsv"
    if os.path.exists(hors_file) and not force_rerun:
        logger.info(f"ArraySplitter output already exists: {hors_file}")
        return True

    # Check if input file exists and has content
    if not os.path.exists(fasta_file) or os.path.getsize(fasta_file) == 0:
        logger.warning(f"Input FASTA file not found or empty: {fasta_file}")
        return False

    # Find arraysplitter binary
    arraysplitter_bin = shutil.which("arraysplitter")
    if not arraysplitter_bin:
        # Try to install via pip
        logger.warning("arraysplitter not found. Installing via pip...")
        try:
            install_result = subprocess.run(
                [sys.executable, "-m", "pip", "install", "arraysplitter"],
                capture_output=True,
                text=True,
                timeout=120
            )
            if install_result.returncode == 0:
                logger.info("✓ arraysplitter installed successfully!")
                arraysplitter_bin = shutil.which("arraysplitter")
                if not arraysplitter_bin:
                    logger.error("arraysplitter installed but binary not found in PATH")
                    return False
            else:
                logger.error(f"Failed to install arraysplitter: {install_result.stderr}")
                return False
        except Exception as e:
            logger.error(f"Failed to install arraysplitter: {e}")
            return False

    # Run ArraySplitter
    logger.info(f"Running ArraySplitter on {os.path.basename(fasta_file)}...")
    command = f"{arraysplitter_bin} -i {fasta_file} -o {output_prefix} -t {threads}"
    logger.debug(f"Command: {command}")

    try:
        # Don't capture output so progress bars are visible
        result = subprocess.run(command, shell=True, timeout=3600)

        if result.returncode == 0:
            logger.info(f"✓ ArraySplitter completed: {output_prefix}")
            return True
        else:
            logger.error(f"ArraySplitter failed with return code {result.returncode}")
            return False
    except subprocess.TimeoutExpired:
        logger.error("ArraySplitter timed out (1 hour limit)")
        return False
    except Exception as e:
        logger.error(f"ArraySplitter execution failed: {e}")
        return False


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
    install_trf_flag = args.get("install_trf", False)
    install_trf_large_flag = args.get("install_trf_large", False)
    install_all_flag = args.get("install_all", False)

    # If no installation commands, return False to continue with main pipeline
    if not (install_fastan_flag or install_tanbed_flag or install_trf_flag or install_trf_large_flag or install_all_flag):
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

    # Install standard TRF (download pre-compiled)
    if install_trf_flag:
        logger.info("Installing standard TRF (pre-compiled binary)...")
        if install_trf_standard(force=True):
            logger.info("✓ Standard TRF installed successfully")
        else:
            logger.error("✗ Standard TRF installation failed")
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
        if taxid is None:
            logger.info(f"No taxid provided, taxon set to '{DEFAULT_TAXON_NAME}'")
            taxon_name = DEFAULT_TAXON_NAME
        else:
            taxon_name = get_taxon_name(taxid)
            if taxon_name is None:
                logger.warning(f"Failed to retrieve taxon name for taxid {taxid} (invalid ID or NCBI connection problem)")
                logger.warning(f"Taxon set to '{DEFAULT_TAXON_NAME}'")
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

    # Extract run mode flags
    # TRF is disabled by default, FasTAN is the default tool
    run_trf = args.get("run_trf", False)  # Note: argparse converts --run-trf to run_trf
    skip_fastan = args.get("nofastan", False)

    # Handle deprecated --notrf flag (warn but ignore since TRF is off by default now)
    if args.get("notrf", False):
        logger.warning("--notrf flag is deprecated. TRF is now disabled by default. Use --run-trf to enable TRF.")

    # Check for existing FasTAN output when --nofastan is used
    use_existing_fastan = False
    if not run_trf and skip_fastan:
        # Check if FasTAN output already exists (now at output_dir level)
        genome_basename = os.path.splitext(os.path.basename(fasta_file))[0]
        if genome_basename.endswith('.gz'):
            genome_basename = os.path.splitext(genome_basename)[0]
        existing_sat_file = os.path.join(output_dir, f"{genome_basename}.sat")

        if os.path.exists(existing_sat_file) and os.path.getsize(existing_sat_file) > 0:
            logger.info(f"Found existing FasTAN output: {existing_sat_file}")
            logger.info("Running downstream steps only (classification, drawing, report)")
            use_existing_fastan = True
            # Set up settings for existing output - all at output_dir level
            settings["trf_prefix"] = os.path.join(output_dir, genome_basename)
            settings["trf_file"] = existing_sat_file
            settings["output_image_dir"] = os.path.join(output_dir, "images")
            settings["distance_file"] = os.path.join(output_dir, "distances.tsv")
            if not os.path.exists(settings["output_image_dir"]):
                os.makedirs(settings["output_image_dir"])
        else:
            logger.error("Cannot skip FasTAN when TRF is not enabled and no existing FasTAN output found.")
            logger.error(f"Expected file: {existing_sat_file}")
            logger.error("Either remove --nofastan flag or add --run-trf flag.")
            sys.exit(1)

    # Step 1: TRF Search (only if --trf flag is set)
    trf_search_result = None
    force_downstream = force_rerun

    if run_trf:
        logger.info(SEPARATOR_LINE)
        logger.info("STEP 1: TRF SEARCH")
        logger.info(SEPARATOR_LINE)
        trf_search_result = run_trf_search(settings, args, force_rerun)

        # If recompute-failed mode was used and TRF was updated, force regeneration of downstream files
        force_downstream = force_rerun or (trf_search_result == "recomputed")
    else:
        logger.info(SEPARATOR_LINE)
        logger.info("STEP 1: TRF SEARCH - SKIPPED (use --trf to enable)")
        logger.info(SEPARATOR_LINE)

    # Step 1b: FasTAN Analysis (unless --nofastan)
    fastan_success = False
    if not skip_fastan:
        logger.info(SEPARATOR_LINE)
        logger.info("STEP 1b: FASTAN ANALYSIS")
        logger.info(SEPARATOR_LINE)
        fastan_success = run_fastan(settings, force_downstream)
        if not fastan_success:
            logger.error("FasTAN analysis failed!")
            logger.error("Please check the error messages above.")
            logger.error("You can try:")
            logger.error("  1. Remove existing binaries: rm -rf ~/.satellome/bin/fastan ~/.satellome/bin/tanbed")
            logger.error("  2. Reinstall: satellome --install-fastan --install-tanbed")
            logger.error("  3. Run again")
            sys.exit(1)

        # Update settings to use FasTAN output for downstream steps
        # Output files are now at output_dir level, not in fastan/
        genome_basename = os.path.splitext(os.path.basename(settings["fasta_file"]))[0]
        if genome_basename.endswith('.gz'):
            genome_basename = os.path.splitext(genome_basename)[0]
        settings["trf_prefix"] = os.path.join(settings["output_dir"], genome_basename)
        settings["trf_file"] = f"{settings['trf_prefix']}.sat"
        settings["output_image_dir"] = os.path.join(settings["output_dir"], "images")
        settings["distance_file"] = os.path.join(settings["output_dir"], "distances.tsv")
        # Create images directory if needed
        if not os.path.exists(settings["output_image_dir"]):
            os.makedirs(settings["output_image_dir"])
        logger.info(f"Using FasTAN output for downstream analysis: {settings['trf_prefix']}")
    else:
        logger.info(SEPARATOR_LINE)
        logger.info("STEP 1b: FASTAN ANALYSIS - SKIPPED (--nofastan flag)")
        logger.info(SEPARATOR_LINE)

    # Remaining steps run if TRF was executed OR FasTAN was successful OR using existing FasTAN output
    if run_trf or fastan_success or use_existing_fastan:
        # Step 2: Add annotations (only if GFF provided)
        if settings.get("gff_file") or settings.get("repeatmasker_file"):
            logger.info(SEPARATOR_LINE)
            logger.info("STEP 2: ADD ANNOTATIONS")
            logger.info(SEPARATOR_LINE)
            add_annotations(settings, force_downstream)
        else:
            logger.info(SEPARATOR_LINE)
            logger.info("STEP 2: ADD ANNOTATIONS - SKIPPED (no GFF/RM files provided)")
            logger.info(SEPARATOR_LINE)

        # Step 3: Classification
        logger.info(SEPARATOR_LINE)
        logger.info("STEP 3: CLASSIFICATION")
        logger.info(SEPARATOR_LINE)
        run_trf_classification(settings, args, force_downstream)

        # Step 4: Drawing and HTML report
        logger.info(SEPARATOR_LINE)
        logger.info("STEP 4: DRAWING AND REPORT")
        logger.info(SEPARATOR_LINE)
        run_trf_drawing(settings, force_downstream)

    # Print summary
    print_summary(project, taxon_name, output_dir, html_report_file)


if __name__ == "__main__":
    main()
