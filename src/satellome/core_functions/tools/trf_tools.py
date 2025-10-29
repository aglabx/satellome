#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2007-2009 Aleksey Komissarov ( ad3002@gmail.com )
# All rights reserved.
#
# This software is licensed as described in the file COPYING, which
# you should have received as part of this distribution.
"""
TRF search wrapper

- trf_search(file_name="")
- trf_search_in_dir(folder, verbose=False, file_suffix=".fa")

Command example: **wgs.AADD.1.gbff.fa 2 5 7 80 10 50 2000 -m -f -d -h**
"""

import logging
import os
import shutil
import tempfile
import subprocess
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor
import sys
import glob
import concurrent.futures

logger = logging.getLogger(__name__)


from satellome.core_functions.io.fasta_file import sc_iter_fasta_brute
from satellome.core_functions.io.file_system import iter_filepath_folder
from satellome.core_functions.io.trf_file import TRFFileIO
from satellome.core_functions.tools.processing import get_genome_size
from satellome.constants import (
    TRF_DEFAULT_PARAMS, TRF_FLAGS,
    KMER_THRESHOLD_DEFAULT,
    TR_CUTOFF_LARGE,
    MIN_SCAFFOLD_LENGTH_FILTER
)

trf_reader = TRFFileIO().iter_parse

def restore_coordinates_in_line(trf_line):
    """
    Restore original coordinates from chunk coordinates in TRF output line.
    
    Handles headers like: chr10__127750000_127925000 80932 91887
    Should restore to: chr10 208682932 208693887
    
    Args:
        trf_line: TRF output line with modified header
    
    Returns:
        TRF line with restored coordinates
    """
    if not trf_line.strip():
        return trf_line
        
    parts = trf_line.strip().split('\t')
    if len(parts) < 3:
        return trf_line
    
    # Check if header contains coordinate info
    header = parts[0]
    if '__' in header:
        # Split the header to get base name and coordinates
        # Format: chr10__127750000_127925000
        base_parts = header.split('__')
        if len(base_parts) == 2:
            base_header = base_parts[0]
            coord_info = base_parts[1]
            
            # Parse the chunk coordinates
            if '_' in coord_info:
                coords = coord_info.split('_')
                if len(coords) >= 2:
                    try:
                        chunk_start = int(coords[0])
                        # Note: coords[1] is the chunk_end, which we don't need
                        
                        # Adjust TRF coordinates
                        trf_start = int(parts[1])
                        trf_end = int(parts[2])
                        
                        # Update the parts with restored coordinates
                        parts[0] = base_header
                        parts[1] = str(trf_start + chunk_start)
                        parts[2] = str(trf_end + chunk_start)
                        
                        return '\t'.join(parts) + '\n'
                    except (ValueError, IndexError) as e:
                        # If parsing fails, log warning and return original line
                        logger.warning(f"Failed to parse chunk coordinates in TRF line: {e}. Returning original line.")
                        return trf_line

    return trf_line

def run_trf(trf_path, fa_file, max_retries=3):
    """Run TRF on a single file with retry logic.

    Args:
        trf_path: Path to TRF binary
        fa_file: Input FASTA file
        max_retries: Maximum number of retry attempts

    Returns:
        True if successful, raises exception on failure
    """
    # Ensure we use absolute paths
    fa_file_abs = os.path.abspath(fa_file)

    command = [trf_path, fa_file_abs] + TRF_DEFAULT_PARAMS + TRF_FLAGS

    last_stdout = ""
    last_stderr = ""
    last_returncode = None

    for attempt in range(max_retries):
        try:
            # Run TRF - NOTE: TRF versions before 4.10.0 return non-zero on success!
            # They return the number of tandem repeats found, not 0
            # Capture stdout and stderr to provide informative error messages
            result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                  check=False, text=True)

            last_stdout = result.stdout
            last_stderr = result.stderr
            last_returncode = result.returncode

            # Check if TRF actually failed by looking for the output .dat file
            trf_params_suffix = ".".join(TRF_DEFAULT_PARAMS)
            expected_dat_file = fa_file_abs + f".{trf_params_suffix}.dat"
            if os.path.exists(expected_dat_file):
                # Success - .dat file was created
                return True

            # Check if TRF was killed by signal
            if result.returncode < 0:
                # Negative return code usually means killed by signal
                error_msg = f"TRF was killed (signal {-result.returncode}) for file {fa_file}"
                if attempt < max_retries - 1:
                    logger.warning(f"{error_msg}, retrying (attempt {attempt + 2}/{max_retries})...")
                    # Small delay before retry
                    import time
                    time.sleep(2)
                    continue
                else:
                    logger.error(f"{error_msg} after {max_retries} attempts")
                    raise subprocess.CalledProcessError(result.returncode, command)

        except FileNotFoundError:
            logger.error(f"TRF binary not found at: {trf_path}")
            raise

    # Build detailed error message with actual TRF output
    error_details = [
        f"TRF failed for {fa_file} after {max_retries} attempts",
        f"Expected output file: {expected_dat_file}",
        f"Last return code: {last_returncode}",
    ]

    if last_stderr.strip():
        error_details.append(f"TRF stderr: {last_stderr.strip()}")

    if last_stdout.strip():
        error_details.append(f"TRF stdout: {last_stdout.strip()}")

    error_details.append(f"Command: {' '.join(command)}")

    error_message = "\n".join(error_details)
    logger.error(error_message)
    raise RuntimeError(error_message)

def trf_search_by_splitting(
    fasta_file,
    threads=30,
    wdir=".",
    project="NaN",
    trf_path="trf",
    parser_program="./trf_parse_raw.py",
    keep_raw=False,
    genome_size=None,
    use_kmer_filter=False,
    kmer_threshold=KMER_THRESHOLD_DEFAULT,
    kmer_bed_file=None,
    abort_on_error=True,
):
    """TRF search by splitting on fasta file in files."""
    folder_path = tempfile.mkdtemp(dir=wdir)

    if genome_size is None:
        genome_size = get_genome_size(fasta_file)

    # Initialize fa_files list
    fa_files = []
    used_smart_splitting = False
    
    # Check if we should use smart k-mer based splitting
    if use_kmer_filter or kmer_bed_file:
        try:
            from satellome.core_functions.tools.kmer_splitting import split_genome_smart
            logger.info("Using k-mer based smart splitting...")
            output_files = split_genome_smart(
                fasta_file,
                folder_path,
                project,
                threads=int(threads),
                kmer_threshold=kmer_threshold,
                use_kmer_filter=use_kmer_filter,
                kmer_bed_file=kmer_bed_file
            )
            # Keep full paths for now, will convert to basenames after changing directory
            fa_files = output_files
            used_smart_splitting = True
        except ImportError:
            logger.warning("kmer_splitting module not available, falling back to standard splitting")
            use_kmer_filter = False
        except (OSError, subprocess.CalledProcessError, ValueError, RuntimeError) as e:
            logger.warning(f"k-mer splitting failed ({type(e).__name__}: {e}), falling back to standard splitting")
            use_kmer_filter = False
    
    # Fall back to standard splitting if k-mer filtering not used or failed
    if not used_smart_splitting:
        ### 1. Split chromosomes into temp file
        total_length = 0
        next_file = 0
        
        logger.info("Splitting genome into ~100kb chunks...")
        with tqdm(total=genome_size, desc="Splitting fasta file", unit=" bp", unit_scale=True, unit_divisor=1000, dynamic_ncols=True) as pbar:
            for i, (header, seq) in enumerate(sc_iter_fasta_brute(fasta_file)):
                file_path = os.path.join(folder_path, "%s.fa" % next_file)
                with open(file_path, "a") as fw:
                    fw.write("%s\n%s\n" % (header, seq))
                total_length += len(seq)
                if total_length > 100000:
                    next_file += 1
                    total_length = 0
                pbar.update(len(seq))
        
        # Get list of created files for processing
        fa_files = [f for f in os.listdir(folder_path) if f.endswith('.fa')]
        logger.info(f"Created {len(fa_files)} chunks")

    ### 2. Run TRF
    fasta_name = ".".join(fasta_file.split("/")[-1].split(".")[:-1])
    output_file = os.path.join(wdir, fasta_name + ".trf")

    current_dir = os.getcwd()

    os.chdir(folder_path)

    # Convert full paths to basenames after changing directory
    if used_smart_splitting and fa_files:
        fa_files = [os.path.basename(f) for f in fa_files]
    
    # If fa_files is empty (no smart splitting), get list of .fa files
    if not fa_files:
        fa_files = [f for f in os.listdir(folder_path) if f.endswith('.fa')]

    # Track failed files
    failed_files = []
    successful_files = 0
    
    # Create a progress bar
    with tqdm(total=len(fa_files), desc="Running TRF", dynamic_ncols=True) as pbar:
        # Using a thread pool to run TRF processes in parallel
        with ThreadPoolExecutor(max_workers=int(threads)) as executor:
            # Submit all tasks and collect futures
            futures = {}
            for fa_file in fa_files:
                future = executor.submit(run_trf, trf_path, fa_file)
                futures[future] = fa_file
            
            # Process completed futures
            from concurrent.futures import as_completed
            for future in as_completed(futures):
                fa_file = futures[future]
                try:
                    # Get result (will raise exception if TRF failed)
                    result = future.result()
                    successful_files += 1
                except Exception as e:
                    logger.error(f"TRF failed for {fa_file}: {e}")
                    failed_files.append(fa_file)
                finally:
                    pbar.update()
    
    # Check if any files failed
    if failed_files:
        logger.error(f"TRF failed for {len(failed_files)} out of {len(fa_files)} files:")
        for f in failed_files[:10]:  # Show first 10 failed files
            logger.error(f"  - {f}")
        if len(failed_files) > 10:
            logger.error(f"  ... and {len(failed_files) - 10} more")
        
        if abort_on_error:
            # Abort the pipeline
            logger.error("Aborting pipeline due to TRF failures. Please investigate the errors and try again.")
            logger.info("Tip: You can try reducing the number of threads or increasing system resources.")
            os.chdir(current_dir)
            if not keep_raw:
                shutil.rmtree(folder_path)
            raise RuntimeError(f"TRF failed for {len(failed_files)} files. Pipeline aborted.")
        else:
            # Continue with partial results
            logger.warning(f"Continuing with partial results ({successful_files}/{len(fa_files)} files processed)")
            logger.warning("The analysis may be incomplete!")
    
    logger.info(f"TRF completed successfully for all {successful_files} files")

    os.chdir(current_dir)

    ### 3. Parse TRF
    
    # Check if TRF produced any .dat files
    dat_files = [f for f in os.listdir(folder_path) if f.endswith('.dat')]
    
    if len(dat_files) == 0:
        logger.warning("No .dat files found! TRF may have failed to run properly.")

    # Use the current Python interpreter to execute the parser script to avoid permission issues
    python_exe = sys.executable

    # Process .dat files in parallel using subprocess for better security and error handling
    
    dat_files = glob.glob(os.path.join(folder_path, "*.dat"))

    def process_dat_file(dat_file):
        """Process a single .dat file."""
        output_file_path = f"{dat_file}.trf"
        cmd = [
            python_exe,
            parser_program,
            "-i", dat_file,
            "-o", output_file_path,
            "-p", project
        ]
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            return dat_file, True, result.stdout
        except subprocess.CalledProcessError as e:
            return dat_file, False, f"Error: {e.stderr}"

    # Process files in parallel
    with concurrent.futures.ThreadPoolExecutor(max_workers=int(threads)) as executor:
        futures = [executor.submit(process_dat_file, dat_file) for dat_file in dat_files]
        for future in concurrent.futures.as_completed(futures):
            dat_file, success, output = future.result()
            if not success:
                logger.warning(f"Failed to process {dat_file}: {output}")

    ### 3. Aggregate data

    with open(output_file, "w") as fw:
        for file_path in iter_filepath_folder(folder_path):
            if file_path.endswith(".trf"):
                with open(file_path) as fh:
                    # If using smart splitting, restore original coordinates
                    if used_smart_splitting:
                        for line in fh:
                            fw.write(restore_coordinates_in_line(line))
                    else:
                        fw.write(fh.read())

    os.chdir(current_dir)

    ## 4. Remove temp folder
    if not keep_raw:
        if folder_path.count("/") <= 3:
            input("Remove: %s ?" % folder_path)
        shutil.rmtree(folder_path)

    return output_file


def _filter_by_bottom_array_length(obj, cutoff):
    if obj.trf_array_length > cutoff:
        return True
    else:
        return False


def _filter_by_bottom_unit_length(obj, cutoff):
    if obj.trf_period > cutoff:
        return True
    else:
        return False


def trf_filter_by_array_length(trf_file, output_file, cutoff):
    """Create output TRF file with tandem repeats with length greater than from input file.
    Function returns number of tandem repeats in output file.
    """
    i = 0
    with open(output_file, "w") as fw:
        for obj in trf_reader(trf_file):
            if _filter_by_bottom_array_length(obj, cutoff):
                i += 1
                fw.write(obj.get_string_repr())
    logger.info(f"Filtered {i} tandem repeats by array length")
    return i


def trf_filter_by_monomer_length(trf_file, output_file, cutoff):
    """Create output TRF file with tandem repeats with unit length greater than from input file.
    Function returns number of tandem repeats in output file.
    """
    i = 0
    with open(output_file, "w") as fw:
        for obj in trf_reader(trf_file):
            if _filter_by_bottom_unit_length(obj, cutoff):
                i += 1
                fw.write(obj.get_string_repr())
    logger.info(f"Filtered {i} tandem repeats by monomer length")
    return i


def trf_filter_exclude_by_gi_list(trf_file, output_file, gi_list_to_exclude):
    """Create output TRF file with tandem repeats with GI that don't match GI_LIST
    List of GI, see TRF and FA specifications, GI is first value in TRF row.
    """
    with open(output_file, "w") as fw:
        for obj in trf_reader(trf_file):
            if not obj.trf_gi in gi_list_to_exclude:
                fw.write(obj.get_string_repr())


def trf_representation(trf_file, trf_output, representation):
    """Write TRF file tab delimited representation.
    representation: numerical|index|agc_apm|with_monomer|family
    """
    with open(trf_output, "w") as fw:
        for obj in trf_reader(trf_file):
            if representation == "numerical":
                line = obj.get_numerical_repr()
            elif representation == "index":
                line = obj.get_index_repr()
            elif representation == "agc_apm":
                line = "%.2f\t%.2f\n" % (obj.trf_array_gc, obj.trf_pmatch)
            elif representation == "with_monomer":
                line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                    obj.trf_id,
                    obj.trf_period,
                    obj.trf_array_length,
                    obj.trf_array_gc,
                    obj.trf_pmatch,
                    obj.trf_gi,
                    obj.trf_consensus,
                )
            elif representation == "family":
                line = obj.get_family_repr()

            fw.write(line)


def trf_write_field_n_data(trf_file, file_output, field, field_format="%s"):
    """Write statistics data: field, N."""
    result = {}
    with open(file_output, "w") as fw:
        for obj in trf_reader(trf_file):
            value = field_format % getattr(obj, field)
            result.setdefault(value)
            result[value] += 1
        result = [(value, n) for value, n in result.items()]
        result.sort()
        for value, n in result:
            line = field_format + "\t%s\n"
            line = line % (value, n)
            fw.write(line)


def trf_write_two_field_data(trf_file, file_output, field_a, field_b):
    """Write statistics data: field_a, field_b."""
    result = []
    with open(file_output, "w") as fw:
        for obj in trf_reader(trf_file):
            value_a = getattr(obj, field_a)
            value_b = getattr(obj, field_b)
            result.append([value_a, value_b])
        result.sort()
        for value_a, value_b in result:
            line = "%s\t%s\n" % (value_a, value_b)
            fw.write(line)


def count_trs_per_chrs(all_trf_file):
    """Function prints chr, all trs, 3000 trs, 10000 trs"""
    chr2n = {}
    chr2n_large = {}
    chr2n_xlarge = {}
    for trf_obj in trf_reader(all_trf_file):
        chr = trf_obj.trf_chr
        chr2n.setdefault(chr, 0)
        chr2n_large.setdefault(chr, 0)
        chr2n_xlarge.setdefault(chr, 0)
        chr2n[chr] += 1
        if trf_obj.trf_array_length > 3000:
            chr2n_large[chr] += 1
        if trf_obj.trf_array_length > TR_CUTOFF_LARGE:
            chr2n_xlarge[chr] += 1
    for chr in chr2n:
        logger.info(f"{chr} {chr2n[chr]} {chr2n_large[chr]} {chr2n_xlarge[chr]}")


def count_trf_subset_by_head(trf_file, head_value):
    """Function prints number of items with given fasta head fragment"""
    total = 0
    n = 0
    total_length = 0
    for trf_obj in trf_reader(trf_file):
        total += 1
        if head_value in trf_obj.trf_head:
            n += 1
            total_length += trf_obj.trf_array_length
    return n, total, total_length


def fix_chr_names(trf_file, temp_file_name=None, case=None):
    """Some fasta heads impossible to parse, so it is simpler to fix them postfactum

    Cases:

    - chromosome names in MSGC genome assembly [MSGC sequences]

    """

    if not temp_file_name:
        temp_file_name = trf_file + ".tmp"
    if case == "MSGC sequences":
        with open(temp_file_name, "a") as fw:
            for obj in trf_reader(trf_file):
                obj.trf_chr = obj.trf_gi
                fw.write(obj.get_string_repr())
    if os.path.isfile(temp_file_name):
        os.remove(trf_file)
        os.rename(temp_file_name, trf_file)


def recompute_failed_chromosomes(
    fasta_file,
    existing_trf_file,
    output_dir,
    project,
    threads=30,
    trf_path="trf",
    parser_program="./trf_parse_raw.py",
    min_scaffold_size=1000000,
    match_first_word=True,
):
    """Recompute TRF only for chromosomes/contigs that are missing or failed.

    This function:
    1. Checks which chromosomes are missing from existing TRF results
    2. Extracts only those chromosomes to a temporary FASTA
    3. Runs TRF only on the missing chromosomes
    4. Merges results back into the existing TRF file

    Args:
        fasta_file: Original FASTA file
        existing_trf_file: Existing TRF file to check and update
        output_dir: Output directory for temporary files
        project: Project name
        threads: Number of threads
        trf_path: Path to TRF binary
        parser_program: Path to TRF parser script
        min_scaffold_size: Minimum scaffold size to check (default 1Mb)
        match_first_word: Match only first word of scaffold names

    Returns:
        True if recomputation was successful, False otherwise
    """
    from collections import defaultdict
    from satellome.core_functions.io.fasta_file import sc_iter_fasta_brute
    from satellome.core_functions.io.tab_file import sc_iter_tab_file
    from satellome.core_functions.models.trf_model import TRModel

    logger.info("="*60)
    logger.info("SMART RECOMPUTE MODE: Checking for failed chromosomes...")
    logger.info("="*60)

    # Step 1: Get scaffold lengths from FASTA
    logger.info(f"Reading scaffold lengths from {fasta_file}...")
    scaffold_lengths = {}
    scaffold_sequences = {}

    for header, sequence in sc_iter_fasta_brute(fasta_file):
        scaffold_name = header.replace(">", "").strip()
        original_name = scaffold_name

        if match_first_word:
            scaffold_name = scaffold_name.split()[0] if scaffold_name else scaffold_name

        scaffold_lengths[scaffold_name] = len(sequence)
        scaffold_sequences[scaffold_name] = (original_name, sequence)

    logger.info(f"Found {len(scaffold_lengths):,} scaffolds in FASTA")

    # Step 2: Get scaffolds that have TRF results
    logger.info(f"Reading existing TRF results from {existing_trf_file}...")
    scaffold_trf_counts = defaultdict(int)

    if os.path.exists(existing_trf_file) and os.path.getsize(existing_trf_file) > 0:
        for trf_obj in sc_iter_tab_file(existing_trf_file, TRModel):
            scaffold_name = trf_obj.trf_head if hasattr(trf_obj, 'trf_head') else 'Unknown'

            if match_first_word and scaffold_name and scaffold_name != 'Unknown':
                scaffold_name = scaffold_name.split()[0]

            scaffold_trf_counts[scaffold_name] += 1

        logger.info(f"Found {sum(scaffold_trf_counts.values()):,} tandem repeats across {len(scaffold_trf_counts):,} scaffolds in existing TRF")
    else:
        logger.warning(f"Existing TRF file not found or empty: {existing_trf_file}")
        logger.info("Will process all chromosomes...")

    # Step 3: Find missing/failed scaffolds (large scaffolds without TRF results)
    missing_scaffolds = []

    for scaffold_name, length in scaffold_lengths.items():
        if length < min_scaffold_size:
            continue

        tr_count = scaffold_trf_counts.get(scaffold_name, 0)

        if tr_count == 0:
            missing_scaffolds.append(scaffold_name)

    if not missing_scaffolds:
        logger.info("✅ All large scaffolds have TRF results! Nothing to recompute.")
        return True

    logger.warning(f"Found {len(missing_scaffolds)} large scaffold(s) with NO tandem repeats:")
    for scaffold in missing_scaffolds[:10]:
        logger.warning(f"  - {scaffold}: {scaffold_lengths[scaffold]:,} bp")
    if len(missing_scaffolds) > 10:
        logger.warning(f"  ... and {len(missing_scaffolds) - 10} more")

    # Step 4: Create temporary FASTA with only missing scaffolds
    temp_fasta = os.path.join(output_dir, f"{project}_missing_scaffolds.fasta")
    logger.info(f"\nCreating temporary FASTA with missing scaffolds: {temp_fasta}")

    with open(temp_fasta, 'w') as fw:
        for scaffold_name in missing_scaffolds:
            if scaffold_name in scaffold_sequences:
                original_name, sequence = scaffold_sequences[scaffold_name]
                fw.write(f">{original_name}\n")
                fw.write(f"{sequence}\n")

    logger.info(f"Wrote {len(missing_scaffolds)} scaffolds to temporary FASTA")

    # Step 5: Run TRF on missing scaffolds
    logger.info("\nRunning TRF on missing scaffolds...")
    temp_trf_file = os.path.join(output_dir, f"{project}_missing_scaffolds.trf")

    try:
        # Use trf_search_by_splitting for the missing scaffolds
        trf_search_by_splitting(
            fasta_file=temp_fasta,
            threads=threads,
            wdir=output_dir,
            project=f"{project}_missing_scaffolds",
            trf_path=trf_path,
            parser_program=parser_program,
            keep_raw=False,
            abort_on_error=True,  # Fail if TRF fails again
        )

        logger.info("✅ TRF completed successfully for missing scaffolds")

    except Exception as e:
        logger.error(f"❌ TRF failed for missing scaffolds: {e}")
        logger.error("Cannot proceed with recomputation. Please investigate the errors.")
        # Clean up temporary files
        if os.path.exists(temp_fasta):
            os.remove(temp_fasta)
        raise

    # Step 6: Merge results
    if os.path.exists(temp_trf_file) and os.path.getsize(temp_trf_file) > 0:
        logger.info("\nMerging TRF results...")

        # Count new TRs
        new_tr_count = 0
        for _ in sc_iter_tab_file(temp_trf_file, TRModel):
            new_tr_count += 1

        logger.info(f"Found {new_tr_count:,} new tandem repeats")

        # Create backup of original TRF file
        backup_file = f"{existing_trf_file}.before_recompute"
        if os.path.exists(existing_trf_file):
            shutil.copy2(existing_trf_file, backup_file)
            logger.info(f"Created backup: {backup_file}")

        # Append new results to existing TRF file
        with open(existing_trf_file, 'a') as fw:
            for trf_obj in sc_iter_tab_file(temp_trf_file, TRModel):
                fw.write(trf_obj.get_string_repr())

        logger.info(f"✅ Merged results into {existing_trf_file}")

        # Clean up temporary files
        logger.info("\nCleaning up temporary files...")
        if os.path.exists(temp_fasta):
            os.remove(temp_fasta)
        if os.path.exists(temp_trf_file):
            os.remove(temp_trf_file)

        logger.info("="*60)
        logger.info("SMART RECOMPUTE COMPLETED SUCCESSFULLY")
        logger.info("="*60)
        return True

    else:
        logger.error(f"❌ Expected TRF output file not found: {temp_trf_file}")
        return False
