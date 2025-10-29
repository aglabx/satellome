#!/usr/bin/env python3
"""
Run Satellome analysis on multiple genome assemblies in parallel.
Handles .fa.gz files and manages parallel execution with resource limits.
"""

import os
import sys
import argparse
import subprocess
import gzip
import shutil
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import time
from datetime import datetime, timedelta
from tqdm import tqdm
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def decompress_genome(gz_file, temp_dir, min_free_gb=500):
    """Decompress a .fa.gz file to a temporary directory.

    Args:
        gz_file: Path to compressed genome file
        temp_dir: Temporary directory for decompression
        min_free_gb: Minimum required free space in GB
    """
    base_name = os.path.basename(gz_file).replace('.gz', '')
    output_path = os.path.join(temp_dir, base_name)

    # Check available space
    stat = os.statvfs(temp_dir)
    free_gb = (stat.f_bavail * stat.f_frsize) / (1024**3)
    if free_gb < min_free_gb:
        raise Exception(f"Insufficient disk space in {temp_dir}: {free_gb:.1f}GB free (need at least {min_free_gb}GB)")
    
    with gzip.open(gz_file, 'rb') as f_in:
        with open(output_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    
    return output_path


def run_satellome(genome_file, output_base_dir, threads_per_job, temp_dir, **satellome_args):
    """Run Satellome analysis on a single genome file.

    Args:
        genome_file: Path to input genome file
        output_base_dir: Base output directory
        threads_per_job: Number of threads per job
        temp_dir: Temporary directory
        **satellome_args: Additional arguments for satellome
    """
    # Get genome ID from filename (e.g., GCF_963932015.1 from GCF_963932015.1.fa.gz)
    base_name = os.path.basename(genome_file)
    # Remove all known FASTA extensions
    for ext in ['.fa.gz', '.fasta.gz', '.fna.gz', '.fa', '.fasta', '.fna']:
        if base_name.endswith(ext):
            genome_id = base_name[:-len(ext)]
            break
    else:
        # Fallback to splitting on first dot
        genome_id = base_name.split('.')[0]
    
    # Create output directory
    output_dir = os.path.join(output_base_dir, genome_id)
    os.makedirs(output_dir, exist_ok=True)
    
    # Check if already completed
    project_name = satellome_args.get('project', 'vgp')
    report_file = os.path.join(output_dir, 'reports', 'satellome_report.html')
    if os.path.exists(report_file) and not satellome_args.get('force', False):
        return genome_id, True, "Already completed (skipped)", 0, "skipped"
    
    # Decompress if needed
    if genome_file.endswith('.gz'):
        min_free_gb = satellome_args.get('min_free_gb', 500)
        decompressed_file = decompress_genome(genome_file, temp_dir, min_free_gb)
    else:
        decompressed_file = genome_file
    
    # Build command
    # Find satellome main.py relative to this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    main_py_path = os.path.join(os.path.dirname(script_dir), 'src', 'satellome', 'main.py')

    # Build base command
    if not os.path.exists(main_py_path):
        # Try using the installed satellome command
        base_cmd = ['satellome']
    else:
        base_cmd = ['python3', main_py_path]

    # Add required arguments
    cmd = base_cmd + [
        '-i', decompressed_file,
        '-o', output_dir,
        '-p', satellome_args.get('project', 'satellome'),
        '-t', str(threads_per_job)
    ]

    # Add optional arguments
    if 'cutoff' in satellome_args:
        cmd.extend(['-c', str(satellome_args['cutoff'])])
    if 'minimal_scaffold_length' in satellome_args:
        cmd.extend(['-l', str(satellome_args['minimal_scaffold_length'])])
    if 'large_file' in satellome_args:
        cmd.extend(['--large_file', satellome_args['large_file']])
    if 'drawing_enhancing' in satellome_args:
        cmd.extend(['-e', str(satellome_args['drawing_enhancing'])])
    if 'taxid' in satellome_args:
        cmd.extend(['--taxid', str(satellome_args['taxid'])])
    if 'taxon' in satellome_args:
        cmd.extend(['--taxon', satellome_args['taxon']])
    if satellome_args.get('force', False):
        cmd.append('--force')
    if satellome_args.get('use_kmer_filter', False):
        cmd.append('--use_kmer_filter')
    if 'kmer_threshold' in satellome_args:
        cmd.extend(['--kmer_threshold', str(satellome_args['kmer_threshold'])])
    if 'kmer_bed' in satellome_args:
        cmd.extend(['--kmer_bed', satellome_args['kmer_bed']])
    if satellome_args.get('continue_on_error', False):
        cmd.append('--continue-on-error')
    if satellome_args.get('keep_trf', False):
        cmd.append('--keep-trf')
    
    try:
        # Run the command
        start_time = time.time()
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        elapsed_time = time.time() - start_time
        
        # Clean up decompressed file if we created one
        if genome_file.endswith('.gz') and os.path.exists(decompressed_file):
            os.remove(decompressed_file)
        
        return genome_id, True, f"Success in {elapsed_time:.1f}s", elapsed_time, "completed"
        
    except subprocess.CalledProcessError as e:
        # Clean up decompressed file even if failed
        if genome_file.endswith('.gz') and os.path.exists(decompressed_file):
            os.remove(decompressed_file)
            
        error_msg = f"Failed with return code {e.returncode}"
        return genome_id, False, error_msg, 0, "failed"
    except Exception as e:
        # Clean up decompressed file even if failed
        if genome_file.endswith('.gz') and os.path.exists(decompressed_file):
            os.remove(decompressed_file)
            
        return genome_id, False, str(e), 0, "failed"


def main():
    parser = argparse.ArgumentParser(description='Run Satellome on multiple genomes in parallel',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Parallel execution arguments
    parallel_group = parser.add_argument_group('Parallel execution options')
    parallel_group.add_argument('-i', '--input_dir', required=True, help='Directory containing .fa.gz files')
    parallel_group.add_argument('-o', '--output_dir', required=True, help='Base output directory')
    parallel_group.add_argument('-j', '--jobs', type=int, default=10, help='Number of parallel jobs')
    parallel_group.add_argument('-t', '--threads', type=int, default=10, help='Number of threads per job')
    parallel_group.add_argument('--temp_dir', default='/tmp/satellome_temp', help='Temporary directory for decompressed files')
    parallel_group.add_argument('--min_free_gb', type=int, default=500, help='Minimum free disk space required (GB)')
    parallel_group.add_argument('--dry_run', action='store_true', help='Show what would be run without executing')

    # Satellome parameters
    satellome_group = parser.add_argument_group('Satellome analysis options')
    satellome_group.add_argument('-p', '--project', default='satellome', help='Project name')
    satellome_group.add_argument('-c', '--cutoff', type=int, default=1000, help='Cutoff for large TRs')
    satellome_group.add_argument('-l', '--minimal_scaffold_length', type=int, default=10000, help='Minimal scaffold length')
    satellome_group.add_argument('-e', '--drawing_enhancing', type=int, default=100000, help='Drawing enhancing')
    satellome_group.add_argument('--large_file', default='1kb', help='Suffix for TR file for analysis (e.g., 1kb, 3kb, 10kb)')
    satellome_group.add_argument('--taxid', type=str, help='NCBI taxid')
    satellome_group.add_argument('--taxon', type=str, help='Taxon name')
    satellome_group.add_argument('--force', action='store_true', help='Force rerun all steps even if output files exist')

    # K-mer filtering options
    kmer_group = parser.add_argument_group('K-mer filtering options')
    kmer_group.add_argument('--use_kmer_filter', action='store_true', help='Use k-mer profiling to filter repeat-poor regions')
    kmer_group.add_argument('--kmer_threshold', type=int, default=90000, help='Unique k-mer threshold for repeat detection')
    kmer_group.add_argument('--kmer_bed', type=str, help='Pre-computed k-mer profile BED file from varprofiler')

    # Error handling options
    error_group = parser.add_argument_group('Error handling options')
    error_group.add_argument('--continue_on_error', action='store_true', help='Continue pipeline even if some TRF runs fail')
    error_group.add_argument('--keep_trf', action='store_true', help='Keep original TRF files before filtering')
    
    args = parser.parse_args()
    
    # Validate input directory
    if not os.path.isdir(args.input_dir):
        logger.error(f"Input directory {args.input_dir} does not exist")
        sys.exit(1)
    
    # Create output and temp directories
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(args.temp_dir, exist_ok=True)
    
    # Find all FASTA files with various extensions
    genome_files = []
    fasta_patterns = ['*.fa.gz', '*.fasta.gz', '*.fna.gz']
    
    for pattern in fasta_patterns:
        for file in Path(args.input_dir).glob(pattern):
            # Skip RepeatModeler files
            if 'repeatModeler' in file.name or 'RepeatModeler' in file.name:
                continue
            genome_files.append(str(file))
    
    # Remove duplicates and sort
    genome_files = sorted(list(set(genome_files)))
    
    if not genome_files:
        logger.error(f"No FASTA files found in {args.input_dir}")
        logger.error(f"Looking for: {', '.join(fasta_patterns)}")
        sys.exit(1)
    
    logger.info(f"Found {len(genome_files)} genome files")
    logger.info(f"Will run {args.jobs} parallel jobs with {args.threads} threads each")
    logger.info(f"Total CPU usage: {args.jobs * args.threads} threads")
    
    if args.dry_run:
        logger.info("\nDry run - files that would be processed:")
        for f in sorted(genome_files):
            logger.info(f"  - {os.path.basename(f)}")
        return
    
    # Run analyses in parallel
    logger.info(f"\nStarting parallel analysis...")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Temporary directory: {args.temp_dir}\n")
    
    results = []
    completed = 0
    failed_count = 0
    total_time = 0
    start_time = time.time()
    
    # Create progress bars
    with tqdm(total=len(genome_files), desc="Overall Progress", position=0, leave=True, 
              bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}] {postfix}') as pbar_main:
        with tqdm(desc="Current Status", position=1, bar_format='{desc}', leave=False) as pbar_status:
            with tqdm(desc=f"Active Jobs (0/{args.jobs})", position=2, bar_format='{desc}', leave=False) as pbar_active:
                
                with ProcessPoolExecutor(max_workers=args.jobs) as executor:
                    # Submit all jobs
                    # Prepare satellome arguments
                    satellome_args = {
                        'project': args.project,
                        'cutoff': args.cutoff,
                        'minimal_scaffold_length': args.minimal_scaffold_length,
                        'drawing_enhancing': args.drawing_enhancing,
                        'large_file': args.large_file,
                        'force': args.force,
                        'use_kmer_filter': args.use_kmer_filter,
                        'kmer_threshold': args.kmer_threshold,
                        'continue_on_error': args.continue_on_error,
                        'keep_trf': args.keep_trf,
                        'min_free_gb': args.min_free_gb
                    }

                    # Add optional arguments
                    if args.taxid:
                        satellome_args['taxid'] = args.taxid
                    if args.taxon:
                        satellome_args['taxon'] = args.taxon
                    if args.kmer_bed:
                        satellome_args['kmer_bed'] = args.kmer_bed

                    future_to_genome = {}
                    for genome_file in genome_files:
                        future = executor.submit(run_satellome, genome_file, args.output_dir,
                                               args.threads, args.temp_dir, **satellome_args)
                        future_to_genome[future] = genome_file
                    
                    # Process completed jobs
                    skipped_count = 0
                    for future in as_completed(future_to_genome):
                        genome_file = future_to_genome[future]
                        genome_id, success, message, elapsed, status = future.result()
                        results.append((genome_id, success, message))
                        
                        completed += 1
                        
                        if status == "skipped":
                            skipped_count += 1
                            pbar_status.set_description(f"⏭️  Skipped: {genome_id[:30]} (already completed)")
                        elif success:
                            if elapsed > 0:  # Only count actual runs for timing
                                total_time += elapsed
                            avg_time = total_time / max(1, (completed - failed_count - skipped_count))
                            remaining = len(genome_files) - completed
                            eta = avg_time * remaining if avg_time > 0 else 0
                            
                            pbar_main.set_postfix({
                                'Success': completed - failed_count - skipped_count,
                                'Failed': failed_count,
                                'Skipped': skipped_count,
                                'Avg': f'{avg_time:.0f}s',
                                'ETA': str(timedelta(seconds=int(eta)))
                            })
                            pbar_status.set_description(f"✓ Completed: {genome_id[:30]} ({elapsed:.0f}s)")
                        else:
                            failed_count += 1
                            avg_time = total_time / max(1, (completed - failed_count - skipped_count))
                            pbar_main.set_postfix({
                                'Success': completed - failed_count - skipped_count,
                                'Failed': failed_count,
                                'Skipped': skipped_count,
                                'Avg': f'{avg_time:.0f}s' if avg_time > 0 else 'N/A'
                            })
                            pbar_status.set_description(f"✗ Failed: {genome_id[:30]}")
                        
                        pbar_main.update(1)
                        
                        # Update active jobs counter
                        active_jobs = len([f for f in future_to_genome if not f.done()])
                        pbar_active.set_description(f"Active Jobs ({active_jobs}/{args.jobs})")
    
    # Summary
    total_elapsed = time.time() - start_time
    logger.info("\n" + "="*70)
    logger.info("SUMMARY")
    logger.info("="*70)
    
    successful = [r for r in results if r[1]]
    failed = [r for r in results if not r[1]]
    
    logger.info(f"Total genomes processed: {len(results)}")
    logger.info(f"Successful: {len(successful)} ({len(successful)/len(results)*100:.1f}%)")
    logger.info(f"Failed: {len(failed)} ({len(failed)/len(results)*100:.1f}%)")
    logger.info(f"Total time: {str(timedelta(seconds=int(total_elapsed)))}")
    logger.info(f"Average time per genome: {total_time/len(successful):.1f}s" if successful else "N/A")
    
    if failed:
        logger.warning("\nFailed genomes:")
        for genome_id, _, error in failed:
            logger.warning(f"  - {genome_id}: {error}")

        # Write failed list to file
        failed_file = os.path.join(args.output_dir, "failed_genomes.txt")
        with open(failed_file, 'w') as f:
            for genome_id, _, error in failed:
                f.write(f"{genome_id}\t{error}\n")
        logger.info(f"\nFailed genomes list saved to: {failed_file}")
    
    # Write success list to file
    if successful:
        success_file = os.path.join(args.output_dir, "successful_genomes.txt")
        with open(success_file, 'w') as f:
            for genome_id, _, message in successful:
                f.write(f"{genome_id}\t{message}\n")
        logger.info(f"Successful genomes list saved to: {success_file}")
    
    # Clean up temp directory
    try:
        shutil.rmtree(args.temp_dir)
    except (OSError, PermissionError) as e:
        logger.warning(f"Could not remove temp directory {args.temp_dir}: {e}")
    
    logger.info("\nAnalysis complete!")


if __name__ == '__main__':
    main()