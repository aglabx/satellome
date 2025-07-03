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


def decompress_genome(gz_file, temp_dir):
    """Decompress a .fa.gz file to a temporary directory."""
    base_name = os.path.basename(gz_file).replace('.gz', '')
    output_path = os.path.join(temp_dir, base_name)
    
    # Check available space (require at least 500GB free)
    stat = os.statvfs(temp_dir)
    free_gb = (stat.f_bavail * stat.f_frsize) / (1024**3)
    if free_gb < 500:
        raise Exception(f"Insufficient disk space in {temp_dir}: {free_gb:.1f}GB free (need at least 500GB)")
    
    with gzip.open(gz_file, 'rb') as f_in:
        with open(output_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    
    return output_path


def run_satellome(genome_file, output_base_dir, threads_per_job, temp_dir):
    """Run Satellome analysis on a single genome file."""
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
    report_file = os.path.join(output_dir, 'vgp_report.html')
    if os.path.exists(report_file):
        return genome_id, True, "Already completed (skipped)", 0, "skipped"
    
    # Decompress if needed
    if genome_file.endswith('.gz'):
        decompressed_file = decompress_genome(genome_file, temp_dir)
    else:
        decompressed_file = genome_file
    
    # Build command
    cmd = [
        'python3',
        os.path.expanduser('/home/akomissarov/Dropbox/workspace/PyBioSnippets/satellome/src/satellome/main.py'),
        '-i', decompressed_file,
        '-o', output_dir,
        '-p', 'vgp',
        '-t', str(threads_per_job),
        '-c', '10000',
        '-l', '1000000',
        '--large_file', '10kb'
    ]
    
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
    parser = argparse.ArgumentParser(description='Run Satellome on multiple genomes in parallel')
    parser.add_argument('-i', '--input_dir', required=True, help='Directory containing .fa.gz files')
    parser.add_argument('-o', '--output_dir', required=True, help='Base output directory')
    parser.add_argument('-j', '--jobs', type=int, default=10, help='Number of parallel jobs (default: 10)')
    parser.add_argument('-t', '--threads', type=int, default=10, help='Number of threads per job (default: 10)')
    parser.add_argument('--temp_dir', default='/tmp/satellome_temp', help='Temporary directory for decompressed files')
    parser.add_argument('--dry_run', action='store_true', help='Show what would be run without executing')
    
    args = parser.parse_args()
    
    # Validate input directory
    if not os.path.isdir(args.input_dir):
        print(f"Error: Input directory {args.input_dir} does not exist")
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
        print(f"No FASTA files found in {args.input_dir}")
        print(f"Looking for: {', '.join(fasta_patterns)}")
        sys.exit(1)
    
    print(f"Found {len(genome_files)} genome files")
    print(f"Will run {args.jobs} parallel jobs with {args.threads} threads each")
    print(f"Total CPU usage: {args.jobs * args.threads} threads")
    
    if args.dry_run:
        print("\nDry run - files that would be processed:")
        for f in sorted(genome_files):
            print(f"  - {os.path.basename(f)}")
        return
    
    # Run analyses in parallel
    print(f"\nStarting parallel analysis...")
    print(f"Output directory: {args.output_dir}")
    print(f"Temporary directory: {args.temp_dir}\n")
    
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
                    future_to_genome = {}
                    for genome_file in genome_files:
                        future = executor.submit(run_satellome, genome_file, args.output_dir, 
                                               args.threads, args.temp_dir)
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
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    
    successful = [r for r in results if r[1]]
    failed = [r for r in results if not r[1]]
    
    print(f"Total genomes processed: {len(results)}")
    print(f"Successful: {len(successful)} ({len(successful)/len(results)*100:.1f}%)")
    print(f"Failed: {len(failed)} ({len(failed)/len(results)*100:.1f}%)")
    print(f"Total time: {str(timedelta(seconds=int(total_elapsed)))}")
    print(f"Average time per genome: {total_time/len(successful):.1f}s" if successful else "N/A")
    
    if failed:
        print("\nFailed genomes:")
        for genome_id, _, error in failed:
            print(f"  - {genome_id}: {error}")
        
        # Write failed list to file
        failed_file = os.path.join(args.output_dir, "failed_genomes.txt")
        with open(failed_file, 'w') as f:
            for genome_id, _, error in failed:
                f.write(f"{genome_id}\t{error}\n")
        print(f"\nFailed genomes list saved to: {failed_file}")
    
    # Write success list to file
    if successful:
        success_file = os.path.join(args.output_dir, "successful_genomes.txt")
        with open(success_file, 'w') as f:
            for genome_id, _, message in successful:
                f.write(f"{genome_id}\t{message}\n")
        print(f"Successful genomes list saved to: {success_file}")
    
    # Clean up temp directory
    try:
        shutil.rmtree(args.temp_dir)
    except:
        pass
    
    print("\nAnalysis complete!")


if __name__ == '__main__':
    main()