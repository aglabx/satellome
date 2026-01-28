#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# @created: 27.12.2024
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

"""
Batch check TRF consistency for multiple genome assemblies.

This script iterates through directories containing genome assemblies and runs
TRF consistency checks on each one. If issues are found, it prompts the user
to skip or delete the TRF directory.
"""

import os
import sys
import glob
import shutil
import subprocess
from pathlib import Path
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

def find_genome_assemblies(base_dir, include_missing=False):
    """Find all genome assembly directories with TRF results.
    
    Args:
        base_dir: Base directory to search (e.g., 'reptiles')
        include_missing: If True, also return assemblies without TRF results
        
    Returns:
        Tuple of (assemblies_with_trf, assemblies_without_trf)
        Each is a list of tuples (assembly_dir, fasta_file, trf_file)
    """
    assemblies_with_trf = []
    assemblies_without_trf = []
    no_fasta_assemblies = []
    
    # Pattern to find assembly directories (both GCA and GCF prefixes)
    gca_dirs = glob.glob(os.path.join(base_dir, "GCA_*"))
    gcf_dirs = glob.glob(os.path.join(base_dir, "GCF_*"))
    assembly_dirs = sorted(gca_dirs + gcf_dirs)
    
    if not assembly_dirs:
        logger.warning(f"  No GCA_* or GCF_* directories found in {base_dir}")
        return assemblies_with_trf, assemblies_without_trf
    
    logger.info(f"  Found {len(assembly_dirs)} assembly directories ({len(gca_dirs)} GCA, {len(gcf_dirs)} GCF)")
    
    for assembly_dir in sorted(assembly_dirs):
        assembly_name = os.path.basename(assembly_dir)
        
        # Expected file paths
        fasta_file = os.path.join(assembly_dir, f"{assembly_name}_genomic.fna")
        trf_file = os.path.join(assembly_dir, "trf", f"{assembly_name}_genomic.sat")
        
        # Check if both files exist
        if os.path.exists(fasta_file) and os.path.exists(trf_file):
            assemblies_with_trf.append((assembly_dir, fasta_file, trf_file))
        elif os.path.exists(fasta_file) and not os.path.exists(trf_file):
            assemblies_without_trf.append((assembly_dir, fasta_file, None))
            # Check if TRF directory exists but file is missing
            trf_dir = os.path.join(assembly_dir, "trf")
            if os.path.exists(trf_dir):
                logger.warning(f"⚠️  TRF directory exists but no .trf file for {assembly_name}")
                # List contents of TRF directory for debugging
                trf_contents = os.listdir(trf_dir)[:5]  # Show first 5 files
                if trf_contents:
                    logger.warning(f"    TRF directory contains: {trf_contents}")
            else:
                logger.warning(f"⚠️  No TRF directory for {assembly_name}")
        elif not os.path.exists(fasta_file):
            no_fasta_assemblies.append(assembly_name)
    
    # Summary of issues
    if assemblies_without_trf:
        logger.warning(f"\n  ⚠️  {len(assemblies_without_trf)} assemblies without TRF results")
    if no_fasta_assemblies:
        logger.warning(f"  ⚠️  {len(no_fasta_assemblies)} assemblies without FASTA files")
        for name in no_fasta_assemblies[:3]:  # Show first 3
            logger.warning(f"    - {name}")
    
    return assemblies_with_trf, assemblies_without_trf


def run_consistency_check(fasta_file, trf_file, min_size=1000000, match_first_word=True, debug=False):
    """Run TRF consistency check for a single assembly.
    
    Args:
        fasta_file: Path to FASTA file
        trf_file: Path to TRF output file
        min_size: Minimum scaffold size to check
        match_first_word: If True, match only first word of scaffold names
        debug: If True, show debug information
        
    Returns:
        Tuple (success, stdout, stderr)
    """
    script_path = os.path.expanduser("~/Dropbox/workspace/new/biology/satellome/scripts/check_trf_consistency.py")
    
    cmd = [
        "python3",
        script_path,
        "-f", fasta_file,
        "-t", trf_file,
        "-s", str(min_size)
    ]
    
    if not match_first_word:
        cmd.append("--no-match-first-word")
    
    if debug:
        cmd.append("--debug")
    
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True
    )
    return result.returncode == 0, result.stdout, result.stderr
    

def prompt_user_action(assembly_dir, failure_reason=""):
    """Prompt user for action when TRF check fails.
    
    Args:
        assembly_dir: Path to assembly directory
        failure_reason: Explanation of why the check failed
        
    Returns:
        User choice: 's' (skip), 'd' (delete), 'q' (quit)
    """
    trf_dir = os.path.join(assembly_dir, "trf")
    
    while True:
        logger.error("\n" + "="*60)
        logger.error(f"❌ ISSUE DETECTED:")
        logger.error(f"   Genome: {os.path.basename(assembly_dir)}")
        logger.error(f"   Path: {assembly_dir}")

        if failure_reason:
            logger.error(f"\n📊 Reason: {failure_reason}")

        logger.info("\nWhat does this mean?")
        logger.info("  • Some large contigs/chromosomes have no tandem repeats")
        logger.info("  • This is suspicious as genomes usually have repeats")
        logger.info("  • TRF may have failed to process these contigs (error/out of memory)")

        logger.info("\nWhat to do?")
        logger.info("  [s] Skip - continue to next assembly")
        logger.info("  [d] Delete - remove TRF directory and continue")
        logger.info("  [v] View - show TRF directory contents")
        logger.info("  [q] Quit - exit the script")
        
        choice = input("\nYour choice [s/d/v/q]: ").lower().strip()
        
        if choice == 's':
            logger.info("➡️  Skipping this assembly...")
            return 's'
        elif choice == 'd':
            confirm = input(f"⚠️  Are you sure you want to delete {trf_dir}? [y/N]: ").lower().strip()
            if confirm == 'y':
                shutil.rmtree(trf_dir)
                logger.info(f"✅ Deleted {trf_dir}")
                return 'd'
            else:
                logger.info("Deletion cancelled")
        elif choice == 'v':
            # Show TRF directory contents
            if os.path.exists(trf_dir):
                logger.info(f"\nContents of {trf_dir}:")
                for item in os.listdir(trf_dir):
                    item_path = os.path.join(trf_dir, item)
                    size = os.path.getsize(item_path) if os.path.isfile(item_path) else 0
                    logger.info(f"  - {item} ({size:,} bytes)")
            else:
                logger.warning(f"Directory {trf_dir} does not exist")
        elif choice == 'q':
            logger.info("Exiting...")
            return 'q'
        else:
            logger.warning("Invalid choice. Please enter 's', 'd', 'v', or 'q'")


def main():
    """Main function to batch check TRF consistency."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Batch check TRF consistency for multiple genome assemblies"
    )
    parser.add_argument(
        "directories",
        nargs="+",
        help="Base directories containing genome assemblies (e.g., reptiles, mammals)"
    )
    parser.add_argument(
        "-s", "--min-size",
        type=int,
        default=1000000,
        help="Minimum scaffold size to check in bp [1000000]"
    )
    parser.add_argument(
        "-o", "--output",
        help="Output summary file [optional]",
        default=None
    )
    parser.add_argument(
        "--auto-skip",
        action="store_true",
        help="Automatically skip failed assemblies without prompting"
    )
    parser.add_argument(
        "--auto-delete",
        action="store_true",
        help="Automatically delete TRF directories for failed assemblies (DANGEROUS!)"
    )
    parser.add_argument(
        "--no-match-first-word",
        action="store_true",
        help="Use full scaffold names instead of just first word (by default uses first word only)"
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Show debug information about scaffold name matching"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Show detailed output for passed assemblies"
    )
    parser.add_argument(
        "--check-missing",
        action="store_true",
        help="Also check assemblies without TRF results (to identify which need TRF analysis)"
    )
    
    args = parser.parse_args()
    
    # Track statistics
    total_assemblies = 0
    passed_assemblies = []
    failed_assemblies = []
    skipped_assemblies = []
    deleted_assemblies = []
    
    # Process each directory
    for base_dir in args.directories:
        if not os.path.exists(base_dir):
            logger.error(f"❌ Directory not found: {base_dir}")
            continue
        
        logger.info(f"\n{'='*60}")
        logger.info(f"Processing directory: {base_dir}")
        logger.info(f"{'='*60}")
        
        # Find genome assemblies
        assemblies_with_trf, assemblies_without_trf = find_genome_assemblies(base_dir, args.check_missing)
        
        if assemblies_without_trf and args.check_missing:
            logger.info(f"\n  📋 Assemblies needing TRF analysis:")
            for assembly_dir, fasta_file, _ in assemblies_without_trf[:10]:  # Show first 10
                logger.info(f"    - {os.path.basename(assembly_dir)}")
            if len(assemblies_without_trf) > 10:
                logger.info(f"    ... and {len(assemblies_without_trf) - 10} more")
        
        if not assemblies_with_trf:
            logger.error(f"\n  ❌ No genome assemblies with complete TRF results found in {base_dir}")
            logger.error(f"     Check if TRF analysis has been run for this directory")
            continue
        
        logger.info(f"\nFound {len(assemblies_with_trf)} assemblies with TRF results to check")
        
        # Rename variable for compatibility
        assemblies = assemblies_with_trf
        
        # Check each assembly
        for idx, (assembly_dir, fasta_file, trf_file) in enumerate(assemblies, 1):
            total_assemblies += 1
            assembly_name = os.path.basename(assembly_dir)
            
            logger.info(f"\n[{idx}/{len(assemblies)}] Checking {assembly_name}...")
            
            # Run consistency check
            success, stdout, stderr = run_consistency_check(
                fasta_file, 
                trf_file, 
                args.min_size,
                match_first_word=not args.no_match_first_word,
                debug=args.debug
            )
            
            if success:
                logger.info(f"✅ PASSED: {assembly_name}")
                passed_assemblies.append(assembly_name)
                # Optionally show stdout for passed checks
                if stdout and args.verbose:
                    logger.info(stdout)
            else:
                failed_assemblies.append(assembly_name)
                
                # Show the consistency check output (contains the actual report)
                if stdout:
                    logger.warning("\n" + stdout)

                # Show error details if any
                if stderr:
                    logger.error(f"\nError details: {stderr}")
                
                # Parse failure reason from stdout
                failure_reason = ""
                if stdout:
                    if "large scaffold(s) with NO tandem repeats" in stdout:
                        import re
                        match = re.search(r"(\d+) large scaffold\(s\) with NO tandem repeats", stdout)
                        if match:
                            failure_reason = f"Found {match.group(1)} large scaffold(s) with NO tandem repeats"
                    elif "large scaffold(s) with suspiciously few" in stdout:
                        match = re.search(r"(\d+) large scaffold\(s\) with suspiciously few", stdout)
                        if match:
                            failure_reason = f"Found {match.group(1)} large scaffold(s) with suspiciously few tandem repeats"
                
                # Decide action
                if args.auto_skip:
                    logger.info(f"➡️  Auto-skipping {assembly_name}")
                    skipped_assemblies.append(assembly_name)
                elif args.auto_delete:
                    trf_dir = os.path.join(assembly_dir, "trf")
                    try:
                        shutil.rmtree(trf_dir)
                        logger.info(f"🗑️  Auto-deleted TRF directory for {assembly_name}")
                        deleted_assemblies.append(assembly_name)
                    except Exception as e:
                        logger.error(f"❌ Error deleting {trf_dir}: {e}")
                        skipped_assemblies.append(assembly_name)
                else:
                    # Interactive mode
                    action = prompt_user_action(assembly_dir, failure_reason)
                    if action == 'q':
                        break
                    elif action == 's':
                        skipped_assemblies.append(assembly_name)
                    elif action == 'd':
                        deleted_assemblies.append(assembly_name)
    
    # Print summary
    logger.info("\n" + "="*60)
    logger.info("BATCH CHECK SUMMARY")
    logger.info("="*60)
    logger.info(f"Total assemblies checked: {total_assemblies}")
    logger.info(f"✅ Passed: {len(passed_assemblies)}")
    logger.info(f"❌ Failed: {len(failed_assemblies)}")
    logger.info(f"➡️  Skipped: {len(skipped_assemblies)}")
    logger.info(f"🗑️  Deleted: {len(deleted_assemblies)}")
    
    # Save summary to file if requested
    if args.output:
        with open(args.output, 'w') as f:
            f.write("TRF BATCH CONSISTENCY CHECK SUMMARY\n")
            f.write("="*60 + "\n\n")
            f.write(f"Total assemblies checked: {total_assemblies}\n")
            f.write(f"Passed: {len(passed_assemblies)}\n")
            f.write(f"Failed: {len(failed_assemblies)}\n")
            f.write(f"Skipped: {len(skipped_assemblies)}\n")
            f.write(f"Deleted: {len(deleted_assemblies)}\n\n")
            
            if passed_assemblies:
                f.write("PASSED ASSEMBLIES:\n")
                for name in sorted(passed_assemblies):
                    f.write(f"  - {name}\n")
                f.write("\n")
            
            if failed_assemblies:
                f.write("FAILED ASSEMBLIES:\n")
                for name in sorted(failed_assemblies):
                    status = "deleted" if name in deleted_assemblies else "skipped"
                    f.write(f"  - {name} ({status})\n")
                f.write("\n")
        
        logger.info(f"\n📄 Summary saved to: {args.output}")
    
    # Exit with appropriate code
    if failed_assemblies and not (args.auto_skip or args.auto_delete):
        sys.exit(1)
    else:
        sys.exit(0)


if __name__ == "__main__":
    main()