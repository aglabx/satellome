#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# @created: 27.12.2024
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

"""
Check TRF results consistency.

This script verifies that TRF analysis completed successfully for all contigs/scaffolds
above a certain size threshold. It's expected that contigs larger than 1Mb should have
at least some tandem repeats detected.
"""

import argparse
import sys
import os
from collections import defaultdict
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Add parent directory to path to import satellome modules
sys.path.append(os.path.join(os.path.dirname(__file__), "../src"))

from satellome.core_functions.io.fasta_file import sc_iter_fasta_brute
from satellome.core_functions.io.tab_file import sc_iter_tab_file
from satellome.core_functions.models.trf_model import TRModel


def get_scaffold_lengths(fasta_file, match_first_word=True):
    """Get lengths of all scaffolds/contigs from FASTA file.
    
    Args:
        fasta_file: Path to input FASTA file
        match_first_word: If True, use only first word of header (default: True)
        
    Returns:
        Dictionary mapping scaffold name to length
    """
    scaffold_lengths = {}
    
    logger.info(f"Reading scaffold lengths from {fasta_file}...")
    for header, sequence in sc_iter_fasta_brute(fasta_file):
        # Remove '>' and strip whitespace
        scaffold_name = header.replace(">", "").strip()
        
        if match_first_word:
            # Take only first word (like: NC_060925.1 from "NC_060925.1 Siniperca chuatsi chromosome 1")
            scaffold_name = scaffold_name.split()[0] if scaffold_name else scaffold_name
        
        scaffold_lengths[scaffold_name] = len(sequence)
    
    return scaffold_lengths


def get_trf_scaffolds(trf_file, match_first_word=True, debug=False):
    """Get all scaffolds that have TRF results.
    
    Args:
        trf_file: Path to TRF output file
        match_first_word: If True, use only first word of scaffold name (default: True)
        debug: If True, show debug information
        
    Returns:
        Dictionary mapping scaffold name to number of tandem repeats
    """
    scaffold_trf_counts = defaultdict(int)
    total_trs = 0
    debug_first_5 = []
    
    logger.info(f"Reading TRF results from {trf_file}...")
    
    # Use TRModel to parse the file properly
    for trf_obj in sc_iter_tab_file(trf_file, TRModel):
        # trf_head contains the header from FASTA (19th column in TRF file)
        # trf_gi might be parsed from trf_head but could be "Unknown" if parsing fails
        scaffold_name = trf_obj.trf_head if hasattr(trf_obj, 'trf_head') else 'Unknown'
        
        if match_first_word and scaffold_name and scaffold_name != 'Unknown':
            # Take only first word (to match FASTA parsing)
            scaffold_name = scaffold_name.split()[0]
        
        scaffold_trf_counts[scaffold_name] += 1
        total_trs += 1
        
        # Collect first 5 scaffold names for debugging
        if debug and len(debug_first_5) < 5 and scaffold_name not in [x[0] for x in debug_first_5]:
            original_name = trf_obj.trf_head if hasattr(trf_obj, 'trf_head') else 'Unknown'
            debug_first_5.append((scaffold_name, original_name))
    
    logger.info(f"Found {total_trs:,} tandem repeats across {len(scaffold_trf_counts):,} scaffolds")
    
    # Debug: show first few scaffold names from TRF
    if debug and debug_first_5:
        logger.debug(f"Debug - First scaffold names from TRF (trf_head field):")
        for processed, original in debug_first_5[:3]:
            logger.debug(f"  - Processed: '{processed}' (Original trf_head: '{original}')")
    
    return scaffold_trf_counts


def check_consistency(fasta_file, trf_file, min_scaffold_size=1000000, min_expected_trs=1, 
                     match_first_word=True, debug=False):
    """Check consistency between FASTA input and TRF output.
    
    Args:
        fasta_file: Path to input FASTA file
        trf_file: Path to TRF output file
        min_scaffold_size: Minimum scaffold size to check (default: 1Mb)
        min_expected_trs: Minimum expected number of TRs for large scaffolds
        match_first_word: If True, match only first word of scaffold names (default: True)
        debug: If True, show debug information
        
    Returns:
        Tuple of (missing_scaffolds, low_tr_scaffolds, statistics)
    """
    # Get scaffold lengths from FASTA
    scaffold_lengths = get_scaffold_lengths(fasta_file, match_first_word=match_first_word)
    
    # Get TRF results
    scaffold_trf_counts = get_trf_scaffolds(trf_file, match_first_word=match_first_word, debug=debug)
    
    # Debug: show first few scaffold names from FASTA
    if debug:
        logger.debug(f"\nDebug - First scaffold names from FASTA:")
        for name in list(scaffold_lengths.keys())[:3]:
            logger.debug(f"  - '{name}' (length: {scaffold_lengths[name]:,} bp)")
    
    # Check if there's a mismatch in scaffold name format
    fasta_scaffolds = set(scaffold_lengths.keys())
    trf_scaffolds = set(scaffold_trf_counts.keys())
    
    # Find scaffolds that are in FASTA but not in TRF (potential naming mismatch)
    only_in_fasta = fasta_scaffolds - trf_scaffolds
    if only_in_fasta and len(only_in_fasta) < len(fasta_scaffolds):
        logger.warning(f"\n⚠️  Warning: {len(only_in_fasta)} scaffolds from FASTA not found in TRF results.")
        logger.warning("  This might be a scaffold naming mismatch issue.")
        logger.warning(f"  Example FASTA scaffolds not in TRF: {list(only_in_fasta)[:3]}")
    
    # Analyze consistency
    missing_scaffolds = []  # Scaffolds with no TRF results
    low_tr_scaffolds = []   # Scaffolds with suspiciously few TRs
    small_scaffolds_skipped = 0
    
    for scaffold_name, length in scaffold_lengths.items():
        # Skip small scaffolds
        if length < min_scaffold_size:
            small_scaffolds_skipped += 1
            continue
        
        tr_count = scaffold_trf_counts.get(scaffold_name, 0)
        
        if tr_count == 0:
            missing_scaffolds.append({
                'name': scaffold_name,
                'length': length,
                'expected_trs': 'at least some'
            })
        elif tr_count < min_expected_trs:
            low_tr_scaffolds.append({
                'name': scaffold_name,
                'length': length,
                'tr_count': tr_count,
                'expected': f'>= {min_expected_trs}'
            })
    
    # Calculate statistics
    total_scaffolds = len(scaffold_lengths)
    large_scaffolds = sum(1 for l in scaffold_lengths.values() if l >= min_scaffold_size)
    scaffolds_with_trs = len(scaffold_trf_counts)
    
    statistics = {
        'total_scaffolds': total_scaffolds,
        'large_scaffolds': large_scaffolds,
        'small_scaffolds_skipped': small_scaffolds_skipped,
        'scaffolds_with_trs': scaffolds_with_trs,
        'min_scaffold_size': min_scaffold_size,
        'total_genome_size': sum(scaffold_lengths.values()),
        'large_genome_size': sum(l for l in scaffold_lengths.values() if l >= min_scaffold_size)
    }
    
    return missing_scaffolds, low_tr_scaffolds, statistics


def print_report(missing_scaffolds, low_tr_scaffolds, statistics):
    """Print consistency check report.
    
    Args:
        missing_scaffolds: List of scaffolds with no TRF results
        low_tr_scaffolds: List of scaffolds with few TRs
        statistics: Dictionary with statistics
    """
    logger.info("\n" + "="*60)
    logger.info("TRF CONSISTENCY CHECK REPORT")
    logger.info("="*60)
    
    # Print statistics
    logger.info("\n📊 Overall Statistics:")
    logger.info(f"  Total scaffolds: {statistics['total_scaffolds']:,}")
    logger.info(f"  Total genome size: {statistics['total_genome_size']:,} bp")
    logger.info(f"  Large scaffolds (>={statistics['min_scaffold_size']:,} bp): {statistics['large_scaffolds']:,}")
    logger.info(f"  Large scaffolds size: {statistics['large_genome_size']:,} bp")
    logger.info(f"  Scaffolds with TRs: {statistics['scaffolds_with_trs']:,}")
    logger.info(f"  Small scaffolds skipped: {statistics['small_scaffolds_skipped']:,}")
    
    # Report issues
    if missing_scaffolds or low_tr_scaffolds:
        logger.warning("\n⚠️  POTENTIAL ISSUES DETECTED:")

        if missing_scaffolds:
            logger.error(f"\n❌ {len(missing_scaffolds)} large scaffold(s) with NO tandem repeats detected:")
            for scaffold in missing_scaffolds[:10]:  # Show first 10
                logger.error(f"    - {scaffold['name']}: {scaffold['length']:,} bp")
            if len(missing_scaffolds) > 10:
                logger.error(f"    ... and {len(missing_scaffolds) - 10} more")

        if low_tr_scaffolds:
            logger.warning(f"\n⚠️  {len(low_tr_scaffolds)} large scaffold(s) with suspiciously few tandem repeats:")
            for scaffold in low_tr_scaffolds[:10]:  # Show first 10
                logger.warning(f"    - {scaffold['name']}: {scaffold['length']:,} bp, only {scaffold['tr_count']} TR(s)")
            if len(low_tr_scaffolds) > 10:
                logger.warning(f"    ... and {len(low_tr_scaffolds) - 10} more")

        logger.info("\n🔍 Possible causes:")
        logger.info("  1. TRF failed to process some files (check for signal 6 errors)")
        logger.info("  2. Some scaffolds genuinely have no/few tandem repeats")
        logger.info("  3. TRF parameters may need adjustment for this genome")
        logger.info("  4. Incomplete TRF run (use --continue-on-error with caution)")

        logger.info("\n💡 Recommendations:")
        logger.info("  1. Check TRF logs for errors")
        logger.info("  2. Re-run TRF for missing scaffolds with --force flag")
        logger.info("  3. Consider adjusting TRF parameters if needed")
        logger.info("  4. Manually inspect problematic scaffolds")
        
        return False  # Consistency check failed
    else:
        logger.info("\n✅ CONSISTENCY CHECK PASSED")
        logger.info("All large scaffolds have tandem repeats detected!")
        return True  # Consistency check passed


def save_report(missing_scaffolds, low_tr_scaffolds, statistics, output_file):
    """Save detailed report to file.
    
    Args:
        missing_scaffolds: List of scaffolds with no TRF results
        low_tr_scaffolds: List of scaffolds with few TRs
        statistics: Dictionary with statistics
        output_file: Output file path
    """
    with open(output_file, 'w') as f:
        f.write("TRF CONSISTENCY CHECK DETAILED REPORT\n")
        f.write("="*60 + "\n\n")
        
        # Write statistics
        f.write("STATISTICS:\n")
        for key, value in statistics.items():
            f.write(f"  {key}: {value}\n")
        f.write("\n")
        
        # Write missing scaffolds
        if missing_scaffolds:
            f.write(f"SCAFFOLDS WITH NO TANDEM REPEATS ({len(missing_scaffolds)}):\n")
            f.write("scaffold_name\tlength_bp\n")
            for scaffold in missing_scaffolds:
                f.write(f"{scaffold['name']}\t{scaffold['length']}\n")
            f.write("\n")
        
        # Write low TR scaffolds
        if low_tr_scaffolds:
            f.write(f"SCAFFOLDS WITH FEW TANDEM REPEATS ({len(low_tr_scaffolds)}):\n")
            f.write("scaffold_name\tlength_bp\ttr_count\n")
            for scaffold in low_tr_scaffolds:
                f.write(f"{scaffold['name']}\t{scaffold['length']}\t{scaffold['tr_count']}\n")
            f.write("\n")
        
        # Write summary
        if missing_scaffolds or low_tr_scaffolds:
            f.write("RESULT: FAILED - Issues detected\n")
        else:
            f.write("RESULT: PASSED - All large scaffolds have tandem repeats\n")
    
    logger.info(f"\n📄 Detailed report saved to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Check TRF results consistency - verify all large scaffolds have tandem repeats"
    )
    parser.add_argument(
        "-f", "--fasta", 
        help="Input FASTA file", 
        required=True
    )
    parser.add_argument(
        "-t", "--trf", 
        help="TRF output file", 
        required=True
    )
    parser.add_argument(
        "-s", "--min-size", 
        help="Minimum scaffold size to check in bp [1000000]", 
        type=int, 
        default=1000000
    )
    parser.add_argument(
        "-m", "--min-trs", 
        help="Minimum expected TRs for large scaffolds [1]", 
        type=int, 
        default=1
    )
    parser.add_argument(
        "-o", "--output", 
        help="Output report file [optional]", 
        default=None
    )
    parser.add_argument(
        "--no-match-first-word",
        help="Use full scaffold names instead of just first word",
        action="store_true",
        default=False
    )
    parser.add_argument(
        "--debug",
        help="Show debug information about scaffold name matching",
        action="store_true",
        default=False
    )
    
    args = parser.parse_args()
    
    # Check input files exist
    if not os.path.exists(args.fasta):
        logger.error(f"FASTA file not found: {args.fasta}")
        sys.exit(1)

    if not os.path.exists(args.trf):
        logger.error(f"TRF file not found: {args.trf}")
        sys.exit(1)
    
    # Run consistency check
    try:
        missing_scaffolds, low_tr_scaffolds, statistics = check_consistency(
            args.fasta, 
            args.trf, 
            args.min_size,
            args.min_trs,
            match_first_word=not args.no_match_first_word,
            debug=args.debug
        )
        
        # Print report to console
        passed = print_report(missing_scaffolds, low_tr_scaffolds, statistics)
        
        # Save detailed report if requested
        if args.output:
            save_report(missing_scaffolds, low_tr_scaffolds, statistics, args.output)
        
        # Exit with appropriate code
        sys.exit(0 if passed else 1)
        
    except Exception as e:
        logger.error(f"Error during consistency check: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(2)


if __name__ == "__main__":
    main()