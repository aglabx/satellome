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

# Add parent directory to path to import satellome modules
sys.path.append(os.path.join(os.path.dirname(__file__), "../src"))

from satellome.core_functions.io.fasta_file import sc_iter_fasta_brute
from satellome.core_functions.io.tab_file import sc_iter_tab_file
from satellome.core_functions.models.trf_model import TRModel


def get_scaffold_lengths(fasta_file):
    """Get lengths of all scaffolds/contigs from FASTA file.
    
    Args:
        fasta_file: Path to input FASTA file
        
    Returns:
        Dictionary mapping scaffold name to length
    """
    scaffold_lengths = {}
    
    print(f"Reading scaffold lengths from {fasta_file}...")
    for header, sequence in sc_iter_fasta_brute(fasta_file):
        # Extract scaffold name from header
        scaffold_name = header.replace(">", "").split()[0]
        scaffold_lengths[scaffold_name] = len(sequence)
    
    return scaffold_lengths


def get_trf_scaffolds(trf_file):
    """Get all scaffolds that have TRF results.
    
    Args:
        trf_file: Path to TRF output file
        
    Returns:
        Dictionary mapping scaffold name to number of tandem repeats
    """
    scaffold_trf_counts = defaultdict(int)
    total_trs = 0
    
    print(f"Reading TRF results from {trf_file}...")
    for trf_obj in sc_iter_tab_file(trf_file, TRModel):
        scaffold_name = trf_obj.trf_gi  # or trf_obj.trf_chr depending on format
        scaffold_trf_counts[scaffold_name] += 1
        total_trs += 1
    
    print(f"Found {total_trs:,} tandem repeats across {len(scaffold_trf_counts):,} scaffolds")
    
    return scaffold_trf_counts


def check_consistency(fasta_file, trf_file, min_scaffold_size=1000000, min_expected_trs=1):
    """Check consistency between FASTA input and TRF output.
    
    Args:
        fasta_file: Path to input FASTA file
        trf_file: Path to TRF output file
        min_scaffold_size: Minimum scaffold size to check (default: 1Mb)
        min_expected_trs: Minimum expected number of TRs for large scaffolds
        
    Returns:
        Tuple of (missing_scaffolds, low_tr_scaffolds, statistics)
    """
    # Get scaffold lengths from FASTA
    scaffold_lengths = get_scaffold_lengths(fasta_file)
    
    # Get TRF results
    scaffold_trf_counts = get_trf_scaffolds(trf_file)
    
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
    print("\n" + "="*60)
    print("TRF CONSISTENCY CHECK REPORT")
    print("="*60)
    
    # Print statistics
    print("\n📊 Overall Statistics:")
    print(f"  Total scaffolds: {statistics['total_scaffolds']:,}")
    print(f"  Total genome size: {statistics['total_genome_size']:,} bp")
    print(f"  Large scaffolds (>={statistics['min_scaffold_size']:,} bp): {statistics['large_scaffolds']:,}")
    print(f"  Large scaffolds size: {statistics['large_genome_size']:,} bp")
    print(f"  Scaffolds with TRs: {statistics['scaffolds_with_trs']:,}")
    print(f"  Small scaffolds skipped: {statistics['small_scaffolds_skipped']:,}")
    
    # Report issues
    if missing_scaffolds or low_tr_scaffolds:
        print("\n⚠️  POTENTIAL ISSUES DETECTED:")
        
        if missing_scaffolds:
            print(f"\n❌ {len(missing_scaffolds)} large scaffold(s) with NO tandem repeats detected:")
            for scaffold in missing_scaffolds[:10]:  # Show first 10
                print(f"    - {scaffold['name']}: {scaffold['length']:,} bp")
            if len(missing_scaffolds) > 10:
                print(f"    ... and {len(missing_scaffolds) - 10} more")
        
        if low_tr_scaffolds:
            print(f"\n⚠️  {len(low_tr_scaffolds)} large scaffold(s) with suspiciously few tandem repeats:")
            for scaffold in low_tr_scaffolds[:10]:  # Show first 10
                print(f"    - {scaffold['name']}: {scaffold['length']:,} bp, only {scaffold['tr_count']} TR(s)")
            if len(low_tr_scaffolds) > 10:
                print(f"    ... and {len(low_tr_scaffolds) - 10} more")
        
        print("\n🔍 Possible causes:")
        print("  1. TRF failed to process some files (check for signal 6 errors)")
        print("  2. Some scaffolds genuinely have no/few tandem repeats")
        print("  3. TRF parameters may need adjustment for this genome")
        print("  4. Incomplete TRF run (use --continue-on-error with caution)")
        
        print("\n💡 Recommendations:")
        print("  1. Check TRF logs for errors")
        print("  2. Re-run TRF for missing scaffolds with --force flag")
        print("  3. Consider adjusting TRF parameters if needed")
        print("  4. Manually inspect problematic scaffolds")
        
        return False  # Consistency check failed
    else:
        print("\n✅ CONSISTENCY CHECK PASSED")
        print("All large scaffolds have tandem repeats detected!")
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
    
    print(f"\n📄 Detailed report saved to: {output_file}")


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
    
    args = parser.parse_args()
    
    # Check input files exist
    if not os.path.exists(args.fasta):
        print(f"Error: FASTA file not found: {args.fasta}")
        sys.exit(1)
    
    if not os.path.exists(args.trf):
        print(f"Error: TRF file not found: {args.trf}")
        sys.exit(1)
    
    # Run consistency check
    try:
        missing_scaffolds, low_tr_scaffolds, statistics = check_consistency(
            args.fasta, 
            args.trf, 
            args.min_size,
            args.min_trs
        )
        
        # Print report to console
        passed = print_report(missing_scaffolds, low_tr_scaffolds, statistics)
        
        # Save detailed report if requested
        if args.output:
            save_report(missing_scaffolds, low_tr_scaffolds, statistics, args.output)
        
        # Exit with appropriate code
        sys.exit(0 if passed else 1)
        
    except Exception as e:
        print(f"Error during consistency check: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(2)


if __name__ == "__main__":
    main()