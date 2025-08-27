#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 27.08.2024
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

"""
Smart genome splitting using varprofiler to identify repeat-rich regions
and optimize TRF processing.
"""

import os
import subprocess
import tempfile
from typing import List, Tuple, Dict
from tqdm import tqdm

from satellome.core_functions.io.fasta_file import sc_iter_fasta_brute


def run_varprofiler(
    fasta_file: str,
    output_bed: str,
    kmer_size: int = 17,
    window_size: int = 100000,
    step_size: int = 25000,
    threads: int = 20
) -> str:
    """
    Run varprofiler to identify unique k-mer density across genome.
    
    Args:
        fasta_file: Input FASTA file path
        output_bed: Output BED file path
        kmer_size: K-mer size (default: 17)
        window_size: Window size in bp (default: 100kb)
        step_size: Step size in bp (default: 25kb)
        threads: Number of threads (default: 20)
    
    Returns:
        Path to output BED file
    """
    command = [
        "varprofiler",
        fasta_file,
        output_bed,
        str(kmer_size),
        str(window_size),
        str(step_size),
        str(threads)
    ]
    
    try:
        subprocess.run(command, check=True, capture_output=True)
        return output_bed
    except subprocess.CalledProcessError as e:
        print(f"Error running varprofiler: {e}")
        print(f"stderr: {e.stderr.decode()}")
        raise


def parse_kmer_bed(bed_file: str) -> Dict[str, List[Tuple[int, int, int]]]:
    """
    Parse BED file from varprofiler.
    
    Args:
        bed_file: Path to BED file
    
    Returns:
        Dictionary with chromosome as key and list of (start, end, unique_kmers) tuples
    """
    regions = {}
    
    with open(bed_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            unique_kmers = int(parts[3])
            
            if chrom not in regions:
                regions[chrom] = []
            
            regions[chrom].append((start, end, unique_kmers))
    
    return regions


def identify_repeat_regions(
    regions: List[Tuple[int, int, int]],
    threshold: int = 90000,
    min_region_size: int = 50000
) -> List[Tuple[int, int]]:
    """
    Identify regions likely to contain repeats based on low unique k-mer count.
    
    Args:
        regions: List of (start, end, unique_kmers) tuples
        threshold: Threshold for unique k-mers (regions below this are repeat-rich)
        min_region_size: Minimum size of region to process
    
    Returns:
        List of (start, end) tuples for repeat-rich regions
    """
    repeat_regions = []
    current_start = None
    current_end = None
    
    for start, end, unique_kmers in sorted(regions):
        if unique_kmers < threshold:
            # This region is repeat-rich
            if current_start is None:
                current_start = start
                current_end = end
            else:
                # Extend current region
                current_end = max(current_end, end)
        else:
            # This region has many unique k-mers (likely no repeats)
            if current_start is not None and (current_end - current_start) >= min_region_size:
                repeat_regions.append((current_start, current_end))
            current_start = None
            current_end = None
    
    # Don't forget the last region
    if current_start is not None and (current_end - current_start) >= min_region_size:
        repeat_regions.append((current_start, current_end))
    
    return repeat_regions


def merge_overlapping_regions(
    regions: List[Tuple[int, int]],
    gap_threshold: int = 100000
) -> List[Tuple[int, int]]:
    """
    Merge overlapping or nearby regions.
    
    Args:
        regions: List of (start, end) tuples
        gap_threshold: Maximum gap to merge regions
    
    Returns:
        List of merged (start, end) tuples
    """
    if not regions:
        return []
    
    sorted_regions = sorted(regions)
    merged = [sorted_regions[0]]
    
    for start, end in sorted_regions[1:]:
        last_start, last_end = merged[-1]
        
        if start <= last_end + gap_threshold:
            # Merge regions
            merged[-1] = (last_start, max(end, last_end))
        else:
            # Add as new region
            merged.append((start, end))
    
    return merged


def split_genome_smart(
    fasta_file: str,
    wdir: str,
    project: str,
    threads: int = 20,
    kmer_threshold: int = 90000,
    chunk_size: int = None,  # Not used - we process entire repeat-rich regions
    overlap_size: int = 0,   # No overlap needed
    use_kmer_filter: bool = True,
    kmer_bed_file: str = None
) -> List[str]:
    """
    Split genome into chunks using varprofiler to skip repeat-poor regions.
    
    Args:
        fasta_file: Input FASTA file
        wdir: Working directory
        project: Project name
        threads: Number of threads
        kmer_threshold: Threshold for unique k-mers
        chunk_size: Size of chunks in bp
        overlap_size: Overlap between chunks in bp
        use_kmer_filter: Whether to use k-mer filtering
    
    Returns:
        List of output file paths
    """
    folder_path = tempfile.mkdtemp(dir=wdir)
    output_files = []
    
    # Step 1: Run k-mer profiling if enabled
    repeat_regions = {}
    if use_kmer_filter:
        # Use provided BED file or generate new one
        if kmer_bed_file and os.path.exists(kmer_bed_file):
            print(f"Using provided k-mer profile: {kmer_bed_file}")
            bed_file = kmer_bed_file
        else:
            print("Running varprofiler to identify repeat-rich regions...")
            bed_file = os.path.join(folder_path, f"{project}.varprofile.bed")
            
            try:
                run_varprofiler(fasta_file, bed_file, threads=threads)
            except Exception as e:
                print(f"Warning: k-mer profiling failed, falling back to standard splitting: {e}")
                use_kmer_filter = False
        
        if use_kmer_filter:
            try:
                kmer_data = parse_kmer_bed(bed_file)
                
                # Identify repeat-rich regions for each chromosome
                for chrom, regions in kmer_data.items():
                    repeat_rich = identify_repeat_regions(regions, threshold=kmer_threshold)
                    merged = merge_overlapping_regions(repeat_rich)
                    repeat_regions[chrom] = merged
                
                print(f"Identified {sum(len(r) for r in repeat_regions.values())} repeat-rich regions")
            except Exception as e:
                print(f"Warning: k-mer profiling failed, falling back to standard splitting: {e}")
                use_kmer_filter = False
    
    # Step 2: Split sequences into chunks
    file_counter = 0
    
    print("Splitting genome into chunks...")
    for header, seq in sc_iter_fasta_brute(fasta_file):
        chrom = header.split()[0].replace(">", "")
        seq_len = len(seq)
        
        if use_kmer_filter and chrom in repeat_regions:
            # Use k-mer guided splitting
            regions_to_process = repeat_regions[chrom]
        else:
            # Process entire chromosome
            regions_to_process = [(0, seq_len)]
        
        for region_start, region_end in regions_to_process:
            # Skip very small regions
            if region_end - region_start < 1000:
                continue
            
            # Extract region sequence (entire repeat-rich region, no chunking)
            region_seq = seq[region_start:region_end]
            
            # Save to file
            output_file = os.path.join(folder_path, f"{file_counter:05d}.fa")
            with open(output_file, 'w') as fw:
                # Modify header to include original coordinates
                new_header = f"{header}__{region_start}_{region_end}"
                fw.write(f"{new_header}\n{region_seq}\n")
            
            output_files.append(output_file)
            file_counter += 1
    
    print(f"Created {len(output_files)} chunks")
    if use_kmer_filter:
        total_bp = sum(end - start for regions in repeat_regions.values() for start, end in regions)
        print(f"Total bp in repeat-rich regions: {total_bp:,}")
    
    return output_files


def restore_coordinates(trf_line: str) -> str:
    """
    Restore original coordinates from chunk coordinates.
    
    Args:
        trf_line: TRF output line with modified header
    
    Returns:
        TRF line with restored coordinates
    """
    parts = trf_line.strip().split('\t')
    if len(parts) < 3:
        return trf_line
    
    # Check if header contains coordinate info
    header = parts[0]
    if '__' in header:
        base_header, coords = header.rsplit('__', 1)
        if '_' in coords:
            try:
                chunk_start, chunk_end = map(int, coords.split('_'))
                
                # Adjust TRF coordinates
                trf_start = int(parts[1])
                trf_end = int(parts[2])
                
                parts[0] = base_header
                parts[1] = str(trf_start + chunk_start)
                parts[2] = str(trf_end + chunk_start)
                
                return '\t'.join(parts) + '\n'
            except:
                pass
    
    return trf_line