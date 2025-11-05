#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 2025-01-05
# @author: Claude (SAT-49)
# @description: Tools for working with BED format files

"""
BED format sequence extraction tools.

This module provides functionality to extract sequences from FASTA files
based on BED format coordinates and create TRF-like output files.
"""

import logging
import os
from satellome.core_functions.io.fasta_file import sc_iter_fasta_brute

logger = logging.getLogger(__name__)


def reverse_complement(seq):
    """
    Return reverse complement of DNA sequence.

    Args:
        seq (str): DNA sequence

    Returns:
        str: Reverse complement of the input sequence

    Examples:
        >>> reverse_complement("ATCG")
        'CGAT'
        >>> reverse_complement("AAAA")
        'TTTT'
    """
    complement = {
        'A': 'T', 'T': 'A',
        'C': 'G', 'G': 'C',
        'a': 't', 't': 'a',
        'c': 'g', 'g': 'c',
        'N': 'N', 'n': 'n'
    }
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


def extract_sequences_from_bed(fasta_file, bed_file, output_file):
    """
    Extract sequences from FASTA file based on BED coordinates and save to TRF-like format.

    Memory-efficient implementation: reads BED first, sorts by chromosome, then processes
    FASTA sequentially chromosome by chromosome.

    BED coordinates are 0-based, half-open intervals [start, end).
    Python slicing works the same way, so seq[start:end] extracts correctly.

    Args:
        fasta_file (str): Path to input FASTA file
        bed_file (str): Path to input BED file from FasTAN/tanbed
        output_file (str): Path to output file (BED + sequence column)

    Returns:
        int: Number of sequences extracted

    Format:
        Input BED: chr  start  end  name  score  strand  [other_columns...]
        Output: chr  start  end  name  score  strand  [other_columns...]  sequence
    """
    logger.info(f"Extracting sequences from {fasta_file} using coordinates from {bed_file}")

    # Step 1: Read and parse BED file, group by chromosome
    logger.info("Reading BED file...")
    bed_entries = {}  # {chr_name: [(line, start, end, strand), ...]}
    skipped_count = 0

    with open(bed_file, 'r') as bed_fh:
        for line in bed_fh:
            line = line.strip()

            # Skip empty lines and comments
            if not line or line.startswith('#'):
                continue

            fields = line.split('\t')

            # BED format requires at least 3 columns: chr, start, end
            if len(fields) < 3:
                logger.warning(f"Skipping invalid BED line (< 3 columns): {line[:50]}")
                skipped_count += 1
                continue

            # Take only first word from chromosome name (to match FASTA parsing)
            chr_name = fields[0].split()[0]
            try:
                start = int(fields[1])  # 0-based
                end = int(fields[2])    # 0-based, exclusive
            except ValueError:
                logger.warning(f"Skipping line with invalid coordinates: {line[:50]}")
                skipped_count += 1
                continue

            # Get strand if available (column 6 in BED6 format)
            strand = fields[5] if len(fields) > 5 else '+'

            # Basic validation
            if start >= end or start < 0:
                logger.warning(f"Invalid coordinates: start={start}, end={end}, skipping")
                skipped_count += 1
                continue

            # Add to chromosome group
            if chr_name not in bed_entries:
                bed_entries[chr_name] = []
            bed_entries[chr_name].append((line, start, end, strand))

    # Sort entries within each chromosome by start position
    for chr_name in bed_entries:
        bed_entries[chr_name].sort(key=lambda x: x[1])  # Sort by start position

    logger.info(f"Loaded {sum(len(v) for v in bed_entries.values())} BED entries for {len(bed_entries)} chromosomes")

    # Step 2: Process FASTA sequentially, extracting sequences for each chromosome
    extracted_count = 0
    seen_chromosomes = set()  # Track chromosome names to detect duplicates

    with open(output_file, 'w') as out_fh:
        # Write header
        out_fh.write("# FasTAN results with extracted sequences\n")
        out_fh.write("# Format: chr\tstart\tend\tname\tscore\tstrand\t[...]\tsequence\n")
        out_fh.write(f"# Source FASTA: {os.path.basename(fasta_file)}\n")
        out_fh.write(f"# Source BED: {os.path.basename(bed_file)}\n")

        logger.info("Processing FASTA file...")
        for header, sequence in sc_iter_fasta_brute(fasta_file):
            # Remove '>' and take first word as chromosome name
            full_header = header.lstrip('>')
            chr_name = full_header.split()[0]

            # Check for duplicate chromosome names (CRITICAL SAFETY CHECK)
            if chr_name in seen_chromosomes:
                logger.error(
                    f"Duplicate chromosome name detected: '{chr_name}'\n"
                    f"Full header: {full_header}\n"
                    f"This indicates ambiguous chromosome naming in FASTA file.\n"
                    f"Extraction results may be incorrect!"
                )
                raise ValueError(
                    f"Duplicate chromosome name '{chr_name}' found in FASTA. "
                    f"First word of FASTA headers must be unique. "
                    f"Please use unique chromosome identifiers."
                )

            seen_chromosomes.add(chr_name)

            # Skip if no BED entries for this chromosome
            if chr_name not in bed_entries:
                logger.debug(f"No BED entries for {chr_name}, skipping")
                continue

            sequence = sequence.upper()
            chr_len = len(sequence)
            logger.debug(f"Processing {chr_name}: {chr_len} bp, {len(bed_entries[chr_name])} entries")

            # Process all BED entries for this chromosome
            for line, start, end, strand in bed_entries[chr_name]:
                # Validate coordinates against actual sequence length
                if end > chr_len:
                    logger.warning(
                        f"Coordinates out of bounds for {chr_name}: end={end} > chr_len={chr_len}, skipping"
                    )
                    skipped_count += 1
                    continue

                # Extract sequence
                # BED coordinates are 0-based, half-open [start, end)
                # Python slice seq[start:end] works directly!
                extracted_seq = sequence[start:end]

                # Apply reverse complement if on negative strand
                if strand == '-':
                    extracted_seq = reverse_complement(extracted_seq)

                # Write BED fields + sequence
                out_fh.write(line + '\t' + extracted_seq + '\n')
                extracted_count += 1

    logger.info(f"✓ Extracted {extracted_count} sequences")
    if skipped_count > 0:
        logger.warning(f"Skipped {skipped_count} entries due to errors")

    return extracted_count
