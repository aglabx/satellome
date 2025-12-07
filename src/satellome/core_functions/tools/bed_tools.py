#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 2025-01-05
# @author: Claude (SAT-49)
# @description: Tools for working with BED format files

"""
BED format sequence extraction and manipulation tools.

Provides utilities for extracting genomic sequences from FASTA files based on
BED format coordinates (from FasTAN/tanbed output). Handles strand orientation,
chromosome name parsing, and creates TRF-compatible output with sequences.

Functions:
    reverse_complement: Compute DNA reverse complement
    extract_sequences_from_bed: Extract and annotate sequences from BED coordinates

Key Features:
    - Memory-efficient chromosome-by-chromosome processing
    - BED6 format support with strand orientation
    - Automatic reverse complement for negative strand
    - Duplicate chromosome name validation
    - Short chromosome name extraction (first word only)
    - Comprehensive logging and error handling
    - TRF-compatible output format

BED Format Notes:
    - Coordinates are 0-based, half-open intervals [start, end)
    - Python slicing seq[start:end] works directly with BED coords
    - Strand column (6th field): '+' forward, '-' reverse complement
    - Chromosome names: Uses first whitespace-delimited word only

Output Format:
    - Original BED columns + repeat_length + extracted_sequence
    - Header comments with source file information
    - Chromosome names normalized to first word

Example:
    >>> # Extract sequences from BED annotations
    >>> from satellome.core_functions.tools.bed_tools import extract_sequences_from_bed
    >>> count = extract_sequences_from_bed(
    ...     "genome.fasta",
    ...     "tandem_repeats.bed",
    ...     "repeats_with_seqs.txt"
    ... )
    INFO:satellome.core_functions.tools.bed_tools:Extracting sequences from genome.fasta...
    INFO:satellome.core_functions.tools.bed_tools:✓ Extracted 1523 sequences
    >>> print(f"Extracted {count} repeat sequences")
    Extracted 1523 repeat sequences

Typical Use Case:
    1. Run FasTAN/tanbed on genome to find tandem repeats (generates BED)
    2. Use extract_sequences_from_bed() to add sequences to BED output
    3. Result file has coordinates + actual repeat sequences for analysis

See Also:
    satellome.core_functions.io.fasta_file: FASTA iteration utilities
    satellome.steps.trf_parse_raw: TRF output parsing
"""

import logging
import os
from satellome.core_functions.io.fasta_file import sc_iter_fasta_brute
from satellome.core_functions.tools.processing import get_gc_content

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


def extract_sequences_from_bed(fasta_file, bed_file, output_file, fasta_output_file=None, project="FasTAN"):
    """
    Extract sequences from FASTA based on BED coordinates and output TRF-compatible format.

    Memory-efficient implementation that reads BED file once, groups by chromosome,
    then processes FASTA sequentially chromosome-by-chromosome.

    Args:
        fasta_file (str): Path to input genome FASTA file
        bed_file (str): Path to BED file with tandem repeat coordinates (FasTAN/tanbed output)
        output_file (str): Path to output TRF-format file (18 tab-separated fields)
        fasta_output_file (str, optional): Path to output FASTA file with extracted sequences.
            Headers format: >chr_start_end_length_period. If None, no FASTA file is created.
        project (str): Project name for TRF output (default: "FasTAN")

    Returns:
        int: Number of sequences successfully extracted

    Raises:
        ValueError: If duplicate chromosome names detected in FASTA (ambiguous mapping)

    Example:
        >>> # Extract repeat sequences from genome
        >>> count = extract_sequences_from_bed(
        ...     "hg38.fasta",
        ...     "tandem_repeats.bed",
        ...     "repeats_annotated.txt"
        ... )
        INFO:...Loaded 1523 BED entries for 24 chromosomes
        INFO:...✓ Extracted 1523 sequences
        >>> print(count)
        1523

    Input BED Format (6+ columns):
        - Column 1: Chromosome name (uses first word only)
        - Column 2: Start position (0-based, inclusive)
        - Column 3: End position (0-based, exclusive)
        - Column 4: Feature name
        - Column 5: Score
        - Column 6: Strand ('+' or '-')
        - Columns 7+: Additional fields (preserved in output)

    Output Format:
        - All original BED columns
        - Column N+1: repeat_length (end - start)
        - Column N+2: extracted_sequence (reverse complemented if strand '-')
        - Header comments with source file information

    Processing Steps:
        1. Load entire BED file, group entries by chromosome
        2. Sort entries within each chromosome by start position
        3. Iterate through FASTA chromosomes sequentially
        4. For each chromosome, extract all BED regions at once
        5. Apply reverse complement if strand is '-'
        6. Write annotated BED lines with sequences

    Coordinate System:
        - BED uses 0-based, half-open intervals: [start, end)
        - Python slicing seq[start:end] works directly (no adjustment needed)
        - Example: BED "chr1 100 105" extracts bases at positions 100,101,102,103,104

    Validation and Error Handling:
        - Skips lines with < 3 columns or invalid coordinates
        - Validates end <= chromosome_length
        - Raises ValueError on duplicate chromosome names (critical safety check)
        - Logs warnings for skipped entries with reasons
        - Strand defaults to '+' if not specified (column 6 missing)

    Note:
        - Uses chromosome short names (first whitespace-delimited word only)
        - Handles both uppercase and lowercase DNA sequences
        - Memory efficient: only one chromosome loaded at a time
        - Output includes header comments for reproducibility
        - Created for SAT-49: BED sequence extraction feature
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
    trf_id_counter = 0

    # Open output files
    out_fh = open(output_file, 'w')
    fasta_fh = open(fasta_output_file, 'w') if fasta_output_file else None

    try:
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

                # Extract sequence (BED coordinates are 0-based, half-open [start, end))
                extracted_seq = sequence[start:end]

                # Apply reverse complement if on negative strand
                if strand == '-':
                    extracted_seq = reverse_complement(extracted_seq)

                # Get period from BED fields (column 4)
                fields = line.split('\t')
                period = int(fields[3]) if len(fields) > 3 else 1

                # Calculate TRF fields
                trf_id_counter += 1
                trf_array_length = end - start
                trf_n_copy = round(trf_array_length / period, 1) if period > 0 else 0
                trf_consensus = extracted_seq[:period] if period <= len(extracted_seq) else extracted_seq
                trf_array_gc = round(get_gc_content(extracted_seq), 2)
                trf_consensus_gc = round(get_gc_content(trf_consensus), 2)

                # TRF format: 18 tab-separated fields
                # project, trf_id, trf_head, trf_l_ind, trf_r_ind, trf_period, trf_n_copy,
                # trf_pmatch, trf_pvar, trf_entropy, trf_consensus, trf_array,
                # trf_array_gc, trf_consensus_gc, trf_array_length, trf_joined, trf_family, trf_ref_annotation
                trf_fields = [
                    project,                    # project
                    str(trf_id_counter),        # trf_id
                    chr_name,                   # trf_head
                    str(start + 1),             # trf_l_ind (1-based for TRF compatibility)
                    str(end),                   # trf_r_ind
                    str(period),                # trf_period
                    str(trf_n_copy),            # trf_n_copy
                    "-1",                       # trf_pmatch (not available from FasTAN)
                    "-1",                       # trf_pvar (not available from FasTAN)
                    "-1",                       # trf_entropy (not available from FasTAN)
                    trf_consensus,              # trf_consensus
                    extracted_seq,              # trf_array
                    str(trf_array_gc),          # trf_array_gc
                    str(trf_consensus_gc),      # trf_consensus_gc
                    str(trf_array_length),      # trf_array_length
                    "",                         # trf_joined
                    "",                         # trf_family
                    "",                         # trf_ref_annotation
                ]
                out_fh.write('\t'.join(trf_fields) + '\n')

                # Write FASTA output if requested
                # Header format: >chr_start_end_length_period
                if fasta_fh:
                    fasta_header = f">{chr_name}_{start}_{end}_{trf_array_length}_{period}"
                    fasta_fh.write(f"{fasta_header}\n{extracted_seq}\n")

                extracted_count += 1

    finally:
        out_fh.close()
        if fasta_fh:
            fasta_fh.close()

    logger.info(f"✓ Extracted {extracted_count} sequences to TRF format")
    if fasta_output_file:
        logger.info(f"✓ FASTA file created: {fasta_output_file}")
    if skipped_count > 0:
        logger.warning(f"Skipped {skipped_count} entries due to errors")

    return extracted_count
