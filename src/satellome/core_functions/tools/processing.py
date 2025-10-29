#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.02.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

from satellome.core_functions.io.fasta_file import sc_iter_fasta_brute
from tqdm import tqdm
import logging

logger = logging.getLogger(__name__)

REVCOMP_DICTIONARY = dict(zip("ATCGNatcgn~[]", "TAGCNtagcn~]["))


def get_gc_content(sequence):
    """Calculate GC content as a fraction (0.0 to 1.0).

    Counts G and C nucleotides (case-insensitive) and returns their
    fraction relative to the total sequence length.

    Args:
        sequence (str): DNA sequence string (case-insensitive)

    Returns:
        float: GC fraction from 0.0 to 1.0. Returns 0.0 for empty sequences.

    Examples:
        >>> get_gc_content("ATGC")
        0.5
        >>> get_gc_content("AAAA")
        0.0
        >>> get_gc_content("GGCC")
        1.0
        >>> get_gc_content("AtGc")
        0.5
        >>> get_gc_content("")
        0.0
    """
    length = len(sequence)
    if not length:
        return 0.0

    gc_count = (sequence.count('G') + sequence.count('g') +
                sequence.count('C') + sequence.count('c'))
    return float(gc_count) / float(length)


def get_gc_percentage(sequence):
    """Calculate GC content as a percentage (0 to 100).

    Args:
        sequence (str): DNA sequence string (case-insensitive)

    Returns:
        float: GC percentage from 0.0 to 100.0

    Example:
        >>> get_gc_percentage("ATGC")
        50.0
    """
    return get_gc_content(sequence) * 100.0


def get_revcomp(sequence):
    """Return reverse complementary sequence.

    >>> complementary('AT CG')
    'CGAT'

    """
    return "".join(
        REVCOMP_DICTIONARY.get(nucleotide, "") for nucleotide in reversed(sequence)
    )

def get_genome_size(fasta_file):
    ''' Compute genome size from fasta file.'''

    logger.info("Calculating genome size...")
    genome_size = 0
    for _, seq in sc_iter_fasta_brute(fasta_file):
        genome_size += len(seq)
    logger.info(f"Genome size: {genome_size} bp")
    return genome_size


def get_genome_size_with_progress(fasta_file):
    '''Compute genome size from fasta file with progress bar.'''
    logger.info(f"Calculating genome size for: {fasta_file}")
    genome_size = 0
    n_seqs = 0
    # First, count number of sequences for tqdm
    seq_headers = [h for h, _ in sc_iter_fasta_brute(fasta_file)]
    with tqdm(total=len(seq_headers), desc="Genome size (scaffolds)") as pbar:
        for _, seq in sc_iter_fasta_brute(fasta_file):
            genome_size += len(seq)
            n_seqs += 1
            pbar.update(1)
    logger.info(f"Total genome size: {genome_size:,} bp in {n_seqs} scaffolds/contigs")
    return genome_size


def count_lines_large_file(filename, chunk_size=1024*1024):
    line_count = 0
    with open(filename, 'rb') as f:
        while chunk := f.read(chunk_size):
            line_count += chunk.count(b'\n')
    return line_count

