#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.02.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

from satellome.core_functions.io.fasta_file import sc_iter_fasta_brute
from tqdm import tqdm

REVCOMP_DICTIONARY = dict(zip("ATCGNatcgn~[]", "TAGCNtagcn~]["))


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

    print("Genome size:", end=" ")
    genome_size = 0
    for _, seq in sc_iter_fasta_brute(fasta_file):
        genome_size += len(seq)
    print(f"{genome_size} bp")
    return genome_size


def get_genome_size_with_progress(fasta_file):
    '''Compute genome size from fasta file with progress bar.'''
    print(f"Calculating genome size for: {fasta_file}")
    genome_size = 0
    n_seqs = 0
    # First, count number of sequences for tqdm
    seq_headers = [h for h, _ in sc_iter_fasta_brute(fasta_file)]
    with tqdm(total=len(seq_headers), desc="Genome size (scaffolds)") as pbar:
        for _, seq in sc_iter_fasta_brute(fasta_file):
            genome_size += len(seq)
            n_seqs += 1
            pbar.update(1)
    print(f"Total genome size: {genome_size:,} bp in {n_seqs} scaffolds/contigs")
    return genome_size


def count_lines_large_file(filename, chunk_size=1024*1024):
    line_count = 0
    with open(filename, 'rb') as f:
        while chunk := f.read(chunk_size):
            line_count += chunk.count(b'\n')
    return line_count

