#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.02.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

from satelome.core_functions.io.fasta_file import sc_iter_fasta_brute

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


def count_lines_large_file(filename, chunk_size=1024*1024):
    line_count = 0
    with open(filename, 'rb') as f:
        while chunk := f.read(chunk_size):
            line_count += chunk.count(b'\n')
    return line_count

