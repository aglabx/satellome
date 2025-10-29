"""Sample data for testing Satellome functions.

This module provides reusable test data that can be imported
directly in tests without using fixtures.
"""

# Sample FASTA headers in various formats
FASTA_HEADERS = {
    'simple': '>chr1',
    'with_description': '>chr1 Homo sapiens chromosome 1',
    'ncbi_format': '>NC_000001.11 Homo sapiens chromosome 1, GRCh38.p14 Primary Assembly',
    'ensembl_format': '>1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF',
    'custom': '>scaffold_1_length_12345',
}

# Sample TRF lines (tab-separated)
TRF_LINES = {
    'microsatellite_di': 'chr1\t1000\t1100\t2\t50.0\t95\t0\t0\t0\tAT\t' + 'AT' * 50,
    'microsatellite_tri': 'chr1\t2000\t2150\t3\t50.0\t90\t5\t5\t5\tATG\t' + 'ATG' * 50,
    'complex': 'chr1\t5000\t5200\t20\t10.0\t85\t10\t10\t10\tATGCATGCATGCATGCATGC\t' + 'ATGCATGCATGCATGCATGC' * 10,
}

# Chromosome name patterns for testing parsing
CHROMOSOME_NAMES = {
    'standard': ['>chr1', '>chr2', '>chrX', '>chrY', '>chrM'],
    'numeric': ['>1', '>2', '>22', '>X', '>Y'],
    'ncbi': ['>NC_000001.11', '>NC_000002.12', '>NC_012920.1'],
    'scaffold': ['>scaffold_1', '>contig_1', '>unplaced_scaffold_1'],
}

# GC content test cases: (sequence, expected_gc_fraction)
GC_TEST_CASES = [
    ('ATGC', 0.5),
    ('AAAA', 0.0),
    ('GGCC', 1.0),
    ('', 0.0),
    ('A', 0.0),
    ('G', 1.0),
    ('ATGATGATG', 1.0/3),  # 3 G's out of 9
]

# Reverse complement test cases: (input, expected_output)
REVCOMP_TEST_CASES = [
    ('ATGC', 'GCAT'),
    ('AAAA', 'TTTT'),
    ('GGCC', 'GGCC'),
    ('GAATTC', 'GAATTC'),  # Palindrome
    ('', ''),
    ('N', 'N'),
]

# Classification test cases: (period, consensus, expected_class)
CLASSIFICATION_TEST_CASES = [
    (1, 'A', 'micro'),
    (2, 'AT', 'micro'),
    (9, 'ATGATGATG', 'micro'),
    (10, 'ATGATGATGA', 'complex'),
    (50, 'A' * 50, 'complex'),
]
