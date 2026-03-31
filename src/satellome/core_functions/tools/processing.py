#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.02.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com
"""
General-purpose sequence processing utilities.

Provides basic DNA sequence operations and file processing functions for
genomic analysis. Includes GC content calculation, reverse complement,
genome size computation, and efficient line counting.

Functions:
    get_gc_content: Calculate GC fraction (0.0-1.0)
    get_gc_percentage: Calculate GC percentage (0-100)
    get_revcomp: Compute reverse complement of DNA sequence
    get_genome_size: Calculate total genome size from FASTA
    get_genome_size_with_progress: Calculate genome size with progress bar
    count_lines_large_file: Efficiently count lines in large files

Constants:
    REVCOMP_DICTIONARY: Nucleotide complement mapping (supports ATCGN + brackets)

Key Features:
    - Case-insensitive DNA operations
    - Memory-efficient FASTA iteration
    - Progress bars for long-running operations
    - Fast binary file reading for line counting
    - Support for ambiguous nucleotides (N)

Example:
    >>> from satellome.core_functions.tools.processing import get_gc_content, get_revcomp
    >>> seq = "ATCGATCG"
    >>> print(f"GC: {get_gc_content(seq):.2%}")
    GC: 50.00%
    >>> print(f"RevComp: {get_revcomp(seq)}")
    RevComp: CGATCGAT
    >>>
    >>> # Calculate genome size
    >>> from satellome.core_functions.tools.processing import get_genome_size
    >>> size = get_genome_size("genome.fasta")
    INFO:...Genome size: 3000000000 bp

Typical Use Case:
    1. Load genome FASTA file
    2. Calculate basic statistics (size, GC content)
    3. Process sequences (reverse complement for negative strand)
    4. Count lines for progress estimation

See Also:
    satellome.core_functions.io.fasta_file: FASTA file iteration
    satellome.core_functions.tools.statistics: Advanced statistics
"""

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
    """
    Compute reverse complement of DNA sequence.

    Returns the reverse complement by reversing the sequence and replacing
    each nucleotide with its Watson-Crick complement (A↔T, C↔G). Supports
    both uppercase and lowercase, plus ambiguous nucleotides (N) and brackets.

    Args:
        sequence (str): DNA sequence (case-insensitive, may include N and brackets)

    Returns:
        str: Reverse complement sequence in same case as input

    Example:
        >>> get_revcomp("ATCG")
        'CGAT'
        >>> get_revcomp("AAAA")
        'TTTT'
        >>> get_revcomp("atcg")
        'cgat'
        >>> get_revcomp("ATNGC")
        'GCNAT'
        >>> get_revcomp("AT CG")  # Spaces removed
        'CGAT'

    Note:
        - Uses REVCOMP_DICTIONARY for efficient lookups
        - Supports: A, T, C, G, N (ambiguous), brackets ~[]
        - Unknown characters are silently removed (empty string in dict)
        - Case-preserving: lowercase input → lowercase output
        - Used extensively in distance calculations for strand comparison
    """
    return "".join(
        REVCOMP_DICTIONARY.get(nucleotide, "") for nucleotide in reversed(sequence)
    )

def get_genome_size(fasta_file):
    """
    Calculate total genome size from FASTA file.

    Iterates through all sequences in FASTA file and sums their lengths.
    Memory-efficient implementation that processes sequences one at a time.

    Args:
        fasta_file (str): Path to input FASTA file (genome assembly)

    Returns:
        int: Total genome size in base pairs (sum of all sequence lengths)

    Example:
        >>> size = get_genome_size("genome.fasta")
        INFO:...Calculating genome size...
        INFO:...Genome size: 3000000000 bp
        >>> print(f"{size:,} bp")
        3,000,000,000 bp

    Note:
        - No progress bar (use get_genome_size_with_progress() for long files)
        - Memory-efficient: sequences not stored, only lengths summed
        - Counts all characters in sequence (including N, gaps if present)
        - Useful for TRF parameter calculation and statistics
        - Logs start and completion messages via logger
    """

    logger.info("Calculating genome size...")
    genome_size = 0
    for _, seq in sc_iter_fasta_brute(fasta_file):
        genome_size += len(seq)
    logger.info(f"Genome size: {genome_size} bp")
    return genome_size


def get_genome_size_with_progress(fasta_file):
    """
    Calculate total genome size from FASTA file.

    Uses the Rust genome-size binary for fast computation (supports .gz).
    Falls back to Python if the binary is not available.

    Args:
        fasta_file (str): Path to input FASTA file (genome assembly)

    Returns:
        int: Total genome size in base pairs (sum of all sequence lengths)
    """
    import shutil
    import subprocess
    import os
    from pathlib import Path

    logger.info(f"Calculating genome size for: {fasta_file}")

    # Try Rust binary first
    genome_size_bin = None
    # Check in satellome bin directory
    bin_dir = Path(__file__).parent.parent.parent / "bin"
    candidate = bin_dir / "genome-size"
    if candidate.exists():
        genome_size_bin = str(candidate)
    else:
        genome_size_bin = shutil.which("genome-size")

    if genome_size_bin:
        try:
            result = subprocess.run(
                [genome_size_bin, fasta_file],
                capture_output=True, text=True, timeout=3600
            )
            if result.returncode == 0:
                for line in result.stdout.strip().split('\n'):
                    parts = line.split('\t')
                    if parts[0] == 'TOTAL' and len(parts) >= 3:
                        genome_size = int(parts[1])
                        n_seqs = int(parts[2])
                        logger.info(f"Total genome size: {genome_size:,} bp in {n_seqs} scaffolds/contigs")
                        return genome_size
        except (subprocess.TimeoutExpired, FileNotFoundError, OSError) as e:
            logger.warning(f"Rust genome-size failed: {e}, falling back to Python")

    # Fallback: Python implementation
    genome_size = 0
    n_seqs = 0
    for _, seq in sc_iter_fasta_brute(fasta_file):
        genome_size += len(seq)
        n_seqs += 1
    logger.info(f"Total genome size: {genome_size:,} bp in {n_seqs} scaffolds/contigs")
    return genome_size


def count_lines_large_file(filename, chunk_size=1024*1024):
    """
    Efficiently count lines in large files using binary chunk reading.

    Fast line counting that avoids loading entire file into memory. Reads
    file in binary mode in fixed-size chunks and counts newline characters.
    Much faster than len(open(file).readlines()) for large files.

    Args:
        filename (str): Path to file to count lines in
        chunk_size (int, optional): Bytes to read per chunk. Defaults to 1MB (1024*1024).

    Returns:
        int: Number of lines in file (newline count)

    Example:
        >>> # Count lines in large TRF output
        >>> count = count_lines_large_file("repeats.tab")
        >>> print(f"{count:,} lines")
        1,523,429 lines
        >>>
        >>> # Use smaller chunks for memory-constrained systems
        >>> count = count_lines_large_file("huge_file.txt", chunk_size=512*1024)

    Performance:
        - Much faster than loading entire file into memory
        - Constant memory usage regardless of file size
        - Typical speed: ~100-500 MB/s depending on disk I/O
        - 1MB chunk size is optimal for most systems

    Note:
        - Opens file in binary mode ('rb') for speed
        - Counts '\\n' bytes (works for Unix/Linux/Mac line endings)
        - May not count final line if file doesn't end with newline
        - Uses walrus operator (Python 3.8+)
        - No progress bar (instant for files < 1GB)
        - Useful for estimating work before processing large datasets
    """
    line_count = 0
    with open(filename, 'rb') as f:
        while chunk := f.read(chunk_size):
            line_count += chunk.count(b'\n')
    return line_count

