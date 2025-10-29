"""Pytest configuration and shared fixtures for Satellome tests."""

import os
import tempfile
import pytest
from pathlib import Path


@pytest.fixture
def temp_output_dir():
    """Create a temporary directory for test outputs.

    The directory is automatically cleaned up after the test.

    Yields:
        Path: Path to temporary directory

    Example:
        def test_something(temp_output_dir):
            output_file = temp_output_dir / "output.txt"
            # ... write to output_file ...
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_sequences():
    """Provide sample DNA sequences for testing.

    Returns:
        dict: Dictionary of test sequences with various characteristics

    Example:
        def test_gc_content(sample_sequences):
            assert get_gc_content(sample_sequences['high_gc']) > 0.6
    """
    return {
        # Basic sequences
        'simple': 'ATGC',
        'empty': '',
        'single': 'A',

        # GC content tests
        'no_gc': 'AAAA',
        'all_gc': 'GGCC',
        'high_gc': 'GCGCGCGC',
        'low_gc': 'ATATATAT',
        'fifty_gc': 'ATGC' * 10,

        # Case sensitivity tests
        'mixed_case': 'AtGc',
        'lowercase': 'atgc',
        'uppercase': 'ATGC',

        # Reverse complement tests
        'palindrome': 'GAATTC',
        'long_seq': 'ATGCATGCATGC',

        # Special characters
        'with_n': 'ATGCNATGC',
        'with_special': 'ATGC~[]',
    }


@pytest.fixture
def sample_trf_data():
    """Provide sample TRF data for testing.

    Returns:
        dict: Dictionary with TRF test data

    Example:
        def test_trf_parsing(sample_trf_data):
            line = sample_trf_data['simple_repeat']
            # ... parse line ...
    """
    return {
        # Simple microsatellite (period=2)
        'simple_repeat': 'chr1\t1000\t2000\t2\t500.0\t85\t10\t20\t5\tAT\t' + 'AT' * 500,

        # Complex repeat (period=10)
        'complex_repeat': 'chr1\t5000\t6000\t10\t100.0\t90\t5\t10\t2\tATGCATGCAT\t' + 'ATGCATGCAT' * 100,

        # Edge case: period=9 (boundary between micro and complex)
        'boundary_repeat': 'chr1\t3000\t4000\t9\t111.0\t88\t8\t15\t3\tATGCATGCA\t' + 'ATGCATGCA' * 111,

        # Chunk coordinates (for restoration testing)
        'chunk_repeat': 'chr1_chunk_0_1000000\t500\t600\t2\t50.0\t85\t10\t20\t5\tAT\t' + 'AT' * 50,

        # Malformed data (for error handling)
        'malformed': 'chr1\tincomplete',
        'empty_line': '',
    }


@pytest.fixture
def sample_fasta_file(temp_output_dir):
    """Create a sample FASTA file for testing.

    Args:
        temp_output_dir: Temporary directory fixture

    Returns:
        Path: Path to the created FASTA file

    Example:
        def test_fasta_parsing(sample_fasta_file):
            sequences = parse_fasta(sample_fasta_file)
            assert len(sequences) == 2
    """
    fasta_file = temp_output_dir / "test.fasta"
    fasta_content = """>chr1 Test chromosome 1
ATGCATGCATGCATGC
GCTAGCTAGCTAGCTA
>chr2 Test chromosome 2
TTTTAAAACCCCGGGG
"""
    fasta_file.write_text(fasta_content)
    return fasta_file


@pytest.fixture
def sample_chunk_headers():
    """Provide sample chunk headers for coordinate restoration testing.

    Returns:
        dict: Dictionary of chunk headers

    Example:
        def test_coordinate_restoration(sample_chunk_headers):
            header = sample_chunk_headers['standard']
            # ... test restoration ...
    """
    return {
        'standard': '>chr1_chunk_0_1000000',
        'zero_offset': '>chr1_chunk_0_0',
        'large_offset': '>chr22_chunk_0_50000000',
        'no_chunk': '>chr1',
        'malformed': '>chr1_chunk',
    }


# Mark all tests in this suite as unit tests by default
def pytest_collection_modifyitems(items):
    """Automatically add 'unit' marker to all tests unless otherwise marked."""
    for item in items:
        if not any(marker in item.keywords for marker in ['integration', 'slow']):
            item.add_marker(pytest.mark.unit)
