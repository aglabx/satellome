"""Test that pytest and fixtures are working correctly."""

import pytest
from pathlib import Path


def test_temp_output_dir(temp_output_dir):
    """Test that temp_output_dir fixture works."""
    assert isinstance(temp_output_dir, Path)
    assert temp_output_dir.exists()
    assert temp_output_dir.is_dir()

    # Test we can write to it
    test_file = temp_output_dir / "test.txt"
    test_file.write_text("test content")
    assert test_file.exists()
    assert test_file.read_text() == "test content"


def test_sample_sequences(sample_sequences):
    """Test that sample_sequences fixture works."""
    assert isinstance(sample_sequences, dict)
    assert 'simple' in sample_sequences
    assert sample_sequences['simple'] == 'ATGC'
    assert sample_sequences['empty'] == ''
    assert sample_sequences['all_gc'] == 'GGCC'


def test_sample_trf_data(sample_trf_data):
    """Test that sample_trf_data fixture works."""
    assert isinstance(sample_trf_data, dict)
    assert 'simple_repeat' in sample_trf_data
    assert 'chr1' in sample_trf_data['simple_repeat']


def test_sample_fasta_file(sample_fasta_file):
    """Test that sample_fasta_file fixture works."""
    assert sample_fasta_file.exists()
    content = sample_fasta_file.read_text()
    assert '>chr1' in content
    assert '>chr2' in content


def test_sample_chunk_headers(sample_chunk_headers):
    """Test that sample_chunk_headers fixture works."""
    assert isinstance(sample_chunk_headers, dict)
    assert 'standard' in sample_chunk_headers
    assert 'chunk' in sample_chunk_headers['standard']


@pytest.mark.unit
def test_pytest_markers():
    """Test that pytest markers are working."""
    # This test should automatically get the 'unit' marker
    pass
