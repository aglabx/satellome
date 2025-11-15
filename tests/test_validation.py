#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 2025-11-15
# @author: Claude Code
# @contact: noreply@anthropic.com

"""
Tests for input validation functions.
"""

import os
import sys
import tempfile
import shutil
import pytest
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from satellome.core_functions.tools.validation import (
    validate_fasta_file,
    validate_gff_file,
    validate_repeatmasker_file,
    validate_trf_binary,
    validate_output_directory,
    FastaValidationError,
    GFFValidationError,
    BinaryValidationError,
    OutputDirValidationError,
    ValidationError,
)


class TestFastaValidation:
    """Test FASTA file validation."""

    def test_valid_fasta_file(self, tmp_path):
        """Test validation of a valid FASTA file."""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(
            ">seq1 test sequence\n"
            "ACGTACGTACGT\n"
            "ACGTACGTACGT\n"
            ">seq2\n"
            "TTTTAAAACCCCGGGG\n"
        )

        stats = validate_fasta_file(str(fasta_file))
        assert stats['num_sequences'] == 2
        assert stats['total_length'] == 40
        assert len(stats['warnings']) == 0

    def test_empty_fasta_file(self, tmp_path):
        """Test validation of empty FASTA file."""
        fasta_file = tmp_path / "empty.fasta"
        fasta_file.write_text("")

        with pytest.raises(FastaValidationError, match="FASTA file is empty"):
            validate_fasta_file(str(fasta_file))

    def test_missing_fasta_file(self):
        """Test validation of non-existent FASTA file."""
        with pytest.raises(FastaValidationError, match="FASTA file not found"):
            validate_fasta_file("/nonexistent/file.fasta")

    def test_fasta_no_headers(self, tmp_path):
        """Test FASTA file without headers."""
        fasta_file = tmp_path / "no_headers.fasta"
        fasta_file.write_text("ACGTACGT\nACGTACGT\n")

        with pytest.raises(FastaValidationError, match="Sequence data before first header"):
            validate_fasta_file(str(fasta_file))

    def test_fasta_sequence_before_header(self, tmp_path):
        """Test FASTA file with sequence before first header."""
        fasta_file = tmp_path / "bad_order.fasta"
        fasta_file.write_text("ACGTACGT\n>seq1\nACGTACGT\n")

        with pytest.raises(FastaValidationError, match="Sequence data before first header"):
            validate_fasta_file(str(fasta_file))

    def test_fasta_invalid_characters(self, tmp_path):
        """Test FASTA file with invalid characters."""
        fasta_file = tmp_path / "invalid_chars.fasta"
        fasta_file.write_text(">seq1\nACGT123XYZ\n")

        stats = validate_fasta_file(str(fasta_file))
        assert len(stats['warnings']) > 0
        assert any('Invalid characters' in w for w in stats['warnings'])

    def test_fasta_empty_sequence(self, tmp_path):
        """Test FASTA file with empty sequence."""
        fasta_file = tmp_path / "empty_seq.fasta"
        fasta_file.write_text(">seq1\n>seq2\nACGT\n")

        stats = validate_fasta_file(str(fasta_file))
        assert any('zero length' in w for w in stats['warnings'])

    def test_fasta_gzipped(self, tmp_path):
        """Test validation of gzipped FASTA file."""
        import gzip
        fasta_file = tmp_path / "test.fasta.gz"
        with gzip.open(fasta_file, 'wt') as f:
            f.write(">seq1\nACGTACGT\n")

        stats = validate_fasta_file(str(fasta_file))
        assert stats['num_sequences'] == 1
        assert stats['total_length'] == 8


class TestGFFValidation:
    """Test GFF file validation."""

    def test_valid_gff_file(self, tmp_path):
        """Test validation of a valid GFF file."""
        gff_file = tmp_path / "test.gff"
        gff_file.write_text(
            "##gff-version 3\n"
            "chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=gene1\n"
            "chr1\ttest\tmRNA\t100\t200\t.\t+\t.\tID=mRNA1;Parent=gene1\n"
        )

        stats = validate_gff_file(str(gff_file))
        assert stats['num_features'] == 2
        assert stats['num_malformed'] == 0

    def test_empty_gff_file(self, tmp_path):
        """Test validation of empty GFF file."""
        gff_file = tmp_path / "empty.gff"
        gff_file.write_text("")

        with pytest.raises(GFFValidationError, match="GFF file is empty"):
            validate_gff_file(str(gff_file))

    def test_missing_gff_file(self):
        """Test validation of non-existent GFF file."""
        with pytest.raises(GFFValidationError, match="GFF file not found"):
            validate_gff_file("/nonexistent/file.gff")

    def test_gff_malformed_lines(self, tmp_path):
        """Test GFF file with malformed lines."""
        gff_file = tmp_path / "malformed.gff"
        gff_file.write_text(
            "##gff-version 3\n"
            "chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=gene1\n"
            "chr1\ttest\tgene\t100\n"  # Missing fields
            "chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=gene2\n"
        )

        stats = validate_gff_file(str(gff_file))
        assert stats['num_features'] == 2
        assert stats['num_malformed'] == 1

    def test_gff_invalid_coordinates(self, tmp_path):
        """Test GFF file with invalid coordinates."""
        gff_file = tmp_path / "invalid_coords.gff"
        gff_file.write_text(
            "##gff-version 3\n"
            "chr1\ttest\tgene\t200\t100\t.\t+\t.\tID=gene1\n"  # start > end
        )

        stats = validate_gff_file(str(gff_file))
        assert any('Start position' in w and '> end position' in w for w in stats['warnings'])


class TestRepeatMaskerValidation:
    """Test RepeatMasker file validation."""

    def test_valid_repeatmasker_file(self, tmp_path):
        """Test validation of a valid RepeatMasker file."""
        rm_file = tmp_path / "test.out"
        # Standard RepeatMasker output format (fields are space-separated)
        # Format: score div del ins chrm start end left strand repeat class begin end left id
        rm_file.write_text(
            "239 23.3 0.0 0.0 chr1 1 100 (1000) C LINE/L1 (100) 500 400 1\n"
            "189 21.1 0.5 0.0 chr1 200 300 (800) + SINE/Alu 1 101 (0) 2\n"
        )

        stats = validate_repeatmasker_file(str(rm_file))
        assert stats['num_features'] == 2
        assert stats['num_malformed'] == 0

    def test_empty_repeatmasker_file(self, tmp_path):
        """Test validation of empty RepeatMasker file."""
        rm_file = tmp_path / "empty.out"
        rm_file.write_text("")

        with pytest.raises(ValidationError, match="RepeatMasker file is empty"):
            validate_repeatmasker_file(str(rm_file))

    def test_missing_repeatmasker_file(self):
        """Test validation of non-existent RepeatMasker file."""
        with pytest.raises(ValidationError, match="RepeatMasker file not found"):
            validate_repeatmasker_file("/nonexistent/file.out")

    def test_repeatmasker_malformed_lines(self, tmp_path):
        """Test RepeatMasker file with malformed lines."""
        rm_file = tmp_path / "malformed.out"
        rm_file.write_text(
            "239 23.3 0.0 0.0 chr1 1 100 (1000) C LINE/L1 (100) 500 400 1\n"
            "189 21.1 0.5 chr1\n"  # Too few fields
        )

        stats = validate_repeatmasker_file(str(rm_file))
        assert stats['num_features'] == 1
        assert stats['num_malformed'] == 1


class TestTRFBinaryValidation:
    """Test TRF binary validation."""

    def test_valid_trf_in_path(self):
        """Test validation of TRF in PATH (if available)."""
        # This test will pass if 'ls' is in PATH (should always be on Unix systems)
        try:
            path = validate_trf_binary("ls")
            assert os.path.exists(path)
            assert os.access(path, os.X_OK)
        except BinaryValidationError:
            pytest.skip("ls not found in PATH")

    def test_invalid_trf_command(self):
        """Test validation of non-existent TRF command."""
        with pytest.raises(BinaryValidationError, match="not found in PATH"):
            validate_trf_binary("nonexistent_trf_command")

    def test_trf_full_path_valid(self, tmp_path):
        """Test validation of TRF with full path."""
        # Create a dummy executable
        trf_path = tmp_path / "trf"
        trf_path.write_text("#!/bin/bash\necho 'fake trf'\n")
        trf_path.chmod(0o755)

        result = validate_trf_binary(str(trf_path))
        assert result == str(trf_path.absolute())

    def test_trf_full_path_not_executable(self, tmp_path):
        """Test validation of TRF with non-executable file."""
        trf_path = tmp_path / "trf"
        trf_path.write_text("fake trf")
        trf_path.chmod(0o644)  # Not executable

        with pytest.raises(BinaryValidationError, match="not executable"):
            validate_trf_binary(str(trf_path))

    def test_trf_full_path_not_found(self):
        """Test validation of TRF with non-existent path."""
        with pytest.raises(BinaryValidationError, match="not found at specified path"):
            validate_trf_binary("/nonexistent/path/to/trf")


class TestOutputDirectoryValidation:
    """Test output directory validation."""

    def test_valid_output_directory(self, tmp_path):
        """Test validation of existing writable directory."""
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        result = validate_output_directory(str(output_dir))
        assert os.path.isabs(result)
        assert os.path.isdir(result)

    def test_create_missing_directory(self, tmp_path):
        """Test creation of missing output directory."""
        output_dir = tmp_path / "new_output"

        result = validate_output_directory(str(output_dir), create_if_missing=True)
        assert os.path.exists(result)
        assert os.path.isdir(result)

    def test_missing_directory_no_create(self, tmp_path):
        """Test validation of missing directory without creation."""
        output_dir = tmp_path / "missing"

        with pytest.raises(OutputDirValidationError, match="does not exist"):
            validate_output_directory(str(output_dir), create_if_missing=False)

    def test_output_is_file_not_directory(self, tmp_path):
        """Test validation when output path is a file, not directory."""
        output_file = tmp_path / "output.txt"
        output_file.write_text("test")

        with pytest.raises(OutputDirValidationError, match="not a directory"):
            validate_output_directory(str(output_file))

    def test_output_directory_no_write_permission(self, tmp_path):
        """Test validation of directory without write permission."""
        output_dir = tmp_path / "readonly"
        output_dir.mkdir()
        output_dir.chmod(0o555)  # Read and execute only, no write

        try:
            with pytest.raises(OutputDirValidationError, match="No write permission"):
                validate_output_directory(str(output_dir))
        finally:
            # Restore permissions for cleanup
            output_dir.chmod(0o755)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
