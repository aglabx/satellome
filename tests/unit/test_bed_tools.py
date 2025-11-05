#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Unit tests for bed_tools.py (SAT-49)

import pytest
import os
import tempfile
from satellome.core_functions.tools.bed_tools import reverse_complement, extract_sequences_from_bed


class TestReverseComplement:
    """Test reverse_complement function."""

    def test_simple_sequence(self):
        """Test reverse complement of simple sequence."""
        assert reverse_complement("ATCG") == "CGAT"
        assert reverse_complement("AAAA") == "TTTT"
        assert reverse_complement("CCCC") == "GGGG"

    def test_mixed_case(self):
        """Test reverse complement handles mixed case."""
        assert reverse_complement("ATC") == "GAT"
        assert reverse_complement("atcg") == "cgat"
        assert reverse_complement("AtCg") == "cGaT"

    def test_with_N(self):
        """Test reverse complement handles N bases."""
        assert reverse_complement("ATNCG") == "CGNAT"
        assert reverse_complement("NNN") == "NNN"

    def test_empty_string(self):
        """Test reverse complement of empty string."""
        assert reverse_complement("") == ""

    def test_single_base(self):
        """Test reverse complement of single base."""
        assert reverse_complement("A") == "T"
        assert reverse_complement("T") == "A"
        assert reverse_complement("C") == "G"
        assert reverse_complement("G") == "C"


class TestExtractSequencesFromBed:
    """Test extract_sequences_from_bed function."""

    @pytest.fixture
    def test_fasta(self, tmp_path):
        """Create a test FASTA file."""
        fasta_file = tmp_path / "test.fasta"
        content = """>chr1 Test chromosome 1
ATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCG
>chr2 Test chromosome 2
GGGGGGGGGGCCCCCCCCCCAAAAAAAAAATTTTTTTTTT
>chr3
ACGTACGTACGTACGT
"""
        fasta_file.write_text(content)
        return str(fasta_file)

    @pytest.fixture
    def test_bed_simple(self, tmp_path):
        """Create a simple test BED file."""
        bed_file = tmp_path / "test.bed"
        # BED format: chr  start  end  name  score  strand
        # chr1 is 64bp long (two lines of 32bp)
        # Extract positions 0-10 (should be "ATCGATCGAT")
        content = """chr1\t0\t10\trepeat1\t100\t+
chr1\t10\t20\trepeat2\t100\t+
chr2\t0\t10\trepeat3\t100\t+
chr3\t0\t16\trepeat4\t100\t+
"""
        bed_file.write_text(content)
        return str(bed_file)

    @pytest.fixture
    def test_bed_reverse_strand(self, tmp_path):
        """Create a BED file with reverse strand."""
        bed_file = tmp_path / "test_reverse.bed"
        # chr1 positions 0-10: "ATCGATCGAT"
        # Reverse complement: "ATCGATCGAT" -> reverse -> "TAGCTAGCTA" -> complement -> "ATCGATCGAT"
        # Wait, that's palindromic! Let me use different coordinates
        # chr2 positions 0-10: "GGGGGGGGGG"
        # Reverse complement: "CCCCCCCCCC"
        content = """chr2\t0\t10\trepeat_minus\t100\t-
chr2\t10\t20\trepeat_plus\t100\t+
"""
        bed_file.write_text(content)
        return str(bed_file)

    @pytest.fixture
    def test_bed_with_comments(self, tmp_path):
        """Create a BED file with comments and empty lines."""
        bed_file = tmp_path / "test_comments.bed"
        content = """# This is a comment
# BED format test file

chr1\t0\t10\trepeat1\t100\t+

# Another comment
chr2\t0\t5\trepeat2\t100\t+
"""
        bed_file.write_text(content)
        return str(bed_file)

    def test_basic_extraction(self, test_fasta, test_bed_simple, tmp_path):
        """Test basic sequence extraction."""
        output_file = tmp_path / "output.trf"
        count = extract_sequences_from_bed(test_fasta, test_bed_simple, str(output_file))

        assert count == 4  # 4 valid BED entries
        assert output_file.exists()

        # Check output content
        lines = output_file.read_text().split('\n')
        # Skip header lines (start with #)
        data_lines = [l for l in lines if l and not l.startswith('#')]

        assert len(data_lines) == 4

        # First entry: chr1 0-10 should extract "ATCGATCGAT"
        fields = data_lines[0].split('\t')
        assert fields[0] == "chr1"
        assert fields[1] == "0"
        assert fields[2] == "10"
        assert fields[-1] == "ATCGATCGAT"

        # Second entry: chr1 10-20 should extract "CGATCGATCG"
        fields = data_lines[1].split('\t')
        assert fields[-1] == "CGATCGATCG"

    def test_reverse_strand(self, test_fasta, test_bed_reverse_strand, tmp_path):
        """Test reverse strand sequence extraction."""
        output_file = tmp_path / "output_reverse.trf"
        count = extract_sequences_from_bed(test_fasta, test_bed_reverse_strand, str(output_file))

        assert count == 2

        lines = output_file.read_text().split('\n')
        data_lines = [l for l in lines if l and not l.startswith('#')]

        # First entry: chr2 0-10 on minus strand
        # Original: "GGGGGGGGGG"
        # Reverse complement: "CCCCCCCCCC"
        fields = data_lines[0].split('\t')
        assert fields[5] == "-"
        assert fields[-1] == "CCCCCCCCCC"

        # Second entry: chr2 10-20 on plus strand
        # Original: "CCCCCCCCCC"
        fields = data_lines[1].split('\t')
        assert fields[5] == "+"
        assert fields[-1] == "CCCCCCCCCC"

    def test_with_comments(self, test_fasta, test_bed_with_comments, tmp_path):
        """Test that comments and empty lines are ignored."""
        output_file = tmp_path / "output_comments.trf"
        count = extract_sequences_from_bed(test_fasta, test_bed_with_comments, str(output_file))

        assert count == 2  # Only 2 valid BED entries

    def test_invalid_coordinates(self, test_fasta, tmp_path):
        """Test handling of invalid coordinates."""
        bed_file = tmp_path / "invalid.bed"
        # chr1 is 64bp, so 100 is out of bounds
        content = """chr1\t100\t110\trepeat1\t100\t+
chr1\t0\t5\trepeat2\t100\t+
chr1\t-1\t10\trepeat3\t100\t+
chr1\t20\t10\trepeat4\t100\t+
"""
        bed_file.write_text(content)

        output_file = tmp_path / "output_invalid.trf"
        count = extract_sequences_from_bed(test_fasta, str(bed_file), str(output_file))

        # Only one valid entry (0-5)
        assert count == 1

    def test_missing_chromosome(self, test_fasta, tmp_path):
        """Test handling of chromosome not in FASTA."""
        bed_file = tmp_path / "missing_chr.bed"
        content = """chrX\t0\t10\trepeat1\t100\t+
chr1\t0\t10\trepeat2\t100\t+
"""
        bed_file.write_text(content)

        output_file = tmp_path / "output_missing.trf"
        count = extract_sequences_from_bed(test_fasta, str(bed_file), str(output_file))

        # Only chr1 should be extracted
        assert count == 1

    def test_bed3_format(self, test_fasta, tmp_path):
        """Test BED3 format (no strand information)."""
        bed_file = tmp_path / "bed3.bed"
        # BED3: only chr, start, end
        content = """chr1\t0\t10
chr2\t0\t5
"""
        bed_file.write_text(content)

        output_file = tmp_path / "output_bed3.trf"
        count = extract_sequences_from_bed(test_fasta, str(bed_file), str(output_file))

        assert count == 2
        # Should default to + strand

    def test_uppercase_conversion(self, test_fasta, test_bed_simple, tmp_path):
        """Test that sequences are converted to uppercase."""
        output_file = tmp_path / "output_uppercase.trf"
        count = extract_sequences_from_bed(test_fasta, test_bed_simple, str(output_file))

        lines = output_file.read_text().split('\n')
        data_lines = [l for l in lines if l and not l.startswith('#')]

        # All sequences should be uppercase
        for line in data_lines:
            seq = line.split('\t')[-1]
            assert seq == seq.upper()
            assert seq.isupper() or len(seq) == 0  # Empty sequence edge case

    def test_chromosome_name_with_spaces(self, tmp_path):
        """Test that chromosome names with spaces are handled correctly."""
        # Create FASTA with full chromosome name in header
        fasta_file = tmp_path / "test_spaces.fasta"
        content = """>NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome
ATCGATCGATCGATCGATCGATCGATCGATCG
"""
        fasta_file.write_text(content)

        # Create BED with full chromosome name (like tanbed outputs)
        bed_file = tmp_path / "test_spaces.bed"
        # BED has full chromosome name in first column
        bed_content = """NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome\t0\t10\trepeat1\t100\t+
"""
        bed_file.write_text(bed_content)

        output_file = tmp_path / "output_spaces.trf"
        count = extract_sequences_from_bed(str(fasta_file), str(bed_file), str(output_file))

        # Should successfully extract sequence
        assert count == 1

        lines = output_file.read_text().split('\n')
        data_lines = [l for l in lines if l and not l.startswith('#')]

        assert len(data_lines) == 1
        fields = data_lines[0].split('\t')
        # Should extract "ATCGATCGAT"
        assert fields[-1] == "ATCGATCGAT"

    def test_duplicate_chromosome_names(self, tmp_path):
        """Test that duplicate chromosome names are detected and raise error."""
        # Create FASTA with duplicate first words in headers
        fasta_file = tmp_path / "test_duplicate.fasta"
        content = """>chr1 First sequence
ATCGATCGATCGATCGATCGATCGATCGATCG
>chr1 Second sequence with same first word
GGGGGGGGGGCCCCCCCCCCAAAAAAAAAATTTTTTTTTT
"""
        fasta_file.write_text(content)

        bed_file = tmp_path / "test_duplicate.bed"
        bed_content = """chr1\t0\t10\trepeat1\t100\t+
"""
        bed_file.write_text(bed_content)

        output_file = tmp_path / "output_duplicate.trf"

        # Should raise ValueError due to duplicate chromosome name
        with pytest.raises(ValueError) as exc_info:
            extract_sequences_from_bed(str(fasta_file), str(bed_file), str(output_file))

        assert "Duplicate chromosome name 'chr1' found in FASTA" in str(exc_info.value)
        assert "First word of FASTA headers must be unique" in str(exc_info.value)
