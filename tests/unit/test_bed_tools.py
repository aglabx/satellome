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

    # TRF format field indices (18 fields total):
    # 0: project, 1: trf_id, 2: trf_head, 3: trf_l_ind, 4: trf_r_ind,
    # 5: trf_period, 6: trf_n_copy, 7: trf_pmatch, 8: trf_pvar, 9: trf_entropy,
    # 10: trf_consensus, 11: trf_array, 12: trf_array_gc, 13: trf_consensus_gc,
    # 14: trf_array_length, 15: trf_joined, 16: trf_family, 17: trf_ref_annotation
    IDX_PROJECT = 0
    IDX_TRF_ID = 1
    IDX_TRF_HEAD = 2
    IDX_TRF_L_IND = 3
    IDX_TRF_R_IND = 4
    IDX_TRF_PERIOD = 5
    IDX_TRF_N_COPY = 6
    IDX_TRF_PMATCH = 7
    IDX_TRF_PVAR = 8
    IDX_TRF_ENTROPY = 9
    IDX_TRF_CONSENSUS = 10
    IDX_TRF_ARRAY = 11
    IDX_TRF_ARRAY_GC = 12
    IDX_TRF_CONSENSUS_GC = 13
    IDX_TRF_ARRAY_LENGTH = 14
    IDX_TRF_JOINED = 15
    IDX_TRF_FAMILY = 16
    IDX_TRF_REF_ANNOTATION = 17

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
        """Create a simple test BED file (tanbed format: chr start end period score)."""
        bed_file = tmp_path / "test.bed"
        # tanbed BED format: chr  start  end  period  score
        # chr1 is 64bp long (two lines of 32bp)
        # Extract positions 0-10 (should be "ATCGATCGAT")
        content = """chr1\t0\t10\t5\t100
chr1\t10\t20\t5\t100
chr2\t0\t10\t2\t100
chr3\t0\t16\t4\t100
"""
        bed_file.write_text(content)
        return str(bed_file)

    @pytest.fixture
    def test_bed_reverse_strand(self, tmp_path):
        """Create a BED file with reverse strand (6th column)."""
        bed_file = tmp_path / "test_reverse.bed"
        # chr2 positions 0-10: "GGGGGGGGGG"
        # Reverse complement: "CCCCCCCCCC"
        content = """chr2\t0\t10\t2\t100\t-
chr2\t10\t20\t2\t100\t+
"""
        bed_file.write_text(content)
        return str(bed_file)

    @pytest.fixture
    def test_bed_with_comments(self, tmp_path):
        """Create a BED file with comments and empty lines."""
        bed_file = tmp_path / "test_comments.bed"
        content = """# This is a comment
# BED format test file

chr1\t0\t10\t5\t100

# Another comment
chr2\t0\t5\t2\t100
"""
        bed_file.write_text(content)
        return str(bed_file)

    def test_basic_extraction(self, test_fasta, test_bed_simple, tmp_path):
        """Test basic sequence extraction in TRF format."""
        output_file = tmp_path / "output.sat"
        count = extract_sequences_from_bed(test_fasta, test_bed_simple, str(output_file), project="test")

        assert count == 4  # 4 valid BED entries
        assert output_file.exists()

        # Check output content (filter out comments starting with # and header row starting with 'project')
        lines = output_file.read_text().strip().split('\n')
        data_lines = [l for l in lines if l and not l.startswith('#') and not l.startswith('project')]

        assert len(data_lines) == 4

        # First entry: chr1 0-10 should extract "ATCGATCGAT"
        fields = data_lines[0].split('\t')
        assert len(fields) == 18, f"Expected 18 TRF fields, got {len(fields)}"
        assert fields[self.IDX_PROJECT] == "test"
        assert fields[self.IDX_TRF_HEAD] == "chr1"
        assert fields[self.IDX_TRF_L_IND] == "1"  # 1-based
        assert fields[self.IDX_TRF_R_IND] == "10"
        assert fields[self.IDX_TRF_PERIOD] == "5"
        assert fields[self.IDX_TRF_ARRAY_LENGTH] == "10"
        assert fields[self.IDX_TRF_ARRAY] == "ATCGATCGAT"
        assert fields[self.IDX_TRF_CONSENSUS] == "ATCGA"  # First 5 chars (period=5)
        assert fields[self.IDX_TRF_PMATCH] == "50"  # Calculated from sequence (ATCGATCGAT with period 5)

        # Second entry: chr1 10-20 should extract "CGATCGATCG"
        fields = data_lines[1].split('\t')
        assert fields[self.IDX_TRF_ARRAY_LENGTH] == "10"
        assert fields[self.IDX_TRF_ARRAY] == "CGATCGATCG"

    def test_reverse_strand(self, test_fasta, test_bed_reverse_strand, tmp_path):
        """Test reverse strand sequence extraction."""
        output_file = tmp_path / "output_reverse.sat"
        count = extract_sequences_from_bed(test_fasta, test_bed_reverse_strand, str(output_file), project="test")

        assert count == 2

        lines = output_file.read_text().strip().split('\n')
        data_lines = [l for l in lines if l and not l.startswith('#') and not l.startswith('project')]

        # First entry: chr2 0-10 on minus strand
        # Original: "GGGGGGGGGG"
        # Reverse complement: "CCCCCCCCCC"
        fields = data_lines[0].split('\t')
        assert fields[self.IDX_TRF_ARRAY] == "CCCCCCCCCC"

        # Second entry: chr2 10-20 on plus strand
        # Original: "CCCCCCCCCC"
        fields = data_lines[1].split('\t')
        assert fields[self.IDX_TRF_ARRAY] == "CCCCCCCCCC"

    def test_with_comments(self, test_fasta, test_bed_with_comments, tmp_path):
        """Test that comments and empty lines are ignored."""
        output_file = tmp_path / "output_comments.sat"
        count = extract_sequences_from_bed(test_fasta, test_bed_with_comments, str(output_file))

        assert count == 2  # Only 2 valid BED entries

    def test_invalid_coordinates(self, test_fasta, tmp_path):
        """Test handling of invalid coordinates."""
        bed_file = tmp_path / "invalid.bed"
        # chr1 is 64bp, so 100 is out of bounds
        # tanbed format: chr start end period score
        content = """chr1\t100\t110\t5\t100
chr1\t0\t5\t2\t100
chr1\t-1\t10\t5\t100
chr1\t20\t10\t5\t100
"""
        bed_file.write_text(content)

        output_file = tmp_path / "output_invalid.sat"
        count = extract_sequences_from_bed(test_fasta, str(bed_file), str(output_file))

        # Only one valid entry (0-5)
        assert count == 1

    def test_missing_chromosome(self, test_fasta, tmp_path):
        """Test handling of chromosome not in FASTA."""
        bed_file = tmp_path / "missing_chr.bed"
        content = """chrX\t0\t10\t5\t100
chr1\t0\t10\t5\t100
"""
        bed_file.write_text(content)

        output_file = tmp_path / "output_missing.sat"
        count = extract_sequences_from_bed(test_fasta, str(bed_file), str(output_file))

        # Only chr1 should be extracted
        assert count == 1

    def test_bed3_format(self, test_fasta, tmp_path):
        """Test BED3 format (no strand information)."""
        bed_file = tmp_path / "bed3.bed"
        # BED3: only chr, start, end (no period - will default to 1)
        content = """chr1\t0\t10
chr2\t0\t5
"""
        bed_file.write_text(content)

        output_file = tmp_path / "output_bed3.sat"
        count = extract_sequences_from_bed(test_fasta, str(bed_file), str(output_file))

        assert count == 2
        # Should default to + strand and period=1

    def test_uppercase_conversion(self, test_fasta, test_bed_simple, tmp_path):
        """Test that sequences are converted to uppercase."""
        output_file = tmp_path / "output_uppercase.sat"
        count = extract_sequences_from_bed(test_fasta, test_bed_simple, str(output_file))

        lines = output_file.read_text().strip().split('\n')
        data_lines = [l for l in lines if l and not l.startswith('#') and not l.startswith('project')]

        # All sequences should be uppercase
        for line in data_lines:
            fields = line.split('\t')
            seq = fields[self.IDX_TRF_ARRAY]
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
        # tanbed format: chr start end period score
        bed_content = """NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome\t0\t10\t5\t100
"""
        bed_file.write_text(bed_content)

        output_file = tmp_path / "output_spaces.sat"
        count = extract_sequences_from_bed(str(fasta_file), str(bed_file), str(output_file), project="test")

        # Should successfully extract sequence
        assert count == 1

        # Use rstrip('\n') instead of strip() to preserve trailing tabs (empty fields)
        lines = output_file.read_text().rstrip('\n').split('\n')
        data_lines = [l for l in lines if l and not l.startswith('#') and not l.startswith('project')]

        assert len(data_lines) == 1
        fields = data_lines[0].split('\t')
        assert len(fields) == 18, f"Expected 18 TRF fields, got {len(fields)}"

        # IMPORTANT: Output should have only first word in trf_head column
        # Input BED: "NC_000913.3 Escherichia coli str. K-12..."
        # Output: "NC_000913.3"
        assert fields[self.IDX_TRF_HEAD] == "NC_000913.3", f"Expected 'NC_000913.3', got '{fields[self.IDX_TRF_HEAD]}'"

        # Check array length (coordinates: 0-10, length should be 10)
        assert fields[self.IDX_TRF_ARRAY_LENGTH] == "10", f"Expected length=10, got {fields[self.IDX_TRF_ARRAY_LENGTH]}"

        # Should extract "ATCGATCGAT"
        assert fields[self.IDX_TRF_ARRAY] == "ATCGATCGAT"

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
        bed_content = """chr1\t0\t10\t5\t100
"""
        bed_file.write_text(bed_content)

        output_file = tmp_path / "output_duplicate.sat"

        # Should raise ValueError due to duplicate chromosome name
        with pytest.raises(ValueError) as exc_info:
            extract_sequences_from_bed(str(fasta_file), str(bed_file), str(output_file))

        assert "Duplicate chromosome name 'chr1' found in FASTA" in str(exc_info.value)
        assert "First word of FASTA headers must be unique" in str(exc_info.value)

    def test_fasta_output(self, test_fasta, test_bed_simple, tmp_path):
        """Test FASTA output file generation."""
        output_file = tmp_path / "output.sat"
        fasta_output = tmp_path / "output.fasta"
        count = extract_sequences_from_bed(
            test_fasta, test_bed_simple, str(output_file),
            fasta_output_file=str(fasta_output), project="test"
        )

        assert count == 4
        assert fasta_output.exists()

        # Check FASTA content
        fasta_content = fasta_output.read_text()
        lines = fasta_content.strip().split('\n')

        # Should have 8 lines (4 headers + 4 sequences)
        assert len(lines) == 8

        # Check first entry header format: >chr_start_end_length_period
        assert lines[0] == ">chr1_0_10_10_5"
        assert lines[1] == "ATCGATCGAT"
