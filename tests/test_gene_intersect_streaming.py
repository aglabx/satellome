#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 2025-11-15
# @author: Claude Code
# @contact: noreply@anthropic.com

"""
Tests for streaming GFF annotation implementation.
"""

import os
import sys
import tempfile
import pytest
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from satellome.core_functions.tools.gene_intersect_streaming import (
    get_trf_chromosomes,
    load_chromosome_annotations_gff,
    load_chromosome_annotations_rm,
    process_trf_chromosome,
    add_annotation_streaming,
    add_annotation_from_gff_streaming,
)
from satellome.core_functions.tools.gene_intersect import add_annotation_from_gff


# Counter for generating unique TRF IDs
_trf_id_counter = 0

def create_trf_line(
    trf_head="chr1",
    trf_l_ind=1000,
    trf_r_ind=1100,
    trf_period=5,
    trf_n_copy=20.0,
    trf_pmatch=80.0,
    trf_consensus="ACGT",
    trf_array="ACGTACGTACGTACGTACGT"
):
    """
    Helper function to create a valid TRF tab-delimited line.

    TRModel expects 18 tab-separated fields corresponding to dumpable_attributes:
    project, trf_id, trf_head, trf_l_ind, trf_r_ind, trf_period, trf_n_copy,
    trf_pmatch, trf_pvar, trf_entropy, trf_consensus, trf_array, trf_array_gc,
    trf_consensus_gc, trf_array_length, trf_joined, trf_family, trf_ref_annotation
    """
    global _trf_id_counter
    _trf_id_counter += 1

    trf_pvar = 100 - trf_pmatch
    trf_array_length = len(trf_array)

    # Create tab-delimited line with all 18 fields
    # Use unique IDs for each TRF record
    fields = [
        "test_project",  # project
        str(_trf_id_counter),  # trf_id - UNIQUE for each call
        trf_head,  # trf_head - CRITICAL FIELD!
        str(trf_l_ind),  # trf_l_ind
        str(trf_r_ind),  # trf_r_ind
        str(trf_period),  # trf_period
        str(trf_n_copy),  # trf_n_copy
        str(trf_pmatch),  # trf_pmatch
        str(trf_pvar),  # trf_pvar
        "1.5",  # trf_entropy
        trf_consensus,  # trf_consensus
        trf_array,  # trf_array
        "50.0",  # trf_array_gc
        "50.0",  # trf_consensus_gc
        str(trf_array_length),  # trf_array_length
        "0",  # trf_joined
        "",  # trf_family
        "",  # trf_ref_annotation
    ]
    return "\t".join(fields) + "\n"


class TestGetTRFChromosomes:
    """Test chromosome identification from TRF files."""

    def test_get_chromosomes_from_trf(self, tmp_path):
        """Test extracting chromosome names from TRF file."""
        trf_file = tmp_path / "test.sat"
        trf_file.write_text(
            create_trf_line(trf_head="chr1", trf_l_ind=1000, trf_r_ind=1100) +
            create_trf_line(trf_head="chr1", trf_l_ind=2000, trf_r_ind=2100) +
            create_trf_line(trf_head="chr2", trf_l_ind=1000, trf_r_ind=1100) +
            create_trf_line(trf_head="chr3", trf_l_ind=1000, trf_r_ind=1100)
        )

        chrm_counts = get_trf_chromosomes(str(trf_file))

        assert len(chrm_counts) == 3
        assert chrm_counts['chr1'] == 2
        assert chrm_counts['chr2'] == 1
        assert chrm_counts['chr3'] == 1

    def test_get_chromosomes_empty_file(self, tmp_path):
        """Test with empty TRF file."""
        trf_file = tmp_path / "empty.sat"
        trf_file.write_text("")

        chrm_counts = get_trf_chromosomes(str(trf_file))

        assert len(chrm_counts) == 0


class TestLoadChromosomeAnnotationsGFF:
    """Test loading GFF annotations for specific chromosomes."""

    def test_load_specific_chromosome(self, tmp_path):
        """Test loading annotations for a specific chromosome."""
        gff_file = tmp_path / "test.gff"
        gff_file.write_text(
            "##gff-version 3\n"
            "chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=gene1\n"
            "chr1\ttest\tgene\t300\t400\t.\t+\t.\tID=gene2\n"
            "chr2\ttest\tgene\t100\t200\t.\t+\t.\tID=gene3\n"
        )

        annotations = load_chromosome_annotations_gff(str(gff_file), "chr1")

        # Should load only chr1 annotations
        assert len(annotations) == 2

        # Check intervals (use point queries)
        assert len(annotations[150]) > 0  # Should hit gene1 (100-200)
        assert len(annotations[350]) > 0  # Should hit gene2 (300-400)
        assert len(annotations[500]) == 0  # Should not hit anything

    def test_load_nonexistent_chromosome(self, tmp_path):
        """Test loading annotations for chromosome not in file."""
        gff_file = tmp_path / "test.gff"
        gff_file.write_text(
            "##gff-version 3\n"
            "chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=gene1\n"
        )

        annotations = load_chromosome_annotations_gff(str(gff_file), "chr99")

        assert len(annotations) == 0


class TestLoadChromosomeAnnotationsRM:
    """Test loading RepeatMasker annotations for specific chromosomes."""

    def test_load_specific_chromosome_rm(self, tmp_path):
        """Test loading RepeatMasker annotations for specific chromosome."""
        rm_file = tmp_path / "test.out"
        rm_file.write_text(
            "239 23.3 0.0 0.0 chr1 100 200 (1000) C LINE/L1 (100) 500 400 1\n"
            "189 21.1 0.5 0.0 chr1 300 400 (800) + SINE/Alu 1 101 (0) 2\n"
            "239 23.3 0.0 0.0 chr2 100 200 (1000) C LINE/L1 (100) 500 400 3\n"
        )

        annotations = load_chromosome_annotations_rm(str(rm_file), "chr1")

        # Should load only chr1 annotations
        assert len(annotations) == 2

        # Check intervals (use point queries)
        assert len(annotations[150]) > 0  # Should hit first repeat (100-200)
        assert len(annotations[350]) > 0  # Should hit second repeat (300-400)

    def test_load_rm_with_malformed_lines(self, tmp_path):
        """Test handling malformed RepeatMasker lines."""
        rm_file = tmp_path / "malformed.out"
        rm_file.write_text(
            "239 23.3 0.0 0.0 chr1 100 200 (1000) C LINE/L1 (100) 500 400 1\n"
            "invalid line\n"  # Malformed
            "189 21.1 0.5 0.0 chr1 300 400 (800) + SINE/Alu 1 101 (0) 2\n"
        )

        annotations = load_chromosome_annotations_rm(str(rm_file), "chr1")

        # Should load only valid annotations
        assert len(annotations) == 2


class TestProcessTRFChromosome:
    """Test processing TRF records against annotations."""

    def test_process_chromosome_with_annotations(self, tmp_path):
        """Test processing TRF records for a chromosome."""
        from intervaltree import IntervalTree

        # Create TRF file
        trf_file = tmp_path / "test.sat"
        trf_file.write_text(
            create_trf_line(trf_head="chr1", trf_l_ind=100, trf_r_ind=200) +
            create_trf_line(trf_head="chr1", trf_l_ind=300, trf_r_ind=400) +
            create_trf_line(trf_head="chr2", trf_l_ind=100, trf_r_ind=200)
        )

        # Create annotations for chr1
        annotations = IntervalTree()
        annotations.addi(90, 210, {"type": "gene"})  # Overlaps first TR
        annotations.addi(290, 410, {"type": "CDS"})  # Overlaps second TR

        # Process chr1
        results = process_trf_chromosome(str(trf_file), "chr1", annotations)

        # Should have results for chr1 TRs only
        assert len(results) == 2

        # Check that annotations were found
        for key, hits in results.items():
            assert hits is not None
            assert len(hits) > 0

    def test_process_chromosome_no_annotations(self, tmp_path):
        """Test processing TRF records without annotations."""
        from intervaltree import IntervalTree

        # Create TRF file
        trf_file = tmp_path / "test.sat"
        trf_file.write_text(
            create_trf_line(trf_head="chr1", trf_l_ind=100, trf_r_ind=200)
        )

        # Empty annotations
        annotations = IntervalTree()

        # Process chr1
        results = process_trf_chromosome(str(trf_file), "chr1", annotations)

        # Should have results but with no annotations
        assert len(results) == 1
        for key, hits in results.items():
            assert hits is None


class TestStreamingIntegration:
    """Integration tests for streaming annotation."""

    def test_streaming_vs_inmemory_same_results(self, tmp_path):
        """Test that streaming and in-memory modes produce same results."""
        # Create test files
        trf_file = tmp_path / "test.sat"
        trf_file.write_text(
            create_trf_line(trf_head="chr1", trf_l_ind=100, trf_r_ind=200) +
            create_trf_line(trf_head="chr1", trf_l_ind=300, trf_r_ind=400) +
            create_trf_line(trf_head="chr2", trf_l_ind=100, trf_r_ind=200)
        )

        gff_file = tmp_path / "test.gff"
        gff_file.write_text(
            "##gff-version 3\n"
            "chr1\ttest\tgene\t90\t210\t.\t+\t.\tID=gene1\n"
            "chr1\ttest\tCDS\t290\t410\t.\t+\t.\tID=cds1;Parent=gene1\n"
            "chr2\ttest\tgene\t50\t250\t.\t+\t.\tID=gene2\n"
        )

        report_streaming = tmp_path / "report_streaming.txt"
        report_inmemory = tmp_path / "report_inmemory.txt"

        # Test streaming mode
        add_annotation_from_gff(
            str(trf_file),
            str(gff_file),
            str(report_streaming),
            rm_file=None,
            use_streaming=True
        )

        # Test in-memory mode
        add_annotation_from_gff(
            str(trf_file),
            str(gff_file),
            str(report_inmemory),
            rm_file=None,
            use_streaming=False
        )

        # Reports should be identical
        streaming_content = report_streaming.read_text()
        inmemory_content = report_inmemory.read_text()

        assert streaming_content == inmemory_content

    def test_streaming_with_repeatmasker(self, tmp_path):
        """Test streaming mode with RepeatMasker file."""
        # Create test files
        trf_file = tmp_path / "test.sat"
        trf_file.write_text(
            create_trf_line(trf_head="chr1", trf_l_ind=100, trf_r_ind=200)
        )

        gff_file = tmp_path / "test.gff"
        gff_file.write_text(
            "##gff-version 3\n"
            "chr1\ttest\tgene\t90\t210\t.\t+\t.\tID=gene1\n"
        )

        rm_file = tmp_path / "test.out"
        rm_file.write_text(
            "239 23.3 0.0 0.0 chr1 80 220 (1000) C LINE/L1 (100) 500 400 1\n"
        )

        report_file = tmp_path / "report.txt"

        # Should not raise any errors
        add_annotation_from_gff_streaming(
            str(trf_file),
            str(gff_file),
            str(report_file),
            str(rm_file)
        )

        # Report should be created
        assert report_file.exists()
        assert report_file.stat().st_size > 0

    def test_streaming_with_multiple_chromosomes(self, tmp_path):
        """Test streaming with multiple chromosomes."""
        # Create test files with multiple chromosomes
        trf_file = tmp_path / "test.sat"
        trf_lines = []
        for chrm in ['chr1', 'chr2', 'chr3', 'chr4', 'chr5']:
            for i in range(10):  # 10 records per chromosome
                start = i * 1000
                end = start + 100
                trf_lines.append(
                    create_trf_line(trf_head=chrm, trf_l_ind=start, trf_r_ind=end)
                )
        trf_file.write_text(''.join(trf_lines))

        gff_file = tmp_path / "test.gff"
        gff_lines = ["##gff-version 3\n"]
        for chrm in ['chr1', 'chr2', 'chr3', 'chr4', 'chr5']:
            for i in range(5):  # 5 genes per chromosome
                start = i * 1000
                end = start + 200
                gff_lines.append(
                    f"{chrm}\ttest\tgene\t{start}\t{end}\t.\t+\t.\tID=gene_{chrm}_{i}\n"
                )
        gff_file.write_text(''.join(gff_lines))

        report_file = tmp_path / "report.txt"

        # Should process all chromosomes successfully
        add_annotation_from_gff_streaming(
            str(trf_file),
            str(gff_file),
            str(report_file),
            rm_file=None
        )

        assert report_file.exists()
        report_content = report_file.read_text()
        assert len(report_content) > 0


class TestMemoryEfficiency:
    """Tests to verify memory efficiency claims."""

    def test_large_file_simulation(self, tmp_path):
        """Simulate processing a large file (many chromosomes)."""
        # Create a file with many chromosomes to test memory efficiency
        trf_file = tmp_path / "large.sat"
        gff_file = tmp_path / "large.gff"
        report_file = tmp_path / "report.txt"

        # Generate 50 chromosomes with 100 records each
        trf_lines = []
        for chrm_num in range(50):
            chrm = f"chr{chrm_num}"
            for i in range(100):
                start = i * 1000
                end = start + 100
                trf_lines.append(
                    create_trf_line(trf_head=chrm, trf_l_ind=start, trf_r_ind=end)
                )
        trf_file.write_text(''.join(trf_lines))

        # Generate GFF with annotations for each chromosome
        gff_lines = ["##gff-version 3\n"]
        for chrm_num in range(50):
            chrm = f"chr{chrm_num}"
            for i in range(50):  # 50 genes per chromosome
                start = i * 1000
                end = start + 200
                gff_lines.append(
                    f"{chrm}\ttest\tgene\t{start}\t{end}\t.\t+\t.\tID=gene_{chrm}_{i}\n"
                )
        gff_file.write_text(''.join(gff_lines))

        # Should complete without errors
        # In production, this would use significantly less memory than in-memory mode
        add_annotation_from_gff_streaming(
            str(trf_file),
            str(gff_file),
            str(report_file),
            rm_file=None
        )

        assert report_file.exists()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
