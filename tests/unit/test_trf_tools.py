"""Unit tests for satellome.core_functions.tools.trf_tools module.

Focus on CRITICAL coordinate restoration logic and filter functions.
"""

import pytest
from satellome.core_functions.tools.trf_tools import (
    restore_coordinates_in_line,
    _filter_by_bottom_array_length,
    _filter_by_bottom_unit_length,
)


class MockTRFObject:
    """Mock TRF object for testing filter functions."""

    def __init__(self, array_length=100, period=2):
        self.trf_array_length = array_length
        self.trf_period = period


class TestRestoreCoordinatesInLine:
    """CRITICAL tests for coordinate restoration from chunk-based TRF output.

    This is the most important function to test as coordinate errors
    would corrupt all downstream analysis.
    """

    def test_restore_basic_coordinates(self):
        """Test basic coordinate restoration from chunk header."""
        # Input: chr10__127750000_127925000 with TRF coordinates 80932 91887
        # Expected: chr10 with absolute coordinates 127830932 127841887
        # (80932 + 127750000 = 127830932, 91887 + 127750000 = 127841887)
        input_line = "chr10__127750000_127925000\t80932\t91887\t2\t50.0\t95\t0\t0\t0\tAT\tATAT\n"
        expected = "chr10\t127830932\t127841887\t2\t50.0\t95\t0\t0\t0\tAT\tATAT\n"

        result = restore_coordinates_in_line(input_line)
        assert result == expected

    def test_restore_coordinates_different_chunks(self):
        """Test coordinate restoration with different chunk offsets."""
        # Chunk starting at 0
        line1 = "chr1__0_1000000\t500\t600\t2\t50.0\t95\t0\t0\t0\tAT\tATAT\n"
        result1 = restore_coordinates_in_line(line1)
        assert result1 == "chr1\t500\t600\t2\t50.0\t95\t0\t0\t0\tAT\tATAT\n"

        # Chunk starting at 1000000
        line2 = "chr1__1000000_2000000\t500\t600\t2\t50.0\t95\t0\t0\t0\tAT\tATAT\n"
        result2 = restore_coordinates_in_line(line2)
        assert result2 == "chr1\t1000500\t1000600\t2\t50.0\t95\t0\t0\t0\tAT\tATAT\n"

    def test_restore_coordinates_large_offset(self):
        """Test coordinate restoration with large chromosome offset."""
        # Very large offset (e.g., end of chromosome)
        line = "chr22__50000000_51000000\t12345\t23456\t10\t100.0\t98\t1\t1\t1\tATGCATGCAT\tseq\n"
        result = restore_coordinates_in_line(line)
        expected = "chr22\t50012345\t50023456\t10\t100.0\t98\t1\t1\t1\tATGCATGCAT\tseq\n"
        assert result == expected

    def test_preserve_non_chunk_header(self):
        """Test that non-chunk headers are preserved unchanged."""
        # Regular header without chunk info
        line = "chr1\t1000\t2000\t2\t50.0\t95\t0\t0\t0\tAT\tATAT\n"
        result = restore_coordinates_in_line(line)
        assert result == line

    def test_preserve_different_format_header(self):
        """Test headers with different formats are preserved."""
        # Header with underscore but not chunk format
        line = "scaffold_123\t1000\t2000\t2\t50.0\t95\t0\t0\t0\tAT\tATAT\n"
        result = restore_coordinates_in_line(line)
        assert result == line

    def test_empty_line(self):
        """Test that empty lines are preserved."""
        result = restore_coordinates_in_line("")
        assert result == ""

        result = restore_coordinates_in_line("   \n")
        assert result == "   \n"

    def test_malformed_line_too_few_fields(self):
        """Test handling of malformed lines with too few fields."""
        # Line with only 2 fields (needs at least 3)
        line = "chr1\t1000\n"
        result = restore_coordinates_in_line(line)
        assert result == line

    def test_malformed_chunk_coordinates(self):
        """Test handling of malformed chunk coordinate format."""
        # Missing second coordinate
        line = "chr1__1000000\t500\t600\t2\t50.0\t95\t0\t0\t0\tAT\tATAT\n"
        result = restore_coordinates_in_line(line)
        # Should return original line if parsing fails
        assert result == line

    def test_invalid_coordinate_values(self):
        """Test handling of non-numeric coordinate values."""
        # Non-numeric TRF coordinates
        line = "chr1__1000000_2000000\tabc\txyz\t2\t50.0\t95\t0\t0\t0\tAT\tATAT\n"
        result = restore_coordinates_in_line(line)
        # Should return original line if parsing fails
        assert result == line

    def test_zero_offset_chunk(self):
        """Test chunk starting at position 0."""
        line = "chr1__0_100000\t500\t600\t2\t50.0\t95\t0\t0\t0\tAT\tATAT\n"
        result = restore_coordinates_in_line(line)
        expected = "chr1\t500\t600\t2\t50.0\t95\t0\t0\t0\tAT\tATAT\n"
        assert result == expected

    def test_preserve_all_fields(self):
        """Test that all fields beyond coordinates are preserved."""
        line = "chr1__1000_2000\t100\t200\t5\t25.5\t88.8\t10\t20\t30\tATGCA\tATGCAATGCA\textra_field\n"
        result = restore_coordinates_in_line(line)
        expected = "chr1\t1100\t1200\t5\t25.5\t88.8\t10\t20\t30\tATGCA\tATGCAATGCA\textra_field\n"
        assert result == expected

    def test_chromosome_name_with_special_chars(self):
        """Test chromosome names with special characters."""
        # Chromosome name with dots and underscores before chunk info
        line = "NC_000001.11__1000_2000\t100\t200\t2\t50.0\t95\t0\t0\t0\tAT\tATAT\n"
        result = restore_coordinates_in_line(line)
        expected = "NC_000001.11\t1100\t1200\t2\t50.0\t95\t0\t0\t0\tAT\tATAT\n"
        assert result == expected

    def test_scaffold_name_restoration(self):
        """Test restoration with scaffold names."""
        line = "scaffold_123__5000_10000\t250\t350\t3\t33.3\t90\t5\t5\t5\tATG\tATGATG\n"
        result = restore_coordinates_in_line(line)
        expected = "scaffold_123\t5250\t5350\t3\t33.3\t90\t5\t5\t5\tATG\tATGATG\n"
        assert result == expected

    def test_boundary_coordinates(self):
        """Test with very small and very large coordinates."""
        # Very small TRF coordinates
        line1 = "chr1__1000_2000\t0\t1\t2\t50.0\t95\t0\t0\t0\tAT\tAT\n"
        result1 = restore_coordinates_in_line(line1)
        expected1 = "chr1\t1000\t1001\t2\t50.0\t95\t0\t0\t0\tAT\tAT\n"
        assert result1 == expected1

        # Large TRF coordinates
        line2 = "chr1__0_1000000\t999998\t999999\t2\t50.0\t95\t0\t0\t0\tAT\tAT\n"
        result2 = restore_coordinates_in_line(line2)
        expected2 = "chr1\t999998\t999999\t2\t50.0\t95\t0\t0\t0\tAT\tAT\n"
        assert result2 == expected2

    def test_multiple_underscores_in_base_name(self):
        """Test chromosome names that contain underscores."""
        # Chromosome name has underscores, then __ chunk separator
        line = "some_long_scaffold_name__1000_2000\t100\t200\t2\t50.0\t95\t0\t0\t0\tAT\tATAT\n"
        result = restore_coordinates_in_line(line)
        expected = "some_long_scaffold_name\t1100\t1200\t2\t50.0\t95\t0\t0\t0\tAT\tATAT\n"
        assert result == expected

    def test_tab_preservation(self):
        """Test that tabs are correctly preserved in output."""
        line = "chr1__1000_2000\t100\t200\t2\t50.0\t95\t0\t0\t0\tAT\tATAT\n"
        result = restore_coordinates_in_line(line)
        # Count tabs in result
        assert result.count('\t') == line.count('\t')

    def test_newline_preservation(self):
        """Test that newline is preserved."""
        line = "chr1__1000_2000\t100\t200\t2\t50.0\t95\t0\t0\t0\tAT\tATAT\n"
        result = restore_coordinates_in_line(line)
        assert result.endswith('\n')

        # Line without newline
        line_no_newline = "chr1__1000_2000\t100\t200\t2\t50.0\t95\t0\t0\t0\tAT\tATAT"
        result_no_newline = restore_coordinates_in_line(line_no_newline)
        assert result_no_newline.endswith('\n')


class TestFilterByBottomArrayLength:
    """Tests for _filter_by_bottom_array_length function."""

    def test_filter_accepts_above_cutoff(self):
        """Test filter accepts TRs above cutoff."""
        obj = MockTRFObject(array_length=1500)
        cutoff = 1000
        assert _filter_by_bottom_array_length(obj, cutoff) is True

    def test_filter_rejects_below_cutoff(self):
        """Test filter rejects TRs below cutoff."""
        obj = MockTRFObject(array_length=500)
        cutoff = 1000
        assert _filter_by_bottom_array_length(obj, cutoff) is False

    def test_filter_rejects_at_cutoff(self):
        """Test filter rejects TRs exactly at cutoff (must be greater)."""
        obj = MockTRFObject(array_length=1000)
        cutoff = 1000
        assert _filter_by_bottom_array_length(obj, cutoff) is False

    def test_filter_boundary_values(self):
        """Test filter at boundary values."""
        cutoff = 1000

        # Just above cutoff
        obj_above = MockTRFObject(array_length=1001)
        assert _filter_by_bottom_array_length(obj_above, cutoff) is True

        # Just below cutoff
        obj_below = MockTRFObject(array_length=999)
        assert _filter_by_bottom_array_length(obj_below, cutoff) is False

    def test_filter_zero_cutoff(self):
        """Test filter with zero cutoff."""
        cutoff = 0

        obj_positive = MockTRFObject(array_length=100)
        assert _filter_by_bottom_array_length(obj_positive, cutoff) is True

        obj_zero = MockTRFObject(array_length=0)
        assert _filter_by_bottom_array_length(obj_zero, cutoff) is False

    def test_filter_large_cutoff(self):
        """Test filter with very large cutoff."""
        cutoff = 1000000

        obj_small = MockTRFObject(array_length=500000)
        assert _filter_by_bottom_array_length(obj_small, cutoff) is False

        obj_large = MockTRFObject(array_length=2000000)
        assert _filter_by_bottom_array_length(obj_large, cutoff) is True

    def test_filter_typical_thresholds(self):
        """Test filter with typical threshold values."""
        # 1kb threshold
        obj = MockTRFObject(array_length=1500)
        assert _filter_by_bottom_array_length(obj, 1000) is True

        # 3kb threshold
        obj = MockTRFObject(array_length=5000)
        assert _filter_by_bottom_array_length(obj, 3000) is True

        # 10kb threshold
        obj = MockTRFObject(array_length=15000)
        assert _filter_by_bottom_array_length(obj, 10000) is True


class TestFilterByBottomUnitLength:
    """Tests for _filter_by_bottom_unit_length function."""

    def test_filter_accepts_above_cutoff(self):
        """Test filter accepts TRs with period above cutoff."""
        obj = MockTRFObject(period=10)
        cutoff = 5
        assert _filter_by_bottom_unit_length(obj, cutoff) is True

    def test_filter_rejects_below_cutoff(self):
        """Test filter rejects TRs with period below cutoff."""
        obj = MockTRFObject(period=3)
        cutoff = 5
        assert _filter_by_bottom_unit_length(obj, cutoff) is False

    def test_filter_rejects_at_cutoff(self):
        """Test filter rejects TRs with period exactly at cutoff."""
        obj = MockTRFObject(period=5)
        cutoff = 5
        assert _filter_by_bottom_unit_length(obj, cutoff) is False

    def test_filter_boundary_values(self):
        """Test filter at boundary values."""
        cutoff = 5

        # Just above
        obj_above = MockTRFObject(period=6)
        assert _filter_by_bottom_unit_length(obj_above, cutoff) is True

        # Just below
        obj_below = MockTRFObject(period=4)
        assert _filter_by_bottom_unit_length(obj_below, cutoff) is False

    def test_filter_microsatellite_threshold(self):
        """Test filter with microsatellite/complex threshold (period=5)."""
        cutoff = 5

        # Microsatellites (period <= 5)
        micro1 = MockTRFObject(period=1)
        micro2 = MockTRFObject(period=2)
        micro5 = MockTRFObject(period=5)

        assert _filter_by_bottom_unit_length(micro1, cutoff) is False
        assert _filter_by_bottom_unit_length(micro2, cutoff) is False
        assert _filter_by_bottom_unit_length(micro5, cutoff) is False

        # Complex (period > 5)
        complex6 = MockTRFObject(period=6)
        complex10 = MockTRFObject(period=10)

        assert _filter_by_bottom_unit_length(complex6, cutoff) is True
        assert _filter_by_bottom_unit_length(complex10, cutoff) is True

    def test_filter_period_ranges(self):
        """Test filter with different period ranges."""
        # Small period cutoff
        obj = MockTRFObject(period=2)
        assert _filter_by_bottom_unit_length(obj, 1) is True
        assert _filter_by_bottom_unit_length(obj, 2) is False

        # Large period cutoff
        obj = MockTRFObject(period=50)
        assert _filter_by_bottom_unit_length(obj, 10) is True
        assert _filter_by_bottom_unit_length(obj, 100) is False


class TestCoordinateRestorationEdgeCases:
    """Additional edge case tests for coordinate restoration."""

    def test_whitespace_in_line(self):
        """Test handling of extra whitespace."""
        # Line with spaces (should still work as split on tabs)
        line = "chr1__1000_2000\t100\t200\t2\t50.0\t95\t0\t0\t0\tAT\tATAT  \n"
        result = restore_coordinates_in_line(line)
        # Should still restore coordinates properly
        assert "chr1\t1100\t1200" in result

    def test_no_newline(self):
        """Test line without trailing newline."""
        line = "chr1__1000_2000\t100\t200\t2\t50.0\t95\t0\t0\t0\tAT\tATAT"
        result = restore_coordinates_in_line(line)
        # Should add newline
        expected = "chr1\t1100\t1200\t2\t50.0\t95\t0\t0\t0\tAT\tATAT\n"
        assert result == expected

    def test_chunk_with_extra_underscores(self):
        """Test chunk coordinate with extra underscores in values."""
        # Chunk coordinates might have extra underscores (edge case)
        line = "chr1__1000_2000_extra\t100\t200\t2\t50.0\t95\t0\t0\t0\tAT\tATAT\n"
        result = restore_coordinates_in_line(line)
        # Should handle first two values
        expected = "chr1\t1100\t1200\t2\t50.0\t95\t0\t0\t0\tAT\tATAT\n"
        assert result == expected

    def test_negative_coordinates(self):
        """Test handling of negative TRF coordinates (should fail gracefully)."""
        line = "chr1__1000_2000\t-100\t200\t2\t50.0\t95\t0\t0\t0\tAT\tATAT\n"
        result = restore_coordinates_in_line(line)
        # Negative start should still work mathematically (1000 + (-100) = 900)
        expected = "chr1\t900\t1200\t2\t50.0\t95\t0\t0\t0\tAT\tATAT\n"
        assert result == expected

    def test_float_coordinates(self):
        """Test handling of float coordinates (should fail during int conversion)."""
        line = "chr1__1000_2000\t100.5\t200.7\t2\t50.0\t95\t0\t0\t0\tAT\tATAT\n"
        result = restore_coordinates_in_line(line)
        # Should return original line if int() conversion fails
        assert result == line

    def test_very_long_chromosome_name(self):
        """Test with very long chromosome names."""
        long_name = "scaffold_" + "x" * 100
        line = f"{long_name}__1000_2000\t100\t200\t2\t50.0\t95\t0\t0\t0\tAT\tATAT\n"
        result = restore_coordinates_in_line(line)
        expected = f"{long_name}\t1100\t1200\t2\t50.0\t95\t0\t0\t0\tAT\tATAT\n"
        assert result == expected

    def test_real_world_example(self):
        """Test with real-world example from documentation."""
        # Example from docstring
        input_line = "chr10__127750000_127925000\t80932\t91887\t2\t50.0\t95\t0\t0\t0\tAT\tATAT\n"
        result = restore_coordinates_in_line(input_line)

        # Verify chromosome name
        assert result.startswith("chr10\t")

        # Verify coordinates are restored
        parts = result.split('\t')
        assert int(parts[1]) == 80932 + 127750000  # 208682932
        assert int(parts[2]) == 91887 + 127750000  # 208693887
