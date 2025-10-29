"""Unit tests for satellome.core_functions.tools.parsers module."""

import pytest
from satellome.core_functions.tools.parsers import (
    parse_fasta_head,
    parse_chromosome_name,
    trf_parse_line,
    trf_parse_param,
    trf_parse_head,
    get_wgs_prefix_from_ref,
    get_wgs_prefix_from_head,
    refine_name,
)
from tests.fixtures.sample_data import FASTA_HEADERS, CHROMOSOME_NAMES


class TestParseFastaHead:
    """Tests for parse_fasta_head() function."""

    def test_parse_ncbi_format(self):
        """Test parsing NCBI format FASTA headers."""
        # >gi|20928815|ref|NW_003237.1|MmUn_WIFeb01_12612 Mus musculus chromosome Un genomic contig
        header = ">gi|20928815|ref|NW_003237.1|MmUn_WIFeb01_12612 Mus musculus chromosome Un genomic contig"
        result = parse_fasta_head(header)
        assert len(result) == 3
        assert result[0] == "20928815"
        assert result[1] == "NW_003237.1"
        assert "Mus musculus" in result[2]

    def test_parse_ncbi_dbj_format(self):
        """Test parsing NCBI DBJ format."""
        # >gi|293886233|dbj|BABO01423189.1|
        header = ">gi|293886233|dbj|BABO01423189.1|"
        result = parse_fasta_head(header)
        assert len(result) == 3
        assert result[0] == "293886233"
        assert result[1] == "BABO01423189.1"

    def test_parse_lcl_format(self):
        """Test parsing local (lcl|) format."""
        # >lcl|HmaUn_WGA106_1 Hydra magnipapillata genomic contig
        header = ">lcl|HmaUn_WGA106_1 Hydra magnipapillata genomic contig"
        result = parse_fasta_head(header)
        assert len(result) == 3
        assert result[0] == "HmaUn_WGA106_1"
        assert result[1] == "HmaUn_WGA106_1"
        assert "Hydra" in result[2]

    def test_parse_psu_format(self):
        """Test parsing psu| format."""
        header = ">psu|scaffold_123 Description here"
        result = parse_fasta_head(header)
        assert len(result) == 3
        assert result[0] == "scaffold_123"
        assert result[1] == "scaffold_123"

    def test_parse_short_format(self):
        """Test parsing short format (>123\tSOMETHING)."""
        header = ">134124\tSOMETHING"
        result = parse_fasta_head(header)
        assert len(result) == 3
        assert result[0] == "134124"
        assert result[1] == "SOMETHING"
        assert result[2] == "Unknown"

    def test_parse_composite_format(self):
        """Test parsing composite format (>123-456)."""
        header = ">134124-14124"
        result = parse_fasta_head(header)
        assert len(result) == 3
        assert result[0] == "134124"
        assert result[1] == "14124"
        assert result[2] == "Unknown"

    def test_parse_number_only(self):
        """Test parsing number-only format."""
        header = ">134124"
        result = parse_fasta_head(header)
        assert len(result) == 3
        assert result[0] == "134124"
        assert result[1] == "Unknown"
        assert result[2] == "Unknown"

    def test_parse_probe_index_format(self):
        """Test parsing probe index format."""
        header = ">probe|misat|ref|CAAA01154094|start|991|end|1019"
        result = parse_fasta_head(header)
        assert len(result) == 3
        # gi should be ref_start_end
        assert "CAAA01154094" in result[0]
        assert "991" in result[0]
        assert "1019" in result[0]

    def test_parse_trace_format(self):
        """Test parsing trace format."""
        header = ">gnl|ti|123456789 Some description"
        result = parse_fasta_head(header)
        assert len(result) == 3
        assert result[0] == "123456789"

    def test_parse_unknown_format(self):
        """Test parsing unknown/unsupported format."""
        header = ">SomeWeirdFormat"
        result = parse_fasta_head(header)
        assert len(result) == 3
        assert result[0] == "Unknown"
        assert result[1] == "Unknown"
        assert result[2] == "Unknown"

    def test_parse_without_leading_gt(self):
        """Test parsing headers without leading >."""
        # Function should handle both with and without >
        header_with = ">gi|20928815|ref|NW_003237.1|desc"
        header_without = "gi|20928815|ref|NW_003237.1|desc"
        result_with = parse_fasta_head(header_with)
        result_without = parse_fasta_head(header_without)
        assert result_with[0] == result_without[0]
        assert result_with[1] == result_without[1]

    @pytest.mark.parametrize("header_key", FASTA_HEADERS.keys())
    def test_parse_sample_headers(self, header_key):
        """Test parsing various sample headers from fixtures."""
        header = FASTA_HEADERS[header_key]
        result = parse_fasta_head(header)
        # Should always return list of 3 elements
        assert isinstance(result, list)
        assert len(result) == 3


class TestParseChromosomeName:
    """Tests for parse_chromosome_name() function."""

    def test_parse_standard_chromosome(self):
        """Test parsing standard 'chromosome X' format."""
        header = ">gi|123|ref|NC_000001.11| Homo sapiens chromosome 1, GRCh38.p14 Primary Assembly"
        result = parse_chromosome_name(header)
        assert result == "1"

    def test_parse_chromosome_with_comma(self):
        """Test parsing chromosome name followed by comma."""
        header = "Mus musculus chromosome 10, complete sequence"
        result = parse_chromosome_name(header)
        assert result == "10"

    def test_parse_chr_prefix(self):
        """Test parsing chr prefix format."""
        header = "chr22 Homo sapiens chromosome 22"
        result = parse_chromosome_name(header)
        assert result == "22"

    def test_parse_sex_chromosomes(self):
        """Test parsing X and Y chromosomes."""
        header_x = "chromosome X, complete sequence"
        header_y = "chromosome Y, complete sequence"
        assert parse_chromosome_name(header_x) == "X"
        assert parse_chromosome_name(header_y) == "Y"

    def test_parse_mitochondrial(self):
        """Test parsing mitochondrial chromosome."""
        header1 = "Homo sapiens mitochondrion, complete genome"
        header2 = "complete mitochondrial genome"
        assert parse_chromosome_name(header1) == "MT"
        assert parse_chromosome_name(header2) == "MT"

    def test_parse_no_chromosome_info(self):
        """Test parsing headers without chromosome information."""
        header = "scaffold_12345 genomic contig"
        result = parse_chromosome_name(header)
        assert result == "?"

    def test_parse_empty_string(self):
        """Test parsing empty string."""
        result = parse_chromosome_name("")
        assert result == "?"

    def test_parse_none(self):
        """Test parsing None value."""
        result = parse_chromosome_name(None)
        assert result == "?"

    @pytest.mark.parametrize("name_type,names", CHROMOSOME_NAMES.items())
    def test_parse_sample_chromosome_names(self, name_type, names):
        """Test parsing various chromosome name formats from fixtures."""
        for name in names:
            # Should either parse successfully or return "?"
            result = parse_chromosome_name(name)
            assert isinstance(result, str)
            assert len(result) > 0


class TestTRFParseLine:
    """Tests for trf_parse_line() function."""

    def test_parse_valid_trf_line(self):
        """Test parsing valid TRF line with 15 fields."""
        line = "1000 1100 2 50.0 2 95 5 0 50 0 50 0 1.5 AT ATAT"
        result = trf_parse_line(line)
        assert len(result) == 15
        assert result[0] == "1000"
        assert result[1] == "1100"
        assert result[2] == "2"
        assert result[13] == "AT"
        assert result[14] == "ATAT"

    def test_parse_trf_line_with_tabs(self):
        """Test parsing TRF line with tabs."""
        line = "1000\t1100\t2\t50.0\t2\t95\t5\t0\t50\t0\t50\t0\t1.5\tAT\tATAT"
        result = trf_parse_line(line)
        assert len(result) == 15

    def test_parse_malformed_trf_line(self):
        """Test parsing malformed TRF line (wrong number of fields)."""
        line = "1000 1100 2"  # Only 3 fields
        result = trf_parse_line(line)
        # Should return default values
        assert len(result) == 15
        assert result[0] == 0
        assert result[13] == "0"
        assert result[14] == "0"

    def test_parse_empty_line(self):
        """Test parsing empty line."""
        result = trf_parse_line("")
        assert len(result) == 15
        assert result[0] == 0

    def test_parse_trf_line_with_whitespace(self):
        """Test parsing TRF line with extra whitespace."""
        # Note: Multiple spaces cause split issues, returns default values
        line = "  1000   1100   2   50.0   2   95   5   0   50   0   50   0   1.5   AT   ATAT  "
        result = trf_parse_line(line)
        assert len(result) == 15
        # Returns default due to incorrect field count from multiple spaces
        assert result[0] == 0

    def test_parse_trf_line_with_long_sequences(self):
        """Test parsing TRF line with long sequence strings."""
        consensus = "ATGCATGCAT"
        array = consensus * 50
        line = f"5000 6000 10 50.0 10 90 5 10 100 100 100 100 2.0 {consensus} {array}"
        result = trf_parse_line(line)
        assert len(result) == 15
        assert result[13] == consensus
        assert result[14] == array


class TestTRFParseParam:
    """Tests for trf_parse_param() function."""

    def test_parse_valid_param_line(self):
        """Test parsing valid TRF parameters line."""
        line = "Parameters: 2 7 7 80 10 50 500"
        result = trf_parse_param(line)
        assert "2 7 7 80 10 50 500" in result

    def test_parse_param_with_newline(self):
        """Test parsing parameters with newlines."""
        line = "Parameters: 2 7 7 80 10 50 500\nSome other text"
        result = trf_parse_param(line)
        assert "2 7 7 80 10 50 500" in result

    def test_parse_param_multiline(self):
        """Test parsing parameters in multiline string."""
        line = "Some header\nParameters: 2 5 7 80 10 50 2000\nMore data"
        result = trf_parse_param(line)
        assert "2 5 7 80 10 50 2000" in result

    def test_parse_invalid_param_line(self):
        """Test parsing line without parameters."""
        line = "No parameters here"
        result = trf_parse_param(line)
        assert result == "Unknown"

    def test_parse_empty_param_line(self):
        """Test parsing empty line."""
        result = trf_parse_param("")
        assert result == "Unknown"


class TestTRFParseHead:
    """Tests for trf_parse_head() function."""

    def test_parse_valid_head_with_newline(self):
        """Test parsing valid TRF head with newline."""
        line = "Sequence: chr1\nSome other data"
        result = trf_parse_head(line)
        assert result == "chr1"

    def test_parse_valid_head_without_newline(self):
        """Test parsing valid TRF head without newline."""
        line = "Sequence: chr2"
        result = trf_parse_head(line)
        assert result == "chr2"

    def test_parse_head_with_long_description(self):
        """Test parsing TRF head with long description."""
        line = "Sequence: NC_000001.11 Homo sapiens chromosome 1\nData..."
        result = trf_parse_head(line)
        assert "NC_000001.11" in result

    def test_parse_invalid_head(self):
        """Test parsing line without Sequence:."""
        line = "No sequence here"
        result = trf_parse_head(line)
        # Function returns None when regex doesn't match
        assert result is None

    def test_parse_empty_head(self):
        """Test parsing empty line."""
        result = trf_parse_head("")
        # Function returns None for empty string
        assert result is None

    def test_parse_head_multiline(self):
        """Test parsing head in multiline string."""
        line = "Some header\nSequence: scaffold_123\nParameters: ..."
        result = trf_parse_head(line)
        assert result == "scaffold_123"


class TestGetWGSPrefixFromRef:
    """Tests for get_wgs_prefix_from_ref() function."""

    def test_extract_prefix_uppercase(self):
        """Test extracting WGS prefix from uppercase ref."""
        assert get_wgs_prefix_from_ref("AADD01") == "AADD"
        assert get_wgs_prefix_from_ref("CAAA01154094") == "CAAA"

    def test_extract_prefix_mixed_case(self):
        """Test extracting WGS prefix with mixed case."""
        # Function looks for uppercase letters
        assert get_wgs_prefix_from_ref("AbCd01") == "A"  # Only first uppercase

    def test_extract_prefix_no_letters(self):
        """Test handling ref without letters."""
        assert get_wgs_prefix_from_ref("012345") == "UNKN"

    def test_extract_prefix_empty_string(self):
        """Test handling empty string."""
        assert get_wgs_prefix_from_ref("") == "UNKN"


class TestGetWGSPrefixFromHead:
    """Tests for get_wgs_prefix_from_head() function."""

    def test_extract_from_ref_format(self):
        """Test extracting WGS prefix from ref| format."""
        header = ">gi|123|ref|AADD01123456|description"
        result = get_wgs_prefix_from_head(header)
        assert result == "AADD"

    def test_extract_from_gb_format(self):
        """Test extracting WGS prefix from gb| format."""
        header = ">gi|123|gb|CAAA01123456|description"
        result = get_wgs_prefix_from_head(header)
        assert result == "CAAA"

    def test_extract_no_prefix(self):
        """Test handling header without WGS prefix."""
        header = ">chr1 Some description"
        result = get_wgs_prefix_from_head(header)
        assert result is None

    def test_extract_short_prefix(self):
        """Test that short prefixes (<4 chars) are not matched."""
        header = ">gi|123|ref|AAD|description"
        result = get_wgs_prefix_from_head(header)
        assert result is None

    def test_extract_from_complex_header(self):
        """Test extracting from complex header."""
        header = ">gi|293886233|ref|BABO01423189.1| Homo sapiens"
        result = get_wgs_prefix_from_head(header)
        assert result == "BABO"


class MockTRFObject:
    """Mock TRF object for testing refine_name."""

    def __init__(self):
        self.trf_head = "chr1 description"
        self.trf_l_ind = 1000
        self.trf_r_ind = 2000
        self.trf_consensus = "atgc"
        self.trf_array = "atgcatgc"
        self.trf_id = None
        self.id = None


class TestRefineName:
    """Tests for refine_name() function."""

    def test_refine_basic_name(self):
        """Test basic name refinement."""
        obj = MockTRFObject()
        result = refine_name(0, obj)

        assert result.trf_id == "chr1_1000_2000"
        assert result.id == "AGT0000000000100"
        assert result.trf_consensus == "ATGC"  # Uppercased
        assert result.trf_array == "ATGCATGC"  # Uppercased

    def test_refine_with_index(self):
        """Test name refinement with different index."""
        obj = MockTRFObject()
        result = refine_name(99, obj)

        assert result.trf_id == "chr1_1000_2000"
        assert result.id == "AGT0000000010000"  # (99+1) * 100

    def test_refine_long_header(self):
        """Test name refinement with long header."""
        obj = MockTRFObject()
        obj.trf_head = "scaffold_12345_extra_info more details"
        result = refine_name(0, obj)

        # Should use only first word
        assert result.trf_id == "scaffold_12345_extra_info_1000_2000"

    def test_refine_empty_header(self):
        """Test name refinement with empty header."""
        obj = MockTRFObject()
        obj.trf_head = ""
        result = refine_name(0, obj)

        # Empty string split() returns empty list, which becomes "[]" in string
        assert result.trf_id == "[]_1000_2000"

    def test_uppercase_sequences(self):
        """Test that sequences are uppercased."""
        obj = MockTRFObject()
        obj.trf_consensus = "atgcNn"
        obj.trf_array = "atgcatgcNN"

        result = refine_name(0, obj)

        assert result.trf_consensus == "ATGCNN"
        assert result.trf_array == "ATGCATGCNN"


class TestEdgeCases:
    """Tests for edge cases and error handling."""

    def test_parse_fasta_head_with_special_chars(self):
        """Test parsing FASTA headers with special characters."""
        header = ">gi|123|ref|NC_000001.11| Test & Special <chars>"
        result = parse_fasta_head(header)
        assert len(result) == 3

    def test_parse_chromosome_with_unicode(self):
        """Test parsing chromosome name with unicode."""
        # Should handle gracefully
        header = "chromosome α beta"
        result = parse_chromosome_name(header)
        # Should not crash

    def test_trf_parse_with_negative_numbers(self):
        """Test TRF parsing with negative numbers."""
        line = "-1000 -1100 2 50.0 2 95 5 0 50 0 50 0 1.5 AT ATAT"
        result = trf_parse_line(line)
        assert len(result) == 15
        assert result[0] == "-1000"

    def test_trf_parse_with_floats(self):
        """Test TRF parsing with float coordinates."""
        line = "1000.5 1100.7 2 50.0 2 95 5 0 50 0 50 0 1.5 AT ATAT"
        result = trf_parse_line(line)
        assert len(result) == 15
        assert "1000.5" in result[0]
