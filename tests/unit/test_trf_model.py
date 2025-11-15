"""Unit tests for satellome.core_functions.models.trf_model module."""

import pytest
from satellome.core_functions.models.trf_model import (
    clear_sequence,
    TRModel,
    TRsClassificationModel,
)


class TestClearSequence:
    """Tests for clear_sequence() utility function."""

    def test_clear_sequence_uppercase(self):
        """Test that sequence is converted to uppercase."""
        assert clear_sequence("atgc") == "ATGC"
        assert clear_sequence("AtGc") == "ATGC"

    def test_clear_sequence_strip_whitespace(self):
        """Test that whitespace is removed."""
        assert clear_sequence(" ATGC ") == "ATGC"
        assert clear_sequence("AT GC") == "ATGC"
        assert clear_sequence("AT\tGC") == "ATGC"
        assert clear_sequence("AT\nGC") == "ATGC"
        assert clear_sequence("A T G C") == "ATGC"

    def test_clear_sequence_remove_invalid_chars(self):
        """Test that invalid characters are removed."""
        # Numbers should be removed
        assert clear_sequence("A1T2G3C4") == "ATGC"
        # Special chars should be removed
        assert clear_sequence("ATGC!@#$%") == "ATGC"
        # Valid extended alphabet preserved
        assert clear_sequence("ATGCNUWSMKRYBDHV-") == "ATGCNUWSMKRYBDHV-"

    def test_clear_sequence_empty(self):
        """Test empty sequence."""
        assert clear_sequence("") == ""
        assert clear_sequence("   ") == ""

    def test_clear_sequence_combined(self):
        """Test combined cleaning operations."""
        assert clear_sequence(" a t g c 123 ") == "ATGC"
        assert clear_sequence("AtGc\n\tXYZ") == "ATGCY"  # Y is valid IUPAC code (pyrimidine)
        assert clear_sequence("at gc !@# 123") == "ATGC"


class TestTRModelBasic:
    """Basic tests for TRModel class."""

    def test_trmodel_creation(self):
        """Test TRModel can be instantiated."""
        model = TRModel()
        assert isinstance(model, TRModel)

    def test_trmodel_dumpable_attributes(self):
        """Test that dumpable_attributes is defined correctly."""
        model = TRModel()
        assert isinstance(model.dumpable_attributes, list)
        assert len(model.dumpable_attributes) > 0
        assert "project" in model.dumpable_attributes
        assert "trf_id" in model.dumpable_attributes
        assert "trf_period" in model.dumpable_attributes

    def test_trmodel_int_attributes(self):
        """Test that int_attributes is defined correctly."""
        model = TRModel()
        assert isinstance(model.int_attributes, list)
        assert "trf_l_ind" in model.int_attributes
        assert "trf_r_ind" in model.int_attributes
        assert "trf_period" in model.int_attributes

    def test_trmodel_float_attributes(self):
        """Test that float_attributes is defined correctly."""
        model = TRModel()
        assert isinstance(model.float_attributes, list)
        assert "trf_n_copy" in model.float_attributes
        assert "trf_pmatch" in model.float_attributes
        assert "trf_array_gc" in model.float_attributes

    def test_set_project_data(self):
        """Test set_project_data method."""
        model = TRModel()
        model.set_project_data("test_project")
        assert model.project == "test_project"


class TestTRModelInitialization:
    """Tests for TRModel initialization from TRF data."""

    @pytest.fixture
    def sample_trf_head(self):
        """Provide sample TRF header in TRF format."""
        return "Sequence: chr1"

    @pytest.fixture
    def sample_trf_body(self):
        """Provide sample TRF body (empty for now)."""
        return ""

    @pytest.fixture
    def sample_trf_line(self):
        """Provide sample TRF line."""
        # TRF format: start end period copies consensus_size %match %indel score A C G T entropy consensus sequence
        # Coordinates are 1-based inclusive, so 1000-1099 = 100 positions
        return "1000 1099 2 50.0 2 95 5 100 0 50 0 50 1.5 AT " + "AT" * 50

    def test_set_raw_trf_basic(self, sample_trf_head, sample_trf_body, sample_trf_line):
        """Test basic initialization from TRF data."""
        model = TRModel()
        model.set_raw_trf(sample_trf_head, sample_trf_body, sample_trf_line)

        assert model.trf_head == "chr1"
        assert model.trf_chr == "?"  # parse_chromosome_name doesn't recognize "chr1" without context
        assert model.trf_l_ind == 1000
        assert model.trf_r_ind == 1099
        assert model.trf_period == 2
        assert model.trf_n_copy == 50.0
        assert model.trf_pmatch == 95.0
        assert model.trf_pvar == 5  # 100 - 95

    def test_set_raw_trf_sequences(self, sample_trf_head, sample_trf_body, sample_trf_line):
        """Test that sequences are cleaned properly."""
        model = TRModel()
        model.set_raw_trf(sample_trf_head, sample_trf_body, sample_trf_line)

        # Sequences should be uppercase and cleaned
        assert model.trf_consensus == "AT"
        assert model.trf_array == "AT" * 50
        assert len(model.trf_array) == 100
        assert model.trf_array_length == 100

    def test_set_raw_trf_gc_content(self, sample_trf_head, sample_trf_body, sample_trf_line):
        """Test that GC content is calculated."""
        model = TRModel()
        model.set_raw_trf(sample_trf_head, sample_trf_body, sample_trf_line)

        # AT sequence has 0% GC
        assert model.trf_array_gc == 0.0
        assert model.trf_consensus_gc == 0.0

    def test_set_raw_trf_with_gc_sequence(self, sample_trf_head, sample_trf_body):
        """Test GC content calculation with GC-rich sequence."""
        # GC sequence (2000-2099 = 100 positions)
        gc_line = "2000 2099 2 50.0 2 95 5 100 50 0 50 0 1.5 GC " + "GC" * 50
        model = TRModel()
        model.set_raw_trf(sample_trf_head, sample_trf_body, gc_line)

        # GC sequence has 100% GC
        assert model.trf_array_gc == 1.0
        assert model.trf_consensus_gc == 1.0


class TestTRModelOutputMethods:
    """Tests for TRModel output formatting methods."""

    @pytest.fixture
    def simple_model(self):
        """Create a simple TRModel instance with basic data."""
        model = TRModel()
        model.trf_id = 1
        model.project = "test_project"
        model.trf_head = "chr1"  # Set first - needed for trf_gi and trf_chr properties
        model.trf_period = 2
        model.trf_array_length = 100
        model.trf_array_gc = 0.5
        model.trf_pvar = 5
        # trf_gi and trf_chr are now computed properties from trf_head
        model.trf_l_ind = 1000
        model.trf_r_ind = 1100
        model.trf_array = "ATGC" * 25
        model.trf_consensus = "ATGC"
        model.trf_family = "Unknown"
        # Removed fields: trf_repbase, trf_superfamily, trf_subfamily
        return model

    def test_get_string_repr(self, simple_model):
        """Test get_string_repr method."""
        result = simple_model.get_string_repr()
        assert isinstance(result, str)

    def test_get_index_repr(self, simple_model):
        """Test get_index_repr method."""
        result = simple_model.get_index_repr()
        assert isinstance(result, str)
        assert "\t" in result
        assert result.endswith("\n")
        # Check that it contains expected values
        assert "1" in result  # trf_id
        assert "2" in result  # period
        assert "100" in result  # array_length

    def test_get_header_string(self, simple_model):
        """Test get_header_string method."""
        result = simple_model.get_header_string()
        assert isinstance(result, str)
        assert result.startswith("#")
        assert result.endswith("\n")
        assert "\t" in result

    def test_get_numerical_repr(self, simple_model):
        """Test get_numerical_repr method."""
        result = simple_model.get_numerical_repr()
        assert isinstance(result, str)
        assert result.endswith("\n")
        parts = result.strip().split("\t")
        assert len(parts) == 3
        assert parts[0] == "2"  # period
        assert parts[1] == "100"  # array_length

    def test_get_fasta_repr_without_project(self, simple_model):
        """Test get_fasta_repr without project."""
        result = simple_model.get_fasta_repr(add_project=False)
        assert result.startswith(">1\n")
        assert "ATGC" in result
        assert result.endswith("\n")

    def test_get_fasta_repr_with_project(self, simple_model):
        """Test get_fasta_repr with project."""
        result = simple_model.get_fasta_repr(add_project=True)
        assert result.startswith(">1_test_project\n")
        assert "ATGC" in result

    def test_get_monomer_fasta_repr(self, simple_model):
        """Test get_monomer_fasta_repr method."""
        result = simple_model.get_monomer_fasta_repr()
        assert result.startswith(">1\n")
        assert "ATGC\n" in result

    def test_get_family_repr(self, simple_model):
        """Test get_family_repr method."""
        result = simple_model.get_family_repr()
        assert isinstance(result, str)
        assert result.endswith("\n")
        assert "\t" in result
        parts = result.strip().split("\t")
        assert len(parts) == 10  # Updated: now returns 10 fields instead of 13

    def test_fasta_property(self, simple_model):
        """Test fasta property."""
        result = simple_model.fasta
        assert isinstance(result, str)
        assert result.startswith(">")

    def test_get_gff3_string_basic(self, simple_model):
        """Test get_gff3_string with basic parameters."""
        result = simple_model.get_gff3_string()
        assert isinstance(result, str)
        assert result.endswith("\n")
        parts = result.strip().split("\t")
        assert len(parts) == 9
        assert parts[0] == "chr1"  # seqid
        assert parts[1] == "PySatDNA"  # source (tool)
        assert parts[2] == "complex_tandem_repeat"  # type
        assert parts[6] == "+"  # strand

    def test_get_gff3_string_with_prefix(self, simple_model):
        """Test get_gff3_string with prefix."""
        result = simple_model.get_gff3_string(prefix="prefix_")
        assert result.startswith("prefix_chr1\t")

    def test_get_gff3_string_custom_type(self, simple_model):
        """Test get_gff3_string with custom type."""
        result = simple_model.get_gff3_string(trs_type="microsatellite")
        parts = result.strip().split("\t")
        assert parts[2] == "microsatellite"

    def test_get_gff3_string_with_properties(self, simple_model):
        """Test get_gff3_string with custom properties."""
        properties = {"Period": "trf_period", "Length": "trf_array_length"}
        result = simple_model.get_gff3_string(properties=properties)
        assert "Period=2" in result
        assert "Length=100" in result

    def test_get_bed_string(self, simple_model):
        """Test get_bed_string method."""
        result = simple_model.get_bed_string()
        assert isinstance(result, str)
        assert result.endswith("\n")
        parts = result.strip().split("\t")
        assert len(parts) == 3
        assert parts[0] == "chr1"
        assert parts[1] == "1000"
        assert parts[2] == "1100"


class TestTRModelOverlap:
    """Tests for set_form_overlap method."""

    @pytest.fixture
    def first_model(self):
        """Create first TRModel for overlap testing."""
        model = TRModel()
        model.trf_l_ind = 1000
        model.trf_r_ind = 1099  # 100 bp (1-based inclusive: 1000-1099 = 100 positions)
        model.trf_period = 2
        model.trf_n_copy = 50.0
        model.trf_pmatch = 95.0
        model.trf_score = 100.0
        model.trf_array = "AT" * 50  # 100 bp
        model.trf_array_length = 100
        model.trf_consensus = "AT"
        model.trf_entropy = 1.0
        return model

    @pytest.fixture
    def second_model(self):
        """Create second TRModel for overlap testing."""
        model = TRModel()
        model.trf_l_ind = 1050  # Overlaps with first (1050-1099 = 50 bp overlap)
        model.trf_r_ind = 1149  # 100 bp (1050-1149 = 100 positions)
        model.trf_period = 2
        model.trf_n_copy = 50.0
        model.trf_pmatch = 90.0
        model.trf_score = 80.0  # Lower score
        model.trf_array = "GC" * 50  # 100 bp
        model.trf_array_length = 100
        model.trf_consensus = "GC"
        model.trf_entropy = 1.5
        return model

    def test_set_form_overlap_basic(self, first_model, second_model):
        """Test basic overlap joining."""
        first_model.set_form_overlap(second_model)

        # Left coordinate should stay the same
        assert first_model.trf_l_ind == 1000
        # Right coordinate should be from second model
        assert first_model.trf_r_ind == 1149
        # Array length should match coordinates
        assert first_model.trf_array_length == 1149 - 1000 + 1

    def test_set_form_overlap_joined_flag(self, first_model, second_model):
        """Test that joined flag is set."""
        first_model.set_form_overlap(second_model)
        assert first_model.trf_joined == 1

    def test_set_form_overlap_pmatch_weighted(self, first_model, second_model):
        """Test that pmatch is weighted average."""
        original_pmatch1 = first_model.trf_pmatch
        original_pmatch2 = second_model.trf_pmatch
        original_len1 = first_model.trf_array_length
        original_len2 = second_model.trf_array_length

        first_model.set_form_overlap(second_model)

        # pmatch should be weighted average
        expected_pmatch = (original_pmatch1 * original_len1 + original_pmatch2 * original_len2) / (original_len1 + original_len2)
        assert abs(first_model.trf_pmatch - expected_pmatch) < 0.01

    def test_set_form_overlap_higher_score_keeps_period(self, first_model, second_model):
        """Test that higher score model keeps its period."""
        # first_model has higher score (100 vs 80)
        original_period = first_model.trf_period
        first_model.set_form_overlap(second_model)
        assert first_model.trf_period == original_period

    def test_set_form_overlap_array_concatenation(self, first_model, second_model):
        """Test that arrays are concatenated correctly."""
        first_model.set_form_overlap(second_model)
        # Array should be longer than either original
        assert len(first_model.trf_array) > 100


class TestTRsClassificationModel:
    """Tests for TRsClassificationModel class."""

    def test_classification_model_creation(self):
        """Test TRsClassificationModel can be instantiated."""
        model = TRsClassificationModel()
        assert isinstance(model, TRsClassificationModel)

    def test_classification_model_dumpable_attributes(self):
        """Test dumpable_attributes is defined."""
        model = TRsClassificationModel()
        assert isinstance(model.dumpable_attributes, list)
        assert "trf_id" in model.dumpable_attributes
        assert "trf_type" in model.dumpable_attributes
        assert "class_micro" in model.dumpable_attributes

    def test_set_with_trs(self):
        """Test set_with_trs method."""
        trf_model = TRModel()
        trf_model.trf_id = 123
        trf_model.trf_period = 3
        trf_model.trf_array_length = 300
        trf_model.project = "test"

        class_model = TRsClassificationModel()
        class_model.set_with_trs(trf_model)

        assert class_model.trf_id == 123
        assert class_model.trf_period == 3
        assert class_model.trf_array_length == 300
        assert class_model.project == "test"

    def test_network_head_property(self):
        """Test network_head property."""
        model = TRsClassificationModel()
        model.trf_id = 1
        model.trf_period = 2
        model.trf_array_length = 100

        result = model.network_head
        assert isinstance(result, str)
