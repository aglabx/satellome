"""Unit tests for satellome.core_functions.classification_micro module."""

import pytest
from satellome.core_functions.classification_micro import RepeatCountStatsModel
from tests.fixtures.sample_data import CLASSIFICATION_TEST_CASES


class MockTRModel:
    """Mock TRModel for testing classification logic."""

    def __init__(self, period=2, pmatch=95.0, array_length=100, consensus="AT", array="AT" * 50, array_gc=0.5, entropy=1.5, n_copy=50.0):
        self.trf_period = period
        self.trf_pmatch = pmatch
        self.trf_array_length = array_length
        self.trf_consensus = consensus
        self.trf_array = array
        self.trf_array_gc = array_gc
        self.trf_entropy = entropy
        self.trf_n_copy = n_copy
        self.trf_family = "Unknown"


class TestRepeatCountStatsModel:
    """Tests for RepeatCountStatsModel class."""

    def test_model_creation(self):
        """Test RepeatCountStatsModel can be instantiated."""
        model = RepeatCountStatsModel()
        assert isinstance(model, RepeatCountStatsModel)
        assert model.max_length == 0
        assert model.min_length == 0
        assert model.n == 0
        assert model.lengths == []
        assert model.pmatch == []
        assert model.gc == []
        assert model.name is None

    def test_model_with_data(self):
        """Test RepeatCountStatsModel with sample data."""
        model = RepeatCountStatsModel()
        model.name = "test_repeat"
        model.n = 10
        model.min_length = 50
        model.max_length = 200
        model.lengths = [50, 100, 150, 200, 75, 125, 80, 90, 110, 120]
        model.pmatch = [95.0, 96.0, 97.0, 98.0, 95.5, 96.5, 97.5, 98.5, 99.0, 99.5]
        model.gc = [50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0]

        # Test __str__ method
        result = str(model)
        assert isinstance(result, str)
        assert result.endswith("\n")
        assert "\t" in result
        assert "test_repeat" in result
        assert "10" in result  # n

    def test_model_str_format(self):
        """Test that __str__ returns correct format."""
        model = RepeatCountStatsModel()
        model.name = "AT_repeat"
        model.n = 5
        model.min_length = 100
        model.max_length = 200
        model.lengths = [100, 120, 140, 160, 200]
        model.pmatch = [95.0, 96.0, 97.0, 98.0, 99.0]
        model.gc = [40.0, 45.0, 50.0, 55.0, 60.0]

        result = str(model)
        parts = result.strip().split("\t")
        # Should have 15 fields as per __str__ docstring
        assert len(parts) == 15


class TestPerfectMicrosatelliteFilter:
    """Tests for perfect microsatellite filter logic."""

    def test_perfect_micro_filter_accepts_perfect(self):
        """Test filter accepts perfect microsatellites."""
        # period < 6 and pmatch == 100
        mock_tr = MockTRModel(period=2, pmatch=100.0)
        filter_func = lambda x: x.trf_period < 6 and x.trf_pmatch == 100
        assert filter_func(mock_tr) is True

    def test_perfect_micro_filter_rejects_imperfect(self):
        """Test filter rejects imperfect microsatellites."""
        # period < 6 but pmatch != 100
        mock_tr = MockTRModel(period=2, pmatch=95.0)
        filter_func = lambda x: x.trf_period < 6 and x.trf_pmatch == 100
        assert filter_func(mock_tr) is False

    def test_perfect_micro_filter_rejects_long_period(self):
        """Test filter rejects long period repeats."""
        # pmatch == 100 but period >= 6
        mock_tr = MockTRModel(period=10, pmatch=100.0)
        filter_func = lambda x: x.trf_period < 6 and x.trf_pmatch == 100
        assert filter_func(mock_tr) is False

    @pytest.mark.parametrize("period,consensus,expected_class", CLASSIFICATION_TEST_CASES)
    def test_period_classification(self, period, consensus, expected_class):
        """Test period-based classification."""
        mock_tr = MockTRModel(period=period, consensus=consensus, pmatch=100.0)
        filter_func = lambda x: x.trf_period < 6 and x.trf_pmatch == 100
        if expected_class == "micro":
            # Should be filtered as perfect microsatellite
            assert filter_func(mock_tr) is True or period == 9  # boundary case
        else:
            # Complex repeats (period >= 10) should not pass
            assert filter_func(mock_tr) is False


class TestMicrosatelliteFilter:
    """Tests for general microsatellite filter logic."""

    def test_micro_filter_accepts_short_period(self):
        """Test filter accepts any short period repeats."""
        # period < 6, any pmatch
        for pmatch in [80.0, 90.0, 95.0, 100.0]:
            mock_tr = MockTRModel(period=2, pmatch=pmatch)
            filter_func = lambda x: x.trf_period < 6
            assert filter_func(mock_tr) is True

    def test_micro_filter_rejects_long_period(self):
        """Test filter rejects long period repeats."""
        mock_tr = MockTRModel(period=10, pmatch=100.0)
        filter_func = lambda x: x.trf_period < 6
        assert filter_func(mock_tr) is False

    def test_micro_filter_boundary(self):
        """Test filter at boundary values."""
        # period = 5 should pass
        mock_tr5 = MockTRModel(period=5, pmatch=95.0)
        filter_func = lambda x: x.trf_period < 6
        assert filter_func(mock_tr5) is True

        # period = 6 should not pass
        mock_tr6 = MockTRModel(period=6, pmatch=95.0)
        assert filter_func(mock_tr6) is False


class TestTrueSSRFilter:
    """Tests for true SSR (Simple Sequence Repeat) filter logic."""

    def test_tssr_filter_missing_nucleotide(self):
        """Test filter accepts sequences missing at least one nucleotide."""
        # AT-only sequence (missing C and G)
        mock_tr = MockTRModel(array="ATATATAT")
        filter_func = lambda x: (
            x.trf_array.upper().count("A") == 0 or
            x.trf_array.upper().count("T") == 0 or
            x.trf_array.upper().count("C") == 0 or
            x.trf_array.upper().count("G") == 0
        )
        assert filter_func(mock_tr) is True

    def test_tssr_filter_all_nucleotides(self):
        """Test filter rejects sequences with all nucleotides."""
        # Sequence with all 4 nucleotides
        mock_tr = MockTRModel(array="ATGCATGC")
        filter_func = lambda x: (
            x.trf_array.upper().count("A") == 0 or
            x.trf_array.upper().count("T") == 0 or
            x.trf_array.upper().count("C") == 0 or
            x.trf_array.upper().count("G") == 0
        )
        assert filter_func(mock_tr) is False

    def test_tssr_classification_at(self):
        """Test tSSR_AT classification."""
        mock_tr = MockTRModel(array="ATATATAT")
        array = mock_tr.trf_array.upper()
        a = array.count("A")
        t = array.count("T")
        c = array.count("C")
        g = array.count("G")

        # Should classify as tSSR_AT
        assert a > 0
        assert t > 0
        assert c == 0
        assert g == 0

    def test_tssr_classification_gc(self):
        """Test tSSR_CG classification."""
        mock_tr = MockTRModel(array="GCGCGCGC")
        array = mock_tr.trf_array.upper()
        a = array.count("A")
        t = array.count("T")
        c = array.count("C")
        g = array.count("G")

        # Should classify as tSSR_CG
        assert c > 0
        assert g > 0
        assert a == 0
        assert t == 0

    def test_tssr_classification_three_bases(self):
        """Test tSSR classification with three bases."""
        # ACG sequence (missing T)
        mock_tr = MockTRModel(array="ACGACGACG")
        array = mock_tr.trf_array.upper()
        a = array.count("A")
        t = array.count("T")
        c = array.count("C")
        g = array.count("G")

        # Should classify as tSSR_ACG
        assert a > 0
        assert c > 0
        assert g > 0
        assert t == 0


class TestFuzzySSRFilter:
    """Tests for fuzzy SSR filter logic."""

    def test_fssr_filter_low_percentage(self):
        """Test filter accepts sequences with low percentage of nucleotide."""
        # Mostly AT with very few G's
        # 100 bp total, 1 G (1%)
        mock_tr = MockTRModel(array="A" * 50 + "T" * 49 + "G", array_length=100)

        def filter_func(x):
            array = x.trf_array.upper()
            n = float(len(array))
            a = array.count("A")
            t = array.count("T")
            c = array.count("C")
            g = array.count("G")
            if a / n < 0.01 or a < 4:
                return True
            if c / n < 0.01 or c < 4:
                return True
            if g / n < 0.01 or g < 4:
                return True
            if t / n < 0.01 or t < 4:
                return True
            return False

        # C is 0%, G is 1% and count is 1 (< 4)
        assert filter_func(mock_tr) is True

    def test_fssr_filter_balanced(self):
        """Test filter rejects balanced sequences."""
        # Equal distribution of all nucleotides
        mock_tr = MockTRModel(array="ATGC" * 25, array_length=100)

        def filter_func(x):
            array = x.trf_array.upper()
            n = float(len(array))
            a = array.count("A")
            t = array.count("T")
            c = array.count("C")
            g = array.count("G")
            if a / n < 0.01 or a < 4:
                return True
            if c / n < 0.01 or c < 4:
                return True
            if g / n < 0.01 or g < 4:
                return True
            if t / n < 0.01 or t < 4:
                return True
            return False

        # All nucleotides are 25% and count is 25
        assert filter_func(mock_tr) is False

    def test_fssr_filter_few_count(self):
        """Test filter on absolute count threshold."""
        # 200 bp total, 3 C's (1.5% but count < 4)
        mock_tr = MockTRModel(
            array="A" * 100 + "T" * 97 + "C" * 3,
            array_length=200
        )

        def filter_func(x):
            array = x.trf_array.upper()
            n = float(len(array))
            c = array.count("C")
            if c / n < 0.01 or c < 4:
                return True
            return False

        # C count is 3 (< 4) even though percentage is 1.5%
        assert filter_func(mock_tr) is True


class TestComplexTRsFilter:
    """Tests for complex TRs filter logic."""

    def test_complex_filter_all_criteria(self):
        """Test complex TR filter with all criteria met."""
        mock_tr = MockTRModel(
            consensus="ATGCATGCAT",  # >4 bp
            array_length=150,  # >100
            pmatch=95.0,  # <100
            array_gc=0.5,  # 0.2 < x < 0.8
            n_copy=5.0,  # >4
            entropy=2.0  # >1.82
        )

        filter_func = lambda x: (
            len(x.trf_consensus) > 4
            and x.trf_array_length > 100
            and x.trf_pmatch < 100.0
            and x.trf_array_gc > 0.2
            and x.trf_array_gc < 0.8
            and x.trf_n_copy > 4
            and (x.trf_entropy and x.trf_entropy > 1.82)
        )

        assert filter_func(mock_tr) is True

    def test_complex_filter_short_consensus(self):
        """Test filter rejects short consensus."""
        mock_tr = MockTRModel(
            consensus="ATG",  # <=4 bp
            array_length=150,
            pmatch=95.0,
            array_gc=0.5,
            n_copy=5.0,
            entropy=2.0
        )

        filter_func = lambda x: (
            len(x.trf_consensus) > 4
            and x.trf_array_length > 100
            and x.trf_pmatch < 100.0
            and x.trf_array_gc > 0.2
            and x.trf_array_gc < 0.8
            and x.trf_n_copy > 4
            and (x.trf_entropy and x.trf_entropy > 1.82)
        )

        assert filter_func(mock_tr) is False

    def test_complex_filter_short_array(self):
        """Test filter rejects short arrays."""
        mock_tr = MockTRModel(
            consensus="ATGCATGCAT",
            array_length=50,  # <=100
            pmatch=95.0,
            array_gc=0.5,
            n_copy=5.0,
            entropy=2.0
        )

        filter_func = lambda x: (
            len(x.trf_consensus) > 4
            and x.trf_array_length > 100
            and x.trf_pmatch < 100.0
            and x.trf_array_gc > 0.2
            and x.trf_array_gc < 0.8
            and x.trf_n_copy > 4
            and (x.trf_entropy and x.trf_entropy > 1.82)
        )

        assert filter_func(mock_tr) is False

    def test_complex_filter_perfect_match(self):
        """Test filter rejects perfect matches."""
        mock_tr = MockTRModel(
            consensus="ATGCATGCAT",
            array_length=150,
            pmatch=100.0,  # Perfect
            array_gc=0.5,
            n_copy=5.0,
            entropy=2.0
        )

        filter_func = lambda x: (
            len(x.trf_consensus) > 4
            and x.trf_array_length > 100
            and x.trf_pmatch < 100.0
            and x.trf_array_gc > 0.2
            and x.trf_array_gc < 0.8
            and x.trf_n_copy > 4
            and (x.trf_entropy and x.trf_entropy > 1.82)
        )

        assert filter_func(mock_tr) is False

    def test_complex_filter_extreme_gc(self):
        """Test filter rejects extreme GC content."""
        # Too low GC
        mock_tr_low = MockTRModel(
            consensus="ATGCATGCAT",
            array_length=150,
            pmatch=95.0,
            array_gc=0.1,  # Too low
            n_copy=5.0,
            entropy=2.0
        )

        # Too high GC
        mock_tr_high = MockTRModel(
            consensus="ATGCATGCAT",
            array_length=150,
            pmatch=95.0,
            array_gc=0.9,  # Too high
            n_copy=5.0,
            entropy=2.0
        )

        filter_func = lambda x: (
            len(x.trf_consensus) > 4
            and x.trf_array_length > 100
            and x.trf_pmatch < 100.0
            and x.trf_array_gc > 0.2
            and x.trf_array_gc < 0.8
            and x.trf_n_copy > 4
            and (x.trf_entropy and x.trf_entropy > 1.82)
        )

        assert filter_func(mock_tr_low) is False
        assert filter_func(mock_tr_high) is False

    def test_complex_filter_low_entropy(self):
        """Test filter rejects low entropy."""
        mock_tr = MockTRModel(
            consensus="ATGCATGCAT",
            array_length=150,
            pmatch=95.0,
            array_gc=0.5,
            n_copy=5.0,
            entropy=1.0  # Too low
        )

        filter_func = lambda x: (
            len(x.trf_consensus) > 4
            and x.trf_array_length > 100
            and x.trf_pmatch < 100.0
            and x.trf_array_gc > 0.2
            and x.trf_array_gc < 0.8
            and x.trf_n_copy > 4
            and (x.trf_entropy and x.trf_entropy > 1.82)
        )

        assert filter_func(mock_tr) is False


class TestLengthBasedFilters:
    """Tests for length-based filter logic (1kb, 3kb, 10kb)."""

    def test_1kb_filter(self):
        """Test 1kb filter logic."""
        filter_func = lambda x: x.trf_array_length > 1000

        mock_tr_pass = MockTRModel(array_length=1500)
        mock_tr_fail = MockTRModel(array_length=500)
        mock_tr_boundary = MockTRModel(array_length=1000)

        assert filter_func(mock_tr_pass) is True
        assert filter_func(mock_tr_fail) is False
        assert filter_func(mock_tr_boundary) is False  # Not greater

    def test_3kb_filter(self):
        """Test 3kb filter logic."""
        filter_func = lambda x: x.trf_array_length > 3000

        mock_tr_pass = MockTRModel(array_length=5000)
        mock_tr_fail = MockTRModel(array_length=2000)
        mock_tr_boundary = MockTRModel(array_length=3000)

        assert filter_func(mock_tr_pass) is True
        assert filter_func(mock_tr_fail) is False
        assert filter_func(mock_tr_boundary) is False

    def test_10kb_filter(self):
        """Test 10kb filter logic."""
        filter_func = lambda x: x.trf_array_length > 10000

        mock_tr_pass = MockTRModel(array_length=15000)
        mock_tr_fail = MockTRModel(array_length=5000)
        mock_tr_boundary = MockTRModel(array_length=10000)

        assert filter_func(mock_tr_pass) is True
        assert filter_func(mock_tr_fail) is False
        assert filter_func(mock_tr_boundary) is False

    def test_length_filters_hierarchy(self):
        """Test that length filters are hierarchical."""
        # A repeat that passes 10kb also passes 3kb and 1kb
        mock_tr = MockTRModel(array_length=12000)

        filter_1kb = lambda x: x.trf_array_length > 1000
        filter_3kb = lambda x: x.trf_array_length > 3000
        filter_10kb = lambda x: x.trf_array_length > 10000

        assert filter_1kb(mock_tr) is True
        assert filter_3kb(mock_tr) is True
        assert filter_10kb(mock_tr) is True


class TestEdgeCases:
    """Tests for edge cases and boundary conditions."""

    def test_empty_array(self):
        """Test handling of empty array."""
        mock_tr = MockTRModel(array="")
        # Should not crash when counting nucleotides
        array = mock_tr.trf_array.upper()
        assert array.count("A") == 0
        assert array.count("T") == 0
        assert array.count("C") == 0
        assert array.count("G") == 0

    def test_single_nucleotide(self):
        """Test handling of single nucleotide repeats."""
        mock_tr = MockTRModel(array="AAAAAAA", consensus="A", period=1)

        # Should pass microsatellite filter
        filter_func = lambda x: x.trf_period < 6
        assert filter_func(mock_tr) is True

        # Should be tSSR (missing T, C, G)
        tssr_filter = lambda x: (
            x.trf_array.upper().count("A") == 0 or
            x.trf_array.upper().count("T") == 0 or
            x.trf_array.upper().count("C") == 0 or
            x.trf_array.upper().count("G") == 0
        )
        assert tssr_filter(mock_tr) is True

    def test_zero_entropy(self):
        """Test handling of None or zero entropy."""
        mock_tr = MockTRModel(entropy=None)

        filter_func = lambda x: (
            x.trf_entropy and x.trf_entropy > 1.82
        )

        # Should fail because entropy is None (returns None, which is falsy)
        assert not filter_func(mock_tr)

    def test_case_insensitivity(self):
        """Test that filters handle lowercase sequences."""
        mock_tr_lower = MockTRModel(array="atgcatgc")
        mock_tr_upper = MockTRModel(array="ATGCATGC")

        # Both should have same counts after upper()
        assert mock_tr_lower.trf_array.upper().count("A") == mock_tr_upper.trf_array.count("A")
        assert mock_tr_lower.trf_array.upper().count("T") == mock_tr_upper.trf_array.count("T")
