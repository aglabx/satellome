"""Tests for per-step wall-clock accounting."""

import time

import pytest

from satellome.core_functions.tools.step_timing import (
    StepTimer,
    format_duration,
    format_timing_table,
)


class TestFormatDuration:
    @pytest.mark.parametrize("seconds,expected", [
        (0.0, "0.0s"), (1.25, "1.2s"), (59.9, "59.9s"),
        (60, "1m00s"), (95, "1m35s"), (3599, "59m59s"),
        (3600, "1h00m"), (7845, "2h10m"),
    ])
    def test_readable_at_every_scale(self, seconds, expected):
        assert format_duration(seconds) == expected


class TestStepTimer:
    def test_records_steps_in_order_of_execution(self):
        timer = StepTimer()
        for name in ("fastan", "classification", "drawing"):
            with timer.step(name):
                pass
        assert list(timer.as_dict()) == ["fastan", "classification", "drawing"]

    def test_a_step_that_raises_is_still_timed(self):
        """How long a failure took is often the interesting part."""
        timer = StepTimer()
        with pytest.raises(RuntimeError):
            with timer.step("drawing"):
                raise RuntimeError("boom")
        assert "drawing" in timer.as_dict()

    def test_repeated_step_accumulates_rather_than_overwrites(self):
        timer = StepTimer()
        timer.record("classification", 2.0)
        timer.record("classification", 3.0)
        assert timer.as_dict()["classification"] == 5.0

    def test_rows_are_sorted_slowest_first(self):
        timer = StepTimer()
        timer.record("drawing", 5.0)
        timer.record("fastan", 100.0)
        timer.record("annotations", 20.0)
        assert [r[0] for r in timer.rows()] == ["fastan", "annotations", "drawing"]

    def test_percentages_are_of_wall_clock_not_of_the_step_sum(self):
        """Otherwise untimed work (validation, genome size) would vanish."""
        timer = StepTimer()
        timer._started = time.time() - 100.0
        timer.record("fastan", 50.0)
        rows = timer.rows()
        assert rows[0][2] == pytest.approx(50.0, abs=1.0)

    def test_unaccounted_time_is_reported(self):
        timer = StepTimer()
        timer._started = time.time() - 100.0
        timer.record("fastan", 60.0)
        rest, percent = timer.unaccounted()
        assert rest == pytest.approx(40.0, abs=1.0)
        assert percent == pytest.approx(40.0, abs=1.0)

    def test_total_is_wall_clock_not_the_sum(self):
        timer = StepTimer()
        timer._started = time.time() - 10.0
        timer.record("a", 1.0)
        assert timer.total == pytest.approx(10.0, abs=0.5)


class TestTable:
    def test_table_shows_time_and_share_per_step(self):
        timer = StepTimer()
        timer._started = time.time() - 1000.0
        timer.record("fastan", 800.0)
        timer.record("classification", 100.0)

        text = "\n".join(format_timing_table(timer, {"fastan": "ok", "classification": "ok"}))

        assert "STEP TIMING" in text
        assert "fastan" in text and "13m20s" in text
        assert "80.0%" in text
        assert "(other)" in text, "untimed work must be visible, not silently dropped"
        assert "TOTAL" in text

    def test_failed_step_is_marked(self):
        timer = StepTimer()
        timer.record("drawing", 3.0)
        text = "\n".join(format_timing_table(timer, {"drawing": "failed"}))
        assert "[failed]" in text

    def test_ok_steps_are_not_annotated(self):
        timer = StepTimer()
        timer.record("drawing", 3.0)
        text = "\n".join(format_timing_table(timer, {"drawing": "ok"}))
        assert "[ok]" not in text

    def test_empty_timer_renders_nothing(self):
        assert format_timing_table(StepTimer(), {}) == []
