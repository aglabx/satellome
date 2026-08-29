"""A raising step must not cost the run its manifest.

On CHM13 a binary built for another platform made a step raise
`OSError: Exec format error`. The exception took the process down before the
manifest was written, so 65 minutes of completed work left an output directory
indistinguishable from an unfinished run — the exact thing run_manifest.json
exists to rule out.
"""

import pytest

from satellome.main import guarded_step
from satellome.core_functions.tools.step_timing import StepTimer


class TestGuardedStep:
    def test_success_is_recorded_and_returned(self):
        steps, timer = {}, StepTimer()
        assert guarded_step("drawing", timer, steps, lambda: True) is True
        assert steps["drawing"] == "ok"
        assert "drawing" in timer.as_dict()

    def test_a_falsy_return_is_a_failure_not_an_error(self):
        steps, timer = {}, StepTimer()
        assert guarded_step("drawing", timer, steps, lambda: False) is False
        assert steps["drawing"] == "failed"

    def test_a_step_returning_none_counts_as_ok(self):
        """add_annotations returns nothing; that is success, not failure."""
        steps, timer = {}, StepTimer()
        guarded_step("annotations", timer, steps, lambda: None)
        assert steps["annotations"] == "ok"

    def test_an_exception_becomes_that_step_s_failure(self, caplog):
        steps, timer = {}, StepTimer()

        def explode():
            raise OSError(8, "Exec format error", "/opt/bin/sat-family")

        with caplog.at_level("ERROR"):
            result = guarded_step("sat_family", timer, steps, explode)

        assert result is False
        assert steps["sat_family"] == "failed", "the manifest must record what happened"
        assert "Exec format error" in caplog.text
        assert "sat_family" in caplog.text

    def test_the_exception_does_not_propagate(self):
        """If it escaped, the manifest write further down would never run."""
        steps, timer = {}, StepTimer()

        def explode():
            raise RuntimeError("boom")

        guarded_step("drawing", timer, steps, explode)  # must not raise
        assert steps["drawing"] == "failed"

    def test_a_failing_step_is_still_timed(self):
        steps, timer = {}, StepTimer()

        def explode():
            raise ValueError("nope")

        guarded_step("classification", timer, steps, explode)
        assert "classification" in timer.as_dict()

    def test_later_steps_still_run_after_one_fails(self):
        steps, timer = {}, StepTimer()

        def explode():
            raise RuntimeError("boom")

        guarded_step("sat_family", timer, steps, explode)
        guarded_step("drawing", timer, steps, lambda: True)

        assert steps == {"sat_family": "failed", "drawing": "ok"}

    def test_arguments_are_forwarded(self):
        steps, timer = {}, StepTimer()
        seen = {}

        def step(settings, force, min_length=0):
            seen.update(settings=settings, force=force, min_length=min_length)
            return True

        guarded_step("ucsc_track", timer, steps, step, {"a": 1}, True, min_length=5)
        assert seen == {"settings": {"a": 1}, "force": True, "min_length": 5}

    def test_a_dict_return_is_passed_through(self):
        """build_ucsc_track returns the files it produced."""
        steps, timer = {}, StepTimer()
        produced = guarded_step("ucsc_track", timer, steps, lambda: {"bed": "/x.bed"})
        assert produced == {"bed": "/x.bed"}
        assert steps["ucsc_track"] == "ok"

    def test_an_empty_dict_return_is_a_failure(self):
        steps, timer = {}, StepTimer()
        assert guarded_step("ucsc_track", timer, steps, lambda: {}) == {}
        assert steps["ucsc_track"] == "failed"
