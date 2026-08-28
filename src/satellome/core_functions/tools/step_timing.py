"""Wall-clock timing for each pipeline step, recorded in the run manifest.

A run takes hours and, until now, produced no record of where those hours went.
That makes every performance question - "is the search the cost, or is it the
Python that consumes its output?" - unanswerable after the fact, and answerable
only by re-running with a stopwatch.

The timings are written into ``run_manifest.json`` as ``steps_timing`` and
printed as a table at the end of the run. A step that fails is still timed: how
long a failure took is often the interesting part.
"""

import logging
import time
from contextlib import contextmanager
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


class StepTimer:
    """Collects per-step wall time in the order the steps ran."""

    def __init__(self):
        self._order: List[str] = []
        self._seconds: Dict[str, float] = {}
        self._started = time.time()

    @contextmanager
    def step(self, name: str):
        """Time a step. Records the elapsed time even if the body raises."""
        start = time.time()
        try:
            yield
        finally:
            self.record(name, time.time() - start)

    def record(self, name: str, seconds: float) -> None:
        if name not in self._seconds:
            self._order.append(name)
            self._seconds[name] = 0.0
        # A step run more than once in a run accumulates rather than overwrites.
        self._seconds[name] += float(seconds)

    @property
    def total(self) -> float:
        """Wall time since the timer was created, not the sum of the steps.

        They differ, and the difference is real work: validation, genome-size
        computation, manifest writing. Reporting the sum as the total would
        hide it.
        """
        return time.time() - self._started

    def as_dict(self) -> Dict[str, float]:
        return {name: round(self._seconds[name], 2) for name in self._order}

    def rows(self) -> List[Tuple[str, float, float]]:
        """(step, seconds, percent-of-wall-clock), slowest first."""
        total = max(self.total, 1e-9)
        rows = [(name, self._seconds[name], 100.0 * self._seconds[name] / total)
                for name in self._order]
        return sorted(rows, key=lambda r: r[1], reverse=True)

    def unaccounted(self) -> Tuple[float, float]:
        """Wall time not inside any timed step, and its percentage."""
        total = max(self.total, 1e-9)
        rest = total - sum(self._seconds.values())
        return rest, 100.0 * rest / total


def format_duration(seconds: float) -> str:
    seconds = float(seconds)
    if seconds < 60:
        return f"{seconds:.1f}s"
    if seconds < 3600:
        return f"{int(seconds // 60)}m{int(seconds % 60):02d}s"
    return f"{int(seconds // 3600)}h{int((seconds % 3600) // 60):02d}m"


def format_timing_table(timer: StepTimer, steps: Optional[Dict[str, str]] = None) -> List[str]:
    """Render the timing summary as log lines."""
    steps = steps or {}
    rows = timer.rows()
    if not rows:
        return []

    width = max(len(name) for name, _, _ in rows)
    lines = ["STEP TIMING"]
    for name, seconds, percent in rows:
        status = steps.get(name, "")
        suffix = f"  [{status}]" if status and status != "ok" else ""
        lines.append(f"  {name:<{width}}  {format_duration(seconds):>8}  {percent:5.1f}%{suffix}")

    rest, rest_percent = timer.unaccounted()
    lines.append(
        f"  {'(other)':<{width}}  {format_duration(rest):>8}  {rest_percent:5.1f}%"
        "   validation, genome size, manifest"
    )
    lines.append(f"  {'TOTAL':<{width}}  {format_duration(timer.total):>8}  100.0%")
    return lines
