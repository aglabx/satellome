"""Reduce a finished satellome output directory to its primary layer.

``satellome compact`` rewrites a run's output so that only the information that
cannot be recomputed from it is kept, and ``satellome expand`` puts the rest
back. The point is not tidiness: the satellite panel produces ~190 MB per
assembly, the roster holds 26,266 of them, and the volume that has to hold the
result is smaller than the result. Without this the panel stops around 77% —
the runners carry disk floors precisely so that they refuse rather than write a
truncated catalogue.

The three ideas the implementation rests on:

* **A policy table, not branches.** :mod:`satellome.compact.policy` maps a file
  *kind* (matched on its trailing suffix — never on a prefix, because the
  basename comes from the input FASTA and routinely disagrees with the
  directory) to a class and an action. Adding an output kind is a row.
* **Never destructive in place.** A dropped file is regenerated and compared to
  itself *before* it is unlinked, a rewritten file is read back and verified
  before the original goes. A crash mid-file costs nothing.
* **A per-directory record.** ``.compact.json`` holds the pre-compaction md5 of
  everything touched and the recipe for everything dropped. Without that record
  compaction is deletion.
"""

from satellome.compact.policy import (  # noqa: F401  (public surface)
    POLICY_VERSION,
    CompactConfig,
    FileClass,
    classify_path,
)

__all__ = ["POLICY_VERSION", "CompactConfig", "FileClass", "classify_path"]
