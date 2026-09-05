"""Stop producing, at the end of a run, what compaction would remove anyway.

The panel has ~10,600 assemblies left to compute and the volume holds 1.4 TB.
Retrofitting the 2,901 GB already on disk is one half of that arithmetic; the
other half is that new runs must land in the shape compaction leaves behind,
rather than being written large and shrunk afterwards.

Which kinds go is read off the same policy table the compactor applies, so the
two cannot drift. What stays regardless is named in
:data:`~satellome.compact.policy.PIPELINE_KEEPS_ANYWAY`: the pipeline reads some
of its own derived views, and a run that deleted those would break ``--rerun``.

Nothing here is silent. The removal is one log line naming every kind and the
bytes it freed, and the run manifest is written afterwards, so the manifest
describes the directory that actually exists.
"""

import logging
import os

from satellome.compact.ledger import human, iter_run_files
from satellome.compact.policy import classify_path, pruned_kind_names

logger = logging.getLogger("satellome")


def prune_run_directory(run_dir, log=logger, dry_run=False):
    """Remove the outputs a default run is not asked to keep.

    Returns ``(removed_paths, freed_bytes)``. Files are removed only when their
    kind is droppable by policy; anything unclassified is left alone, because a
    file nobody classified must not disappear because nobody classified it.
    """
    droppable = set(pruned_kind_names())
    removed = []
    freed = 0

    for rel_path in iter_run_files(run_dir):
        match = classify_path(rel_path)
        if match.kind is None or match.kind.name not in droppable:
            continue
        abs_path = os.path.join(run_dir, rel_path)
        try:
            size = os.path.getsize(abs_path)
        except OSError as e:
            log.warning(f"Cannot stat {rel_path} while pruning: {e}")
            continue
        if not dry_run:
            try:
                os.unlink(abs_path)
            except OSError as e:
                # Failing to remove an output is not a reason to fail the run,
                # but it must be visible: the directory then holds more than the
                # summary says it does.
                log.warning(f"Could not remove {rel_path}: {e}")
                continue
        removed.append(rel_path)
        freed += size

    if removed:
        by_kind = {}
        for rel_path in removed:
            by_kind.setdefault(classify_path(rel_path).kind.name, 0)
            by_kind[classify_path(rel_path).kind.name] += 1
        summary = ", ".join(f"{name} x{n}" for name, n in sorted(by_kind.items()))
        verb = "would remove" if dry_run else "removed"
        log.info(
            f"Compact output: {verb} {len(removed)} derived file(s), "
            f"{human(freed)} freed ({summary}). "
            f"Run with --extended-output to keep them, or "
            f"'satellome expand {run_dir}' after a compaction to rebuild them."
        )
    return removed, freed
