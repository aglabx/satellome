"""A cheap check that the recipes can actually run here, before promising a saving.

``--dry-run`` used to price the per-copy cut as if the decomposer would always
reproduce what it produced months ago. On a machine carrying a different
arraysplitter build it does not — 16.5 MB was predicted to become 2.8 MB and
became 12.1 MB, because three kinds were correctly kept back. An estimate that
assumes the happy path is not an estimate; it is the best case wearing the word.

So the ledger runs the decomposer over a handful of sub-threshold arrays first
and compares the rows it gets against the rows already on disk. That costs a
fraction of a second and answers the only question that moves the number by
5x: will the per-copy cut happen at all.

This is a probe, not the gate. ``satellome compact --check-recipes`` is the
gate: it checks every kind over the whole file. The probe checks a few arrays,
which is enough to catch a missing binary, a version with a different schema, or
a decomposer that emits different bytes.
"""

import logging
import os
from typing import NamedTuple

from satellome.compact import formats, recipes

logger = logging.getLogger("satellome")

#: Arrays fed to the probe.
#:
#: Sized against a measured failure, not guessed: a locally installed
#: arraysplitter 0.1.0 decomposes 3.4% of this corpus's arrays differently from
#: the build that produced it — a different period, a different monomer count.
#: Five arrays would miss that 84% of the time. Two hundred miss it under 1% of
#: the time, and cost about a tenth of a second.
PROBE_ARRAYS = 200


class ProbeResult(NamedTuple):
    #: Kind names the per-copy cut can be applied to here.
    usable: frozenset
    #: {kind name: why not}
    refused: dict

    def can_filter(self, kind_name):
        return kind_name in self.usable


ALL_USABLE = ProbeResult(frozenset({"monomers", "hors", "decomposed"}), {})


def _corpus_file(run_dir, suffix):
    fastan = os.path.join(run_dir, "fastan")
    if not os.path.isdir(fastan):
        return None, "none"
    for name in sorted(os.listdir(fastan)):
        if name.endswith(suffix):
            return os.path.join(fastan, name), "none"
        if name.endswith(suffix + ".gz"):
            return os.path.join(fastan, name), "gzip"
    return None, "none"


def _rows_for(path, compression, wanted):
    """Rows of a per-copy TSV whose array id is in *wanted*, keyed by id."""
    found = {}
    with formats.open_maybe_gzip(path, compression) as fh:
        fh.readline()
        for raw in fh:
            key = raw.split(b"\t", 1)[0]
            if key in wanted:
                found.setdefault(key, []).append(raw)
            elif found and len(found) == len(wanted):
                break
    return found


def _records_for(path, compression, wanted):
    found = {}
    with formats.open_maybe_gzip(path, compression) as fh:
        for head, body in formats.iter_fasta_records(fh):
            key = head[1:].split(None, 1)[0] if head.startswith(b">") else b""
            if key in wanted:
                found.setdefault(key, []).append((head, body))
    return found


def probe_decomposer(run_dir, master_rows, prefix, min_array_length, workdir):
    """Can the per-copy layer be rebuilt on this machine, for this catalogue?

    Returns a :class:`ProbeResult`; every kind it refuses carries a reason,
    which is what the ledger prints instead of a saving it will not deliver.
    """
    kinds = {
        "monomers": ".monomers.tsv",
        "hors": ".hors.tsv",
        "decomposed": ".decomposed.fasta",
    }
    present = {
        name: _corpus_file(run_dir, suffix)
        for name, suffix in kinds.items()
    }
    present = {name: value for name, value in present.items() if value[0]}
    if not present:
        return ProbeResult(frozenset(), {})

    try:
        binary = recipes.arraysplitter_binary()
    except recipes.RecipeError as e:
        return ProbeResult(frozenset(), {name: str(e) for name in present})

    os.makedirs(workdir, exist_ok=True)
    arrays_fasta = os.path.join(workdir, "probe.fasta")
    chosen = []
    with open(arrays_fasta, "wb") as fh:
        for line in master_rows:
            fields = line.split(b"\t")
            if len(fields) < 15:
                continue
            try:
                length = int(fields[14])
                l_ind, r_ind = int(fields[3]), int(fields[4])
            except ValueError:
                continue
            if length >= min_array_length:
                continue
            array_id = b"%s_%d_%d_%d_%s" % (
                fields[2], l_ind - 1, r_ind, length, fields[5]
            )
            fh.write(b">%s\n%s\n" % (array_id, fields[11]))
            chosen.append(array_id)
            if len(chosen) >= PROBE_ARRAYS:
                break

    if not chosen:
        # Nothing is below the threshold, so the cut removes nothing and there
        # is nothing to be unable to rebuild.
        return ProbeResult(frozenset(present), {})

    try:
        recipes.run_decomposer(
            arrays_fasta, os.path.join(workdir, prefix), threads=1, binary=binary
        )
    except recipes.RecipeError as e:
        return ProbeResult(frozenset(), {name: str(e) for name in present})

    wanted = set(chosen)
    usable = set()
    refused = {}
    for name, suffix in kinds.items():
        if name not in present:
            continue
        produced = os.path.join(workdir, prefix + suffix)
        if not os.path.exists(produced):
            refused[name] = (
                f"this arraysplitter build does not write {suffix}; a different "
                f"version cannot reproduce this catalogue"
            )
            continue
        path, compression = present[name]
        if name == "decomposed":
            theirs = _records_for(path, compression, wanted)
            mine = _records_for(produced, "none", wanted)
        else:
            theirs = _rows_for(path, compression, wanted)
            mine = _rows_for(produced, "none", wanted)
        mismatched = [k for k in wanted if theirs.get(k) != mine.get(k)]
        if mismatched:
            refused[name] = (
                f"the decomposer here produces different rows for "
                f"{len(mismatched)} of {len(wanted)} probed arrays "
                f"(e.g. {mismatched[0].decode('utf-8', 'replace')})"
            )
            continue
        usable.add(name)
    return ProbeResult(frozenset(usable), refused)
