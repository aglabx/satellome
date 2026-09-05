"""Is this directory finished, and is it safe to touch?

The obvious guards do not work, and that was checked on the live tree rather
than assumed:

* **``.rc_complete`` is not a usable completeness marker.** 1,020 finished,
  filed catalogues do not carry it — they predate the marker. Gating on it would
  exclude them from compaction permanently.
* **A ``fasta/`` subdirectory is not an in-flight signal either.** 1,015 tree
  directories carry one left over from an earlier campaign. The sign of "still
  computing" and the sign of "old run" coincide.

What actually distinguishes a directory under a live worker, checked mid-run on
``GCA_947311925.1``: no ``.rc_complete``, a ``fasta/`` directory, and **``.sat``
files sitting uncompressed** — compression happens at the end. So the test is
positive and about content, not about markers.

The structural guarantee is stronger than any of these checks: compaction runs
only on the finished tree, never on ``_incoming/``. Staging is where computation
happens; the tree is where finished work lives. A guarantee by construction
cannot be written wrong, and a check can.

Every refusal is a named line in the report. A silent skip is indistinguishable
from success, and at 15,620 directories that difference is the whole point.
"""

import gzip
import os
from typing import List, NamedTuple

#: Path components that mean "this is staging, not the tree".
STAGING_COMPONENTS = ("_incoming", "_superseded")


class Readiness(NamedTuple):
    ok: bool
    reasons: List[str]
    #: The run's basename prefix, once we are confident enough to name it.
    prefix: str = ""


def in_staging(run_dir):
    """True when *run_dir* lies under a staging directory."""
    parts = os.path.abspath(run_dir).split(os.sep)
    return any(component in STAGING_COMPONENTS for component in parts)


def _gzip_readable_to_end(path):
    """Read a gzip member to its end. Returns (ok, detail).

    ``gzip -t`` in Python form. A master table that stops mid-stream is the one
    input that would let compaction produce a smaller broken output and call it
    a success.
    """
    try:
        with gzip.open(path, "rb") as fh:
            for _ in iter(lambda: fh.read(1 << 20), b""):
                pass
    except (OSError, EOFError) as e:
        return False, str(e)
    return True, ""


def find_master(run_dir):
    """Locate the master ``.sat``/``.sat.gz`` and return ``(path, prefix)``.

    Classified by suffix, never by matching the directory name: the basename
    comes from the input FASTA, and directory ``GCA_029289425.3`` legitimately
    holds ``GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic.sat.gz``.
    """
    # Keyed by base name so that ``X.sat`` and ``X.sat.gz`` — which is what an
    # interrupted run leaves behind — count as one master, not two candidates.
    # Otherwise the in-flight signal would be reported as "cannot tell which
    # .sat is the master", which is true but not the useful thing to say.
    candidates = {}
    for name in sorted(os.listdir(run_dir)):
        path = os.path.join(run_dir, name)
        if not os.path.isfile(path):
            continue
        stem = name[:-3] if name.endswith(".gz") else name
        if not stem.endswith(".sat"):
            continue
        base = stem[: -len(".sat")]
        # A view, not the master: anything with a second extension in front.
        if base.rsplit(".", 1)[-1] in {
            "1kb", "10kb", "100kb", "1000kb", "micro", "pmicro", "tssr",
            "complex", "original",
        }:
            continue
        # Prefer the compressed copy when both are present.
        if base not in candidates or name.endswith(".gz"):
            candidates[base] = path
    if len(candidates) == 1:
        base, path = next(iter(candidates.items()))
        return path, base
    return None, ""


def check(run_dir, require_bed=True):
    """Positive completeness test for one tree directory.

    Returns a :class:`Readiness`. Anything false comes with the reason, by name.
    """
    reasons = []

    if not os.path.isdir(run_dir):
        return Readiness(False, [f"not a directory: {run_dir}"])

    if in_staging(run_dir):
        return Readiness(
            False,
            [
                f"lies under {'/'.join(STAGING_COMPONENTS)} — compaction runs on "
                f"the finished tree only, never on staging"
            ],
        )

    master, prefix = find_master(run_dir)
    if master is None:
        found = [n for n in sorted(os.listdir(run_dir)) if n.endswith((".sat", ".sat.gz"))]
        if not found:
            reasons.append("no master .sat/.sat.gz present")
        else:
            reasons.append(
                "cannot tell which .sat is the master (candidates: "
                + ", ".join(found) + ")"
            )
        return Readiness(False, reasons)

    if os.path.getsize(master) == 0:
        reasons.append(f"{os.path.basename(master)} is empty")
    elif master.endswith(".gz"):
        ok, detail = _gzip_readable_to_end(master)
        if not ok:
            reasons.append(
                f"{os.path.basename(master)} does not decompress to the last "
                f"byte ({detail}) — compacting it would produce a smaller "
                f"broken file"
            )

    # An uncompressed .sat sitting beside a .sat.gz of the same kind is the
    # signature of a run still in progress: compression happens at the end.
    for name in sorted(os.listdir(run_dir)):
        if name.endswith(".sat") and os.path.exists(os.path.join(run_dir, name + ".gz")):
            reasons.append(
                f"{name} sits uncompressed beside {name}.gz — a run is probably "
                f"still writing here"
            )

    fastan_dir = os.path.join(run_dir, "fastan")
    if require_bed:
        beds = (
            [n for n in os.listdir(fastan_dir) if n.endswith(".bed")]
            if os.path.isdir(fastan_dir)
            else []
        )
        if not beds:
            reasons.append("no fastan/*.bed present")

    if not os.path.exists(os.path.join(run_dir, "results.yaml")):
        reasons.append("no results.yaml present")

    partials = [
        os.path.join(root, name)
        for root, _dirs, files in os.walk(run_dir)
        for name in files
        if name.endswith(".partial")
    ]
    if partials:
        reasons.append(
            f"{len(partials)} in-progress .partial file(s) present, e.g. "
            f"{os.path.relpath(partials[0], run_dir)}"
        )

    return Readiness(not reasons, reasons, prefix)
