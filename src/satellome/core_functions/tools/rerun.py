"""Re-run individual steps against a finished output directory.

A full run is expensive: the tandem-repeat search dominates it and its result
does not change when you decide you also want a UCSC track, a redrawn report,
or annotations from a GFF you only just received. This module lets those later
steps run again on what is already there.

How the original options are recovered: the run manifest records the run's own
argv, project, input FASTA, taxon and genome size. Re-parsing that argv gives
back the cutoffs, `--gff` and the rest, so a re-run does not depend on the user
remembering the exact command — and cannot silently apply *different*
parameters than the run it is amending, which would produce an output directory
whose files disagree with each other.

The manifest stays the record of what the directory contains: a re-run updates
the status of the steps it ran, refreshes the file inventory, and appends to a
``reruns`` list, so `--verify-run` keeps working and the history is visible.
"""

import datetime
import logging
import os
from typing import Callable, Dict, List, NamedTuple, Optional, Sequence, Tuple

from satellome.core_functions.tools.run_manifest import (
    MANIFEST_NAME,
    ManifestError,
    load_manifest,
)

logger = logging.getLogger(__name__)


class RerunError(Exception):
    """The re-run cannot proceed, with a reason the user can act on."""


class StepSpec(NamedTuple):
    name: str
    description: str
    # Files that must already exist, as {label: path-builder(settings)}.
    requires: Callable[[dict], Dict[str, str]]


def _sat_file(settings) -> Dict[str, str]:
    return {"tandem repeat file": settings["trf_file"]}


def _sat_and_fasta(settings) -> Dict[str, str]:
    return {
        "tandem repeat file": settings["trf_file"],
        "genome FASTA": settings["fasta_file"],
    }


# Steps that can be re-run cheaply because they consume the search results
# rather than producing them. The expensive search steps are deliberately not
# here: re-running those is what --force is for, and quietly offering them
# under a flag named "rerun" would hide hours of work behind a one-word option.
RERUNNABLE: Tuple[StepSpec, ...] = (
    StepSpec("classification", "Classify repeats and write the classified .sat", _sat_file),
    StepSpec("annotations", "Intersect repeats with a GFF/RepeatMasker annotation", _sat_file),
    StepSpec("sat_family", "Cluster satellite DNA families", _sat_file),
    StepSpec("drawing", "Redraw the plots and the HTML report", _sat_file),
    StepSpec("ucsc_track", "Write the UCSC browser track (BED9, chrom.sizes, bigBed)", _sat_and_fasta),
)

RERUNNABLE_BY_NAME: Dict[str, StepSpec] = {s.name: s for s in RERUNNABLE}


def format_step_list() -> List[str]:
    lines = ["Steps accepted by --rerun:"]
    width = max(len(s.name) for s in RERUNNABLE)
    for spec in RERUNNABLE:
        lines.append(f"  {spec.name:<{width}}  {spec.description}")
    lines.append("")
    lines.append("Example: satellome --rerun ucsc_track -o /path/to/output")
    lines.append(
        "The tandem-repeat search is not re-runnable here on purpose; use "
        "--force on a normal run to recompute it."
    )
    return lines


def parse_step_names(value: str) -> List[str]:
    """Split and validate a comma-separated step list."""
    names = [part.strip() for part in str(value).split(",") if part.strip()]
    if not names:
        raise RerunError("--rerun needs at least one step name (see --list-steps)")

    unknown = [n for n in names if n not in RERUNNABLE_BY_NAME]
    if unknown:
        raise RerunError(
            f"unknown step(s): {', '.join(unknown)}. "
            f"Known steps: {', '.join(RERUNNABLE_BY_NAME)}"
        )

    # Preserve the pipeline's own order, so re-running two steps applies them
    # in the same sequence a full run would.
    ordered = [s.name for s in RERUNNABLE if s.name in names]
    return ordered


def load_run_context(output_dir: str) -> dict:
    """Recover a finished run's parameters from its manifest.

    Raises RerunError with an actionable message when the directory is not a
    finished run - the alternative, guessing a project name and an input path,
    would write files that silently belong to a different analysis.
    """
    if not os.path.isdir(output_dir):
        raise RerunError(f"not a directory: {output_dir}")

    try:
        manifest = load_manifest(output_dir)
    except ManifestError as error:
        raise RerunError(
            f"{output_dir} has no usable {MANIFEST_NAME} ({error}). "
            "A re-run amends a finished run; run the full pipeline first."
        ) from error

    argv = manifest.get("argv")
    if not isinstance(argv, list):
        raise RerunError(
            f"{MANIFEST_NAME} does not record the original command line, so the "
            "run's options cannot be recovered. Re-run the full pipeline once "
            "with this version to write a manifest that does."
        )

    project = manifest.get("project")
    if not project:
        raise RerunError(f"{MANIFEST_NAME} does not name a project")

    return {
        "manifest": manifest,
        "argv": [str(a) for a in argv],
        "project": project,
        "input_fasta": manifest.get("input_fasta"),
        "taxon_name": manifest.get("taxon_name"),
        "taxid": manifest.get("taxid"),
        "genome_size": manifest.get("genome_size") or 0,
        "steps": dict(manifest.get("steps") or {}),
    }


def check_requirements(step_names: Sequence[str], settings: dict) -> None:
    """Fail before doing any work if a step's inputs are not there."""
    missing: List[str] = []
    for name in step_names:
        for label, path in RERUNNABLE_BY_NAME[name].requires(settings).items():
            if not path or not os.path.exists(path):
                missing.append(f"{name} needs the {label}: {path or '(unknown path)'}")
    if missing:
        raise RerunError(
            "the existing run does not have what these steps need:\n  "
            + "\n  ".join(missing)
        )


def merge_steps(previous: Dict[str, str], just_run: Dict[str, str]) -> Dict[str, str]:
    """Statuses after a re-run: previous ones kept, re-run ones replaced."""
    merged = dict(previous)
    merged.update(just_run)
    return merged


def record_rerun(manifest: dict, step_names: Sequence[str], version: str) -> List[dict]:
    """Append this re-run to the manifest's history and return the list."""
    history = list(manifest.get("reruns") or [])
    history.append({
        "at": datetime.datetime.now().isoformat(timespec="seconds"),
        "steps": list(step_names),
        "satellome_version": version,
    })
    return history
