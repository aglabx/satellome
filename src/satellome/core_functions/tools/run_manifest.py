"""Run manifest: a verifiable statement of what a satellome run produced.

Why this exists
---------------
Before this, the only thing a downstream consumer could ask about an output
directory was "does ``<project>.sat`` exist, and does ``fastan/*.monomers.tsv``
exist?". That check cannot distinguish a complete run from a directory whose
files were read while still being written — a truncated file exists just as hard
as a complete one, and ``gzip -t`` on a ``.gz`` made from a truncated read passes
too, because the *container* is intact. The result is a partial dataset entering
a dashboard as data, unmarked. That is exactly the silent corruption the project
rules forbid.

The manifest is written **last**, after every step, and records the byte size of
every file the run produced plus the status of each pipeline step. It turns two
unanswerable questions into answerable ones:

* "Is this directory a finished run?"  → is there a manifest?
* "Are these files the ones the run wrote?" → do the sizes still match?

A run whose drawing step crashed still gets a manifest, with that step marked
``failed``: the data files are complete and usable, and the failure is visible
in the artifact itself rather than only in a return code nobody kept.

Compressed copies
-----------------
Consumers routinely gzip the directory afterwards. If a recorded file is gone
but ``<file>.gz`` is there, verification reads the gzip ISIZE trailer — the
uncompressed length the compressor actually consumed — and compares it to the
recorded size. A gzip of a truncated (or still-growing) file disagrees there,
which is precisely the case that used to pass unnoticed. ISIZE is stored modulo
2**32, so for files ≥ 4 GiB the comparison is modular; it can therefore miss a
truncation that happens to be an exact multiple of 4 GiB, and the report says so
rather than claiming a stronger guarantee than it has.
"""

import json
import logging
import os
import platform
import socket
import struct
import sys
from collections import namedtuple

from satellome.core_functions.tools.atomic_io import PARTIAL_SUFFIX, atomic_output

logger = logging.getLogger("satellome")

MANIFEST_NAME = "run_manifest.json"
MANIFEST_VERSION = 1
GZ_ISIZE_MODULUS = 2 ** 32

#: Files that describe the run rather than being its output.
_EXCLUDED_NAMES = {MANIFEST_NAME, ".satellome.lock"}

VerifyResult = namedtuple(
    "VerifyResult", ["ok", "problems", "warnings", "manifest", "checked"]
)


class ManifestError(Exception):
    """The manifest exists but cannot be used (missing, corrupt, wrong shape)."""


def _iter_output_files(output_dir):
    """Yield (relative_path, size) for every regular file under ``output_dir``.

    In-progress ``*.partial`` files are skipped: by construction they are not
    outputs, and their presence at manifest time is reported separately.
    """
    for root, _dirs, files in os.walk(output_dir):
        for name in sorted(files):
            if name in _EXCLUDED_NAMES:
                continue
            abs_path = os.path.join(root, name)
            rel_path = os.path.relpath(abs_path, output_dir)
            if name.endswith(PARTIAL_SUFFIX):
                continue
            if os.path.islink(abs_path) or not os.path.isfile(abs_path):
                continue
            try:
                size = os.path.getsize(abs_path)
            except OSError as e:
                # An unreadable file inside our own output directory is not an
                # expected absence — surface it instead of dropping the entry.
                logger.error(f"Cannot stat output file {abs_path}: {e}")
                raise
            yield rel_path, size


def _leftover_partials(output_dir):
    """Relative paths of ``*.partial`` files still present under ``output_dir``."""
    leftovers = []
    for root, _dirs, files in os.walk(output_dir):
        for name in files:
            if name.endswith(PARTIAL_SUFFIX):
                leftovers.append(
                    os.path.relpath(os.path.join(root, name), output_dir)
                )
    return sorted(leftovers)


def build_manifest(output_dir, project, version, steps, extra=None, records=None):
    """Build the manifest dict for a finished run.

    Args:
        output_dir: run directory to inventory.
        project: project name of the run.
        version: satellome version string.
        steps: ``{step_name: "ok" | "failed" | "skipped"}``.
        extra: optional dict merged into the manifest (input file, taxon, ...).
        records: optional ``{relative_path: record_count}`` for files whose
            record count the pipeline already knows (e.g. extracted arrays).
    """
    records = records or {}
    files = []
    for rel_path, size in _iter_output_files(output_dir):
        entry = {"path": rel_path, "bytes": size}
        if rel_path in records and records[rel_path] is not None:
            entry["records"] = int(records[rel_path])
        files.append(entry)

    manifest = {
        "manifest_version": MANIFEST_VERSION,
        "satellome_version": version,
        "project": project,
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "pid": os.getpid(),
        "argv": sys.argv[1:],
        "steps": dict(steps),
        "complete": all(s != "failed" for s in steps.values()),
        "files": files,
        "leftover_partials": _leftover_partials(output_dir),
    }
    if extra:
        manifest.update(extra)
    return manifest


def write_run_manifest(output_dir, manifest):
    """Write the manifest atomically as the very last act of a run."""
    path = os.path.join(output_dir, MANIFEST_NAME)
    with atomic_output(path) as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
        fh.write("\n")
    logger.info(f"✓ Run manifest written: {path} ({len(manifest['files'])} files)")
    if manifest.get("leftover_partials"):
        logger.warning(
            "Leftover in-progress files found at manifest time: "
            + ", ".join(manifest["leftover_partials"])
        )
    return path


def load_manifest(output_dir):
    """Load and shape-check the manifest of ``output_dir``.

    Raises:
        ManifestError: no manifest, unreadable, invalid JSON, or wrong shape.
            A corrupt manifest is never treated as "no manifest" — the caller
            must be able to tell "this run never finished" from "this run's
            record of itself is damaged".
    """
    path = os.path.join(output_dir, MANIFEST_NAME)
    if not os.path.exists(path):
        raise ManifestError(
            f"No {MANIFEST_NAME} in {output_dir}: the run did not finish, or it "
            f"predates run manifests"
        )
    try:
        with open(path, "r") as fh:
            manifest = json.load(fh)
    except json.JSONDecodeError as e:
        raise ManifestError(f"Corrupt {MANIFEST_NAME} in {output_dir}: {e}") from e
    except OSError as e:
        raise ManifestError(f"Unreadable {MANIFEST_NAME} in {output_dir}: {e}") from e

    if not isinstance(manifest, dict) or not isinstance(manifest.get("files"), list):
        raise ManifestError(
            f"Corrupt {MANIFEST_NAME} in {output_dir}: expected an object with a "
            f"'files' list, got {type(manifest).__name__}"
        )
    return manifest


def gz_uncompressed_size(gz_path):
    """Uncompressed length recorded in a gzip file's ISIZE trailer.

    Returns None if the file is too short to hold a trailer or cannot be read
    (i.e. it is not a usable gzip member at all).
    """
    try:
        size = os.path.getsize(gz_path)
        if size < 18:  # 10-byte header + 8-byte trailer minimum
            return None
        with open(gz_path, "rb") as fh:
            fh.seek(-4, os.SEEK_END)
            return struct.unpack("<I", fh.read(4))[0]
    except OSError as e:
        logger.warning(f"Could not read gzip trailer of {gz_path}: {e}")
        return None


def verify_run(output_dir):
    """Check an output directory against its manifest.

    Returns a :class:`VerifyResult`. ``ok`` is False if the run is not verifiably
    complete: no/corrupt manifest, a failed step, a missing recorded file, or a
    size that no longer matches what the run wrote (the truncated-while-copying
    case). ``warnings`` holds things a consumer should know but that do not by
    themselves condemn the directory (modular ISIZE comparison on huge files,
    files present that the manifest does not list).
    """
    problems = []
    warnings = []

    try:
        manifest = load_manifest(output_dir)
    except ManifestError as e:
        return VerifyResult(False, [str(e)], warnings, None, 0)

    failed_steps = [
        name for name, status in (manifest.get("steps") or {}).items()
        if status == "failed"
    ]
    for name in sorted(failed_steps):
        problems.append(f"step '{name}' failed during the run")

    for leftover in manifest.get("leftover_partials") or []:
        problems.append(
            f"in-progress file left behind: {leftover} (a write did not complete)"
        )

    checked = 0
    for entry in manifest["files"]:
        rel_path = entry.get("path")
        expected = entry.get("bytes")
        if rel_path is None or expected is None:
            problems.append(f"malformed manifest entry: {entry!r}")
            continue

        abs_path = os.path.join(output_dir, rel_path)
        gz_path = f"{abs_path}.gz"

        if os.path.isfile(abs_path):
            found = os.path.getsize(abs_path)
            checked += 1
            if found != expected:
                verb = "truncated" if found < expected else "grown"
                problems.append(
                    f"{rel_path}: {verb} since the run wrote it "
                    f"(manifest {expected} bytes, on disk {found})"
                )
        elif os.path.isfile(gz_path):
            isize = gz_uncompressed_size(gz_path)
            checked += 1
            if isize is None:
                problems.append(
                    f"{rel_path}.gz: not a readable gzip member "
                    f"(cannot confirm it holds the {expected} bytes of {rel_path})"
                )
                continue
            if isize != expected % GZ_ISIZE_MODULUS:
                problems.append(
                    f"{rel_path}.gz: compressed copy holds {isize} uncompressed "
                    f"bytes, the run wrote {expected} — the file was read while "
                    f"incomplete, or changed before compression"
                )
            elif expected >= GZ_ISIZE_MODULUS:
                warnings.append(
                    f"{rel_path}.gz: size confirmed modulo 4 GiB only "
                    f"(gzip ISIZE is 32-bit and the file is {expected} bytes)"
                )
        else:
            problems.append(
                f"{rel_path}: recorded by the run but missing (no plain or .gz copy)"
            )

    recorded = {entry.get("path") for entry in manifest["files"]}
    try:
        on_disk = {rel for rel, _ in _iter_output_files(output_dir)}
    except OSError as e:
        problems.append(f"cannot inventory {output_dir}: {e}")
        on_disk = set()
    for rel in sorted(on_disk - recorded):
        if rel.endswith(".gz") and rel[: -len(".gz")] in recorded:
            continue
        warnings.append(f"{rel}: present but not listed in the manifest")

    return VerifyResult(not problems, problems, warnings, manifest, checked)


def format_verify_result(output_dir, result):
    """Render a :class:`VerifyResult` as human-readable lines."""
    lines = []
    if result.ok:
        lines.append(f"OK: {output_dir} matches its manifest ({result.checked} files verified)")
    else:
        lines.append(f"FAILED: {output_dir} is not a verifiably complete run")
        for problem in result.problems:
            lines.append(f"  ! {problem}")
    for warning in result.warnings:
        lines.append(f"  - {warning}")
    return lines
