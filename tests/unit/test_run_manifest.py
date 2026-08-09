"""Tests for the run manifest — the regression guard for silent truncation.

The production incident these tests encode: a driver script compressed output
directories while files were still being written. gzip warned `file size changed
while zipping` on stderr and produced a *valid* .gz of a *truncated* file. Every
check the consumer had then passed —

    * `<project>.sat` exists                        → yes
    * `fastan/*.monomers.tsv` exists                → yes
    * `gzip -t` on the archives                     → yes

— and the truncated directory entered a dashboard as data, unmarked. Each test
below therefore asserts *both* that the old signals still look fine and that
verify_run reports the specific problem, so the tests fail if verification ever
regresses to "the files are there, ship it".
"""

import gzip
import json
import os

import pytest

from satellome.core_functions.tools.run_manifest import (
    MANIFEST_NAME,
    ManifestError,
    build_manifest,
    gz_uncompressed_size,
    load_manifest,
    verify_run,
    write_run_manifest,
)


STEPS_OK = {"fastan": "ok", "classification": "ok", "drawing": "ok"}


def _make_run(tmp_path, sat_bytes=b"project\ttrf_id\n" + b"row\n" * 100):
    """Create a plausible finished run directory (unmanifested)."""
    (tmp_path / "fastan").mkdir()
    sat = tmp_path / "genome.sat"
    sat.write_bytes(sat_bytes)
    monomers = tmp_path / "fastan" / "genome.monomers.tsv"
    monomers.write_text("monomer\tcount\nACGT\t12\n")
    return sat, monomers


def _finish(tmp_path, steps=None, records=None):
    manifest = build_manifest(
        str(tmp_path), project="genome", version="test",
        steps=steps or STEPS_OK, records=records,
    )
    write_run_manifest(str(tmp_path), manifest)
    return manifest


def test_complete_run_verifies(tmp_path):
    _make_run(tmp_path)
    _finish(tmp_path)

    result = verify_run(str(tmp_path))

    assert result.ok, result.problems
    assert result.problems == []
    # .sat, monomers.tsv — the manifest itself is not one of its own entries.
    assert result.checked == 2
    assert MANIFEST_NAME not in {e["path"] for e in result.manifest["files"]}


def test_truncated_sat_is_caught_although_the_file_exists(tmp_path):
    """The exact silence: file present, plausible, and 40% shorter."""
    sat, monomers = _make_run(tmp_path)
    _finish(tmp_path)
    full_size = sat.stat().st_size

    with open(sat, "r+b") as fh:  # simulate a reader-side truncated copy
        fh.truncate(full_size // 2)

    # Everything the old check looked at still passes:
    assert sat.exists() and monomers.exists()

    result = verify_run(str(tmp_path))

    assert not result.ok
    assert len(result.problems) == 1
    problem = result.problems[0]
    assert "genome.sat" in problem
    assert "truncated" in problem
    assert str(full_size) in problem and str(full_size // 2) in problem


def test_grown_file_is_caught(tmp_path):
    """A file that kept growing after the manifest means the run was raced."""
    sat, _ = _make_run(tmp_path)
    _finish(tmp_path)
    with open(sat, "ab") as fh:
        fh.write(b"row\n" * 10)

    result = verify_run(str(tmp_path))

    assert not result.ok
    assert any("grown" in p and "genome.sat" in p for p in result.problems)


def test_missing_file_is_caught(tmp_path):
    sat, _ = _make_run(tmp_path)
    _finish(tmp_path)
    os.unlink(sat)

    result = verify_run(str(tmp_path))

    assert not result.ok
    assert any(
        "genome.sat" in p and "missing" in p for p in result.problems
    ), result.problems


def test_gzip_of_truncated_file_is_caught_though_gzip_itself_is_valid(tmp_path):
    """gzip -t passes on this archive; the ISIZE trailer still gives it away."""
    sat, _ = _make_run(tmp_path)
    _finish(tmp_path)
    full_size = sat.stat().st_size
    partial = sat.read_bytes()[: full_size // 3]

    with gzip.open(str(sat) + ".gz", "wb") as fh:
        fh.write(partial)
    os.unlink(sat)

    # The archive is a perfectly valid gzip member — decompressing works.
    with gzip.open(str(sat) + ".gz", "rb") as fh:
        assert fh.read() == partial

    result = verify_run(str(tmp_path))

    assert not result.ok
    problem = next(p for p in result.problems if "genome.sat.gz" in p)
    assert str(len(partial)) in problem
    assert str(full_size) in problem
    assert "incomplete" in problem


def test_gzip_of_complete_file_verifies(tmp_path):
    """The legitimate workflow — verify, then compress — must stay green."""
    sat, monomers = _make_run(tmp_path)
    _finish(tmp_path)
    data = sat.read_bytes()
    with gzip.open(str(sat) + ".gz", "wb") as fh:
        fh.write(data)
    os.unlink(sat)

    result = verify_run(str(tmp_path))

    assert result.ok, result.problems
    assert result.warnings == []


def test_failed_step_is_visible_in_verification(tmp_path):
    """Drawing crashed, data is complete: the run is still flagged, not silent."""
    _make_run(tmp_path)
    _finish(tmp_path, steps={"fastan": "ok", "classification": "ok", "drawing": "failed"})

    result = verify_run(str(tmp_path))

    assert not result.ok
    assert result.problems == ["step 'drawing' failed during the run"]
    assert result.manifest["complete"] is False


def test_skipped_steps_do_not_flag_the_run(tmp_path):
    _make_run(tmp_path)
    _finish(tmp_path, steps={"trf_search": "skipped", "fastan": "ok", "drawing": "ok"})

    assert verify_run(str(tmp_path)).ok


def test_no_manifest_means_not_verifiable(tmp_path):
    _make_run(tmp_path)  # files exist, run never finished

    result = verify_run(str(tmp_path))

    assert not result.ok
    assert MANIFEST_NAME in result.problems[0]
    assert "did not finish" in result.problems[0]
    assert result.manifest is None


def test_corrupt_manifest_is_distinguished_from_missing(tmp_path):
    _make_run(tmp_path)
    _finish(tmp_path)
    (tmp_path / MANIFEST_NAME).write_text("{not json at all")

    with pytest.raises(ManifestError, match="Corrupt"):
        load_manifest(str(tmp_path))

    result = verify_run(str(tmp_path))
    assert not result.ok
    assert "Corrupt" in result.problems[0]


def test_manifest_with_wrong_shape_is_corrupt(tmp_path):
    _make_run(tmp_path)
    (tmp_path / MANIFEST_NAME).write_text(json.dumps(["not", "an", "object"]))

    with pytest.raises(ManifestError, match="expected an object with a 'files' list"):
        load_manifest(str(tmp_path))


def test_leftover_partial_is_recorded_and_flagged(tmp_path):
    """A stray .partial means some write never completed — never a silent skip."""
    _make_run(tmp_path)
    (tmp_path / "genome.10kb.sat.partial").write_text("half filtered\n")
    manifest = _finish(tmp_path)

    assert manifest["leftover_partials"] == ["genome.10kb.sat.partial"]
    # The partial is not offered as an output file of the run.
    assert "genome.10kb.sat.partial" not in {e["path"] for e in manifest["files"]}

    result = verify_run(str(tmp_path))
    assert not result.ok
    assert any("in-progress file left behind" in p for p in result.problems)


def test_unlisted_extra_file_is_a_warning_not_a_failure(tmp_path):
    _make_run(tmp_path)
    _finish(tmp_path)
    (tmp_path / "notes.txt").write_text("added by hand later\n")

    result = verify_run(str(tmp_path))

    assert result.ok
    assert any("notes.txt" in w and "not listed" in w for w in result.warnings)


def test_record_counts_are_carried_into_the_manifest(tmp_path):
    _make_run(tmp_path)
    manifest = _finish(tmp_path, records={"genome.sat": 1523})

    entry = next(e for e in manifest["files"] if e["path"] == "genome.sat")
    assert entry["records"] == 1523


def test_nested_files_are_recorded_with_relative_paths(tmp_path):
    _make_run(tmp_path)
    manifest = _finish(tmp_path)

    paths = {e["path"] for e in manifest["files"]}
    assert paths == {"genome.sat", os.path.join("fastan", "genome.monomers.tsv")}


def test_gz_uncompressed_size_matches_payload(tmp_path):
    payload = b"x" * 5000
    gz = tmp_path / "p.gz"
    with gzip.open(gz, "wb") as fh:
        fh.write(payload)

    assert gz_uncompressed_size(str(gz)) == 5000


def test_gz_uncompressed_size_none_for_non_gzip(tmp_path):
    tiny = tmp_path / "t.gz"
    tiny.write_bytes(b"short")
    assert gz_uncompressed_size(str(tiny)) is None
