"""Tests for atomic publication of output files.

The behaviour under test is what a *concurrent reader* can observe. Outputs used
to be streamed straight into their final path, so a driver script that listed the
output directory could open a file that was still growing, gzip it, and get a
structurally valid .gz of a truncated file — invisible to `gzip -t` and to any
existence check. So these tests assert on the state of the final path *during*
the write, not merely that the write eventually produced the right bytes.
"""

import os

import pytest

from satellome.core_functions.tools.atomic_io import (
    PARTIAL_SUFFIX,
    atomic_output,
    atomic_outputs,
    discard,
    partial_path,
    publish,
)


def test_final_path_absent_until_write_completes(tmp_path):
    target = tmp_path / "out.sat"

    with atomic_output(str(target)) as fh:
        fh.write("partial data")
        fh.flush()
        # A reader listing the directory now must not see the final name at all.
        assert not target.exists()
        assert (tmp_path / f"out.sat{PARTIAL_SUFFIX}").exists()

    assert target.read_text() == "partial data"
    assert not (tmp_path / f"out.sat{PARTIAL_SUFFIX}").exists()


def test_failed_write_publishes_nothing_and_leaves_no_partial(tmp_path):
    target = tmp_path / "out.sat"

    with pytest.raises(ValueError, match="extraction blew up"):
        with atomic_output(str(target)) as fh:
            fh.write("half a record")
            raise ValueError("extraction blew up")

    assert not target.exists(), "a failed write must not publish a file"
    assert not (tmp_path / f"out.sat{PARTIAL_SUFFIX}").exists()


def test_failed_rewrite_keeps_previous_complete_version(tmp_path):
    """A rerun that crashes must not destroy the artifact already on disk."""
    target = tmp_path / "out.sat"
    target.write_text("complete output of the previous run\n")

    with pytest.raises(RuntimeError):
        with atomic_output(str(target)) as fh:
            fh.write("truncated rewrite")
            raise RuntimeError("boom")

    assert target.read_text() == "complete output of the previous run\n"


def test_sibling_outputs_are_published_together(tmp_path):
    sat = tmp_path / "out.sat"
    fasta = tmp_path / "out.fasta"

    with atomic_outputs(str(sat), str(fasta)) as (sat_fh, fasta_fh):
        sat_fh.write("record\n")
        fasta_fh.write(">seq\nACGT\n")
        assert not sat.exists()
        assert not fasta.exists()

    assert sat.read_text() == "record\n"
    assert fasta.read_text() == ">seq\nACGT\n"


def test_sibling_outputs_fail_together(tmp_path):
    """The .sat and .fasta must never disagree: either both appear or neither."""
    sat = tmp_path / "out.sat"
    fasta = tmp_path / "out.fasta"

    with pytest.raises(ValueError):
        with atomic_outputs(str(sat), str(fasta)) as (sat_fh, fasta_fh):
            sat_fh.write("record\n")
            fasta_fh.write(">seq\nACGT\n")
            raise ValueError("duplicate chromosome name")

    assert not sat.exists()
    assert not fasta.exists()
    assert list(tmp_path.iterdir()) == []


def test_optional_output_yields_none_handle(tmp_path):
    sat = tmp_path / "out.sat"

    with atomic_outputs(str(sat), None) as (sat_fh, fasta_fh):
        assert fasta_fh is None
        sat_fh.write("record\n")

    assert sat.read_text() == "record\n"


def test_publish_moves_external_binary_output(tmp_path):
    """The Rust bed-extract path writes the .partial itself; we publish it."""
    target = tmp_path / "out.sat"
    (tmp_path / f"out.sat{PARTIAL_SUFFIX}").write_text("written by bed-extract\n")

    publish(str(target))

    assert target.read_text() == "written by bed-extract\n"
    assert not (tmp_path / f"out.sat{PARTIAL_SUFFIX}").exists()


def test_publish_without_partial_raises(tmp_path):
    target = tmp_path / "missing.sat"
    with pytest.raises(FileNotFoundError, match="missing.sat.partial"):
        publish(str(target))
    assert not target.exists()


def test_discard_removes_partial_and_tolerates_absence(tmp_path):
    target = tmp_path / "out.sat"
    (tmp_path / f"out.sat{PARTIAL_SUFFIX}").write_text("junk")

    discard(str(target))
    assert not (tmp_path / f"out.sat{PARTIAL_SUFFIX}").exists()

    discard(str(target))  # idempotent, no raise


def test_partial_path_shape(tmp_path):
    assert partial_path("/a/b/c.sat") == f"/a/b/c.sat{PARTIAL_SUFFIX}"


def test_read_mode_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="for writing"):
        with atomic_output(str(tmp_path / "x"), mode="r"):
            pass


def test_creates_missing_parent_directory(tmp_path):
    target = tmp_path / "fastan" / "deep" / "out.tsv"

    with atomic_output(str(target)) as fh:
        fh.write("data\n")

    assert target.read_text() == "data\n"
    assert os.path.isdir(tmp_path / "fastan" / "deep")
