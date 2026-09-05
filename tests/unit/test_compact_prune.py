"""A default run must not leave behind what compaction would remove anyway."""

import os

import pytest

from satellome.compact.prune import prune_run_directory
from tests.unit.conftest_compact import build_run_dir, write_stub_decomposer

PREFIX = "TESTASM_v1"


@pytest.fixture
def run_dir(tmp_path, monkeypatch):
    decomposer = write_stub_decomposer(tmp_path)
    monkeypatch.setenv("SATELLOME_ARRAYSPLITTER", decomposer)
    return build_run_dir(tmp_path, prefix=PREFIX, decomposer=decomposer)


def rel_files(run_dir):
    return {
        os.path.relpath(os.path.join(root, name), run_dir)
        for root, _dirs, files in os.walk(run_dir)
        for name in files
    }


def test_prune_removes_the_derived_views_and_the_gff(run_dir):
    removed, freed = prune_run_directory(run_dir)
    remaining = rel_files(run_dir)
    assert freed > 0
    for name in [f"{PREFIX}.micro.sat.gz", f"{PREFIX}.pmicro.sat.gz",
                 f"{PREFIX}.tssr.sat.gz", f"{PREFIX}.complex.sat.gz",
                 f"{PREFIX}.100kb.sat.gz", f"{PREFIX}.1000kb.sat.gz"]:
        assert name not in remaining
    assert not any(rel.startswith("gff3/") for rel in remaining)
    assert f"fastan/{PREFIX}.lengths" not in remaining
    assert f"fastan/{PREFIX}.1aln" not in remaining


def test_prune_keeps_what_the_pipeline_itself_reads(run_dir):
    prune_run_directory(run_dir)
    remaining = rel_files(run_dir)
    # --large_file defaults to 1kb; drawing and annotation read that file, and
    # --rerun drawing reads it again later.
    assert f"{PREFIX}.1kb.sat.gz" in remaining
    assert f"{PREFIX}.10kb.sat.gz" in remaining


def test_prune_keeps_every_primary_output(run_dir):
    prune_run_directory(run_dir)
    remaining = rel_files(run_dir)
    for name in [f"{PREFIX}.sat.gz", f"fastan/{PREFIX}.bed",
                 f"fastan/{PREFIX}.monomers.tsv.gz", f"fastan/{PREFIX}.hors.tsv.gz",
                 f"fastan/{PREFIX}.decomposed.fasta.gz",
                 f"fastan/{PREFIX}.summary.tsv.gz", "results.yaml"]:
        assert name in remaining, name


def test_prune_leaves_unclassified_files_alone(run_dir):
    stray = os.path.join(run_dir, "notes.parquet")
    open(stray, "wb").write(b"opaque\n")
    prune_run_directory(run_dir)
    assert os.path.exists(stray)


def test_dry_run_reports_without_removing(run_dir):
    before = rel_files(run_dir)
    removed, freed = prune_run_directory(run_dir, dry_run=True)
    assert removed and freed > 0
    assert rel_files(run_dir) == before


def test_prune_is_idempotent(run_dir):
    prune_run_directory(run_dir)
    after_first = rel_files(run_dir)
    removed, freed = prune_run_directory(run_dir)
    assert removed == [] and freed == 0
    assert rel_files(run_dir) == after_first


def test_prune_names_what_it_removed(run_dir, caplog):
    import logging

    with caplog.at_level(logging.INFO, logger="satellome"):
        prune_run_directory(run_dir)
    text = caplog.text
    assert "Compact output" in text
    assert "--extended-output" in text
    assert "gff" in text and "lengths" in text


def test_the_flag_exists_and_defaults_to_pruning():
    """A default run prunes; --extended-output is the opt-out, not the default."""
    from satellome.main import build_parser

    parser = build_parser()
    assert parser.parse_args([]).extended_output is False
    assert parser.parse_args(["--extended-output"]).extended_output is True


def test_a_recorded_argv_still_says_whether_to_prune_on_a_rerun():
    """--rerun recovers the flag from the run's own argv, not the current one."""
    from satellome.main import build_parser

    extended = vars(build_parser().parse_args(["--extended-output", "-t", "8"]))
    plain = vars(build_parser().parse_args(["-t", "8"]))
    assert extended["extended_output"] is True
    assert plain["extended_output"] is False
