"""The subcommands, their exit codes, and that they do not disturb the old CLI."""

import os

import pytest

from satellome.compact import cli, record
from tests.unit.conftest_compact import build_run_dir, snapshot, write_stub_decomposer

PREFIX = "TESTASM_v1"


@pytest.fixture
def run_dir(tmp_path, monkeypatch):
    decomposer = write_stub_decomposer(tmp_path)
    monkeypatch.setenv("SATELLOME_ARRAYSPLITTER", decomposer)
    return build_run_dir(tmp_path, prefix=PREFIX, decomposer=decomposer)


def test_dispatch_ignores_the_pipeline_command_line():
    # 'satellome -i genome.fa -o out -t 8' must reach the flag parser untouched.
    assert cli.dispatch(["-i", "genome.fa", "-o", "out", "-t", "8"]) is None
    assert cli.dispatch([]) is None
    assert cli.dispatch(["--doctor"]) is None


def test_dispatch_claims_its_own_subcommands(run_dir):
    assert cli.dispatch(["compact", "--dry-run", run_dir]) == 0


def test_compact_then_expand_through_the_cli(run_dir):
    before = snapshot(run_dir)
    assert cli.dispatch(["compact", "--level", "3", run_dir]) == 0
    assert record.is_compacted(run_dir)
    assert cli.dispatch(["expand", run_dir]) == 0
    after = snapshot(run_dir)
    ignore = {record.RECORD_NAME, "run_manifest.json"}
    for rel, digest in before.items():
        if rel.endswith(".1aln") or rel in ignore:
            continue
        assert after[rel] == digest, rel


def test_a_refused_directory_exits_non_zero(tmp_path):
    empty = tmp_path / "not_a_run"
    empty.mkdir()
    assert cli.dispatch(["compact", str(empty)]) == 1


def test_expand_on_a_directory_that_was_never_compacted_exits_non_zero(run_dir):
    assert cli.dispatch(["expand", run_dir]) == 1


def test_missing_directory_argument_is_a_usage_error(capsys):
    assert cli.dispatch(["compact"]) == 2
    assert cli.dispatch(["expand"]) == 2


def test_explain_prints_the_policy_and_writes_nothing(capsys):
    assert cli.dispatch(["compact", "--explain"]) == 0
    out = capsys.readouterr().out
    assert ".monomers.tsv" in out and "primary_filtered" in out
    assert ".lengths" in out and "projection" in out
    # Every row carries its justification.
    assert "recomputable from the master" in out


def test_from_file_reads_a_directory_list(tmp_path, run_dir):
    listing = tmp_path / "dirs.txt"
    listing.write_text(f"# a comment\n{run_dir}\n\n")
    assert cli.dispatch(["compact", "--dry-run", "--from-file", str(listing)]) == 0


def test_from_file_that_does_not_exist_is_a_usage_error(tmp_path):
    assert cli.dispatch(
        ["compact", "--from-file", str(tmp_path / "nope.txt")]
    ) == 2


def test_check_recipes_exits_zero_when_every_recipe_reproduces(run_dir):
    assert cli.dispatch(["compact", "--check-recipes", run_dir]) == 0
    assert not record.is_compacted(run_dir)


def test_continue_on_error_processes_the_rest(tmp_path, run_dir):
    broken = tmp_path / "broken"
    broken.mkdir()
    code = cli.dispatch(
        ["compact", "--level", "3", "--continue-on-error", str(broken), run_dir]
    )
    assert code == 1  # one directory was refused
    assert record.is_compacted(run_dir)  # the other was still done


def test_min_array_length_is_a_parameter_not_a_constant(run_dir):
    from satellome.compact import columnar

    assert cli.dispatch(
        ["compact", "--level", "3", "--min-array-length", "5000", run_dir]
    ) == 0
    rec = record.load_record(run_dir)
    assert rec["config"]["min_array_length"] == 5000
    container = os.path.join(run_dir, "fastan", f"{PREFIX}.monomers.tsv.satz")
    kept = {line.split(b"\t")[0] for line in columnar.read_table(container)}
    # The 2400 bp array is above 1 kb but below 5 kb, so this threshold drops it.
    assert b"NC_073249.2_100_2500_2400_12" not in kept
    assert b"CM002639.2_900_12900_12000_24" in kept
