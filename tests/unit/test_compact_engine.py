"""End to end: compact, expand, and the promises made in between.

Every assertion here is about *content* — the rows and bytes a consumer would
read — not about exit codes. A compaction that returns 0 and loses a column is
the failure mode this whole feature exists to prevent.
"""

import gzip
import json
import os

import pytest

from satellome.compact import engine, readiness, record
from satellome.compact.policy import CompactConfig
from tests.unit.conftest_compact import (
    build_run_dir,
    content_md5,
    snapshot,
    write_stub_decomposer,
)

PREFIX = "TESTASM_v1"


@pytest.fixture
def decomposer(tmp_path, monkeypatch):
    path = write_stub_decomposer(tmp_path)
    monkeypatch.setenv("SATELLOME_ARRAYSPLITTER", path)
    return path


@pytest.fixture
def run_dir(tmp_path, decomposer):
    return build_run_dir(tmp_path, prefix=PREFIX, decomposer=decomposer)


@pytest.fixture
def config():
    return CompactConfig(min_array_length=1000, level=3, verify_drops=True)


def test_the_fixture_is_a_complete_run(run_dir):
    ready = readiness.check(run_dir)
    assert ready.ok, ready.reasons
    assert ready.prefix == PREFIX


def test_compact_then_expand_restores_every_file_content_identically(
    run_dir, config, decomposer
):
    before = snapshot(run_dir)
    assert engine.compact(run_dir, config).ok
    assert engine.expand(run_dir).ok
    after = snapshot(run_dir)

    # The alignment intermediate is dropped by decision and cannot be rebuilt
    # without the genome; everything else must come back exactly.
    # run_manifest.json is deliberately not restored byte for byte: it gains a
    # compaction_history entry, which is the point of keeping --verify-run
    # meaningful. Its own tests cover that.
    ignore = {record.RECORD_NAME, "run_manifest.json"}
    expected = {k: v for k, v in before.items()
                if not k.endswith(".1aln") and k not in ignore}
    restored = {k: v for k, v in after.items() if k not in ignore}
    assert set(restored) == set(expected)
    differing = {k for k in expected if restored[k] != expected[k]}
    assert not differing, f"content changed for {sorted(differing)}"


def test_expand_reports_the_one_file_it_cannot_rebuild(run_dir, config, decomposer):
    engine.compact(run_dir, config)
    outcome = engine.expand(run_dir)
    assert any(".1aln" in line and "not restorable" in line for line in outcome.lines)
    assert not os.path.exists(os.path.join(run_dir, "fastan", f"{PREFIX}.1aln"))


def test_the_derived_views_are_actually_gone_after_compaction(run_dir, config):
    engine.compact(run_dir, config)
    for name in [f"{PREFIX}.micro.sat.gz", f"{PREFIX}.1kb.sat.gz",
                 f"{PREFIX}.100kb.sat.gz", f"{PREFIX}.tssr.sat.gz"]:
        assert not os.path.exists(os.path.join(run_dir, name)), name
    assert not os.path.isdir(os.path.join(run_dir, "gff3")) or not os.listdir(
        os.path.join(run_dir, "gff3")
    )
    assert not os.path.exists(os.path.join(run_dir, "fastan", f"{PREFIX}.lengths"))


def test_the_master_survives_as_a_readable_table(run_dir, config):
    from satellome.compact import columnar

    original = content_md5(os.path.join(run_dir, f"{PREFIX}.sat.gz"))
    engine.compact(run_dir, config)
    container = os.path.join(run_dir, f"{PREFIX}.sat.satz")
    assert os.path.exists(container)
    assert columnar.read_footer(container)["md5"] == original
    rows = list(columnar.read_table(container))
    # Every array of the run, with the header row stored separately.
    assert len(rows) == 405
    assert rows[0].split(b"\t")[2] == b"NC_073249.2"
    assert rows[0].split(b"\t")[14] == b"2400"


def test_the_per_copy_cut_keeps_rows_of_long_arrays_and_drops_the_rest(
    run_dir, config
):
    from satellome.compact import columnar

    engine.compact(run_dir, config)
    container = os.path.join(run_dir, "fastan", f"{PREFIX}.monomers.tsv.satz")
    kept_ids = {line.split(b"\t")[0] for line in columnar.read_table(container)}
    # 2400, 1600 and 12000 bp arrays stay; 480 and 300 bp arrays go.
    assert b"NC_073249.2_100_2500_2400_12" in kept_ids
    assert b"CM002639.2_900_12900_12000_24" in kept_ids
    assert b"NC_073249.2_5000_5480_480_8" not in kept_ids
    assert b"CM002639.2_77_377_300_10" not in kept_ids


def test_sub_threshold_arrays_keep_their_master_row(run_dir, config):
    from satellome.compact import columnar

    engine.compact(run_dir, config)
    rows = list(columnar.read_table(os.path.join(run_dir, f"{PREFIX}.sat.satz")))
    lengths = {int(row.split(b"\t")[14]) for row in rows}
    assert {300, 480} <= lengths
    # And the array sequence is still there, which is what makes the dropped
    # copy rows recomputable at all.
    short = [r for r in rows if int(r.split(b"\t")[14]) == 300][0]
    assert len(short.split(b"\t")[11]) == 300


def test_compaction_is_a_no_op_the_second_time(run_dir, config):
    engine.compact(run_dir, config)
    before = snapshot(run_dir)
    outcome = engine.compact(run_dir, config)
    assert outcome.ok and outcome.status == "noop"
    assert snapshot(run_dir) == before


def test_a_truncated_master_is_refused_and_nothing_is_touched(run_dir, config):
    master = os.path.join(run_dir, f"{PREFIX}.sat.gz")
    data = open(master, "rb").read()
    open(master, "wb").write(data[: len(data) // 2])
    before = sorted(os.listdir(run_dir))

    outcome = engine.compact(run_dir, config)
    assert not outcome.ok and outcome.status == "refused"
    assert any("does not decompress to the last byte" in line for line in outcome.lines)
    assert sorted(os.listdir(run_dir)) == before
    assert not record.is_compacted(run_dir)


def test_a_staging_directory_is_refused_by_construction(tmp_path, decomposer, config):
    staging = tmp_path / "_incoming"
    staging.mkdir()
    run_dir = build_run_dir(staging, prefix=PREFIX, decomposer=decomposer)
    outcome = engine.compact(run_dir, config)
    assert not outcome.ok
    assert any("staging" in line for line in outcome.lines)


def test_an_in_flight_run_is_refused(run_dir, config):
    # Compression happens at the end, so an uncompressed .sat beside its .gz is
    # the signature of a worker still writing here.
    open(os.path.join(run_dir, f"{PREFIX}.sat"), "w").write("partial\n")
    outcome = engine.compact(run_dir, config)
    assert not outcome.ok
    assert any("sits uncompressed beside" in line for line in outcome.lines)


def test_a_leftover_fasta_directory_is_not_mistaken_for_an_active_run(run_dir, config):
    # 1,015 finished directories carry one. Refusing them would be wrong.
    assert os.path.isdir(os.path.join(run_dir, "fasta"))
    assert readiness.check(run_dir).ok


def test_a_missing_rc_complete_marker_does_not_block_compaction(run_dir, config):
    # 1,020 finished catalogues predate the marker.
    assert not os.path.exists(os.path.join(run_dir, ".rc_complete"))
    assert engine.compact(run_dir, config).ok


def test_the_record_carries_a_digest_and_a_recipe_for_everything_it_removed(
    run_dir, config
):
    engine.compact(run_dir, config)
    rec = record.load_record(run_dir)
    dropped = [e for e in rec["files"] if e["action"] == "drop"]
    assert dropped
    for entry in dropped:
        assert entry["before"]["md5"], entry["path"]
        assert entry["recipe"]["kind"], entry["path"]
    reencoded = [e for e in rec["files"] if e["action"].startswith("reencode")]
    for entry in reencoded:
        assert entry["before"]["md5"] and entry["after"]["md5"], entry["path"]
    assert rec["totals"]["freed"] > 0


def test_expand_refuses_a_directory_that_was_never_compacted(run_dir):
    outcome = engine.expand(run_dir)
    assert not outcome.ok
    assert any(record.RECORD_NAME in line for line in outcome.lines)


def test_expand_fails_loudly_when_a_container_is_missing(run_dir, config):
    engine.compact(run_dir, config)
    os.unlink(os.path.join(run_dir, "fastan", f"{PREFIX}.bed.satz"))
    outcome = engine.expand(run_dir)
    assert not outcome.ok
    assert any("bed" in line and "gone" in line for line in outcome.lines)
    # The record stays, so the directory is not silently left half-restored.
    assert record.is_compacted(run_dir)


def test_expand_detects_a_decomposer_that_no_longer_reproduces_the_rows(
    run_dir, config, tmp_path, monkeypatch
):
    engine.compact(run_dir, config)
    # Swap in a decomposer that produces different rows: this is the version
    # drift the corpus is exposed to, and it must surface as a failure rather
    # than as a plausible-looking file.
    other = tmp_path / "wrong_arraysplitter"
    other.write_text(
        "#!/bin/sh\n"
        'out=""\n'
        'while [ $# -gt 0 ]; do if [ "$1" = "-o" ]; then out="$2"; fi; shift; done\n'
        'printf "array_id\\ttype\\tidx\\tlength\\tsequence\\n" > "$out.monomers.tsv"\n'
        'printf "array_id\\ttype\\tidx\\tlength\\tsequence\\n" > "$out.hors.tsv"\n'
        ': > "$out.decomposed.fasta"\n'
    )
    other.chmod(0o755)
    monkeypatch.setenv("SATELLOME_ARRAYSPLITTER", str(other))

    outcome = engine.expand(run_dir)
    assert not outcome.ok
    assert any("no longer reproduces" in line for line in outcome.lines)


def test_a_cut_that_cannot_be_proven_is_not_made_but_the_file_is_still_encoded(
    run_dir, config, tmp_path, monkeypatch
):
    from satellome.compact import columnar

    broken = tmp_path / "broken_arraysplitter"
    broken.write_text("#!/bin/sh\nexit 3\n")
    broken.chmod(0o755)
    monkeypatch.setenv("SATELLOME_ARRAYSPLITTER", str(broken))

    before = content_md5(os.path.join(run_dir, "fastan", f"{PREFIX}.monomers.tsv.gz"))
    outcome = engine.compact(run_dir, config)
    assert outcome.ok

    # No rows were dropped, because dropping them could not be shown reversible.
    container = os.path.join(run_dir, "fastan", f"{PREFIX}.monomers.tsv.satz")
    assert columnar.read_footer(container)["md5"] == before
    kept_ids = {line.split(b"\t")[0] for line in columnar.read_table(container)}
    assert b"NC_073249.2_5000_5480_480_8" in kept_ids

    # And it is recorded as an unfiltered re-encode, with the reason named.
    entry = [e for e in record.load_record(run_dir)["files"]
             if e["kind"] == "monomers"][0]
    assert entry["action"] == "reencode"
    assert entry["filter"] is None
    assert any("monomers" in note for note in outcome.kept_back)


def test_expand_restores_a_file_whose_cut_was_not_made(
    run_dir, config, tmp_path, monkeypatch
):
    broken = tmp_path / "broken_arraysplitter"
    broken.write_text("#!/bin/sh\nexit 3\n")
    broken.chmod(0o755)
    monkeypatch.setenv("SATELLOME_ARRAYSPLITTER", str(broken))

    before = content_md5(os.path.join(run_dir, "fastan", f"{PREFIX}.monomers.tsv.gz"))
    engine.compact(run_dir, config)
    assert engine.expand(run_dir).ok
    restored = os.path.join(run_dir, "fastan", f"{PREFIX}.monomers.tsv.gz")
    assert content_md5(restored) == before


def test_an_unclassified_file_is_kept_and_named(run_dir, config):
    stray = os.path.join(run_dir, "notes.parquet")
    open(stray, "wb").write(b"opaque\n")
    outcome = engine.compact(run_dir, config)
    assert os.path.exists(stray)
    rec = record.load_record(run_dir)
    assert "notes.parquet" in rec["kept_unclassified"]
    assert any("unclassified" in line for line in outcome.lines)


def test_the_dry_run_ledger_predicts_the_realised_saving_within_ten_percent(
    run_dir, config
):
    planned = engine.compact(run_dir, config, dry_run=True)
    assert planned.status == "planned"
    # Nothing was written.
    assert not record.is_compacted(run_dir)

    actual = engine.compact(run_dir, config)
    predicted_freed = planned.before - planned.after
    realised_freed = actual.before - actual.after
    assert realised_freed > 0
    assert abs(predicted_freed - realised_freed) <= 0.10 * max(
        predicted_freed, realised_freed
    ), f"predicted {predicted_freed}, realised {realised_freed}"


def test_verify_run_still_passes_on_a_compacted_directory(run_dir, config):
    from satellome.core_functions.tools.run_manifest import verify_run

    assert verify_run(run_dir).ok
    engine.compact(run_dir, config)
    result = verify_run(run_dir)
    assert result.ok, result.problems
    manifest = json.load(open(os.path.join(run_dir, "run_manifest.json")))
    history = manifest["compaction_history"]
    assert history[-1]["action"] == "compact"
    assert history[-1]["before_bytes"] > history[-1]["after_bytes"]


def test_verify_run_still_passes_after_expansion(run_dir, config):
    from satellome.core_functions.tools.run_manifest import verify_run

    engine.compact(run_dir, config)
    engine.expand(run_dir)
    result = verify_run(run_dir)
    assert result.ok, result.problems
    manifest = json.load(open(os.path.join(run_dir, "run_manifest.json")))
    assert manifest["compaction_history"][-1]["action"] == "expand"


def test_a_directory_without_a_manifest_is_left_alone(tmp_path, decomposer, config):
    run_dir = build_run_dir(tmp_path, prefix=PREFIX, decomposer=decomposer,
                            with_manifest=False)
    assert engine.compact(run_dir, config).ok
    assert not os.path.exists(os.path.join(run_dir, "run_manifest.json"))


def test_uncompressed_tables_compact_the_same_way(tmp_path, decomposer, config):
    # 1,012 assemblies carry the per-copy kinds uncompressed; both forms must work.
    run_dir = build_run_dir(tmp_path, prefix=PREFIX, decomposer=decomposer,
                            gzip_tables=False)
    before = snapshot(run_dir)
    assert engine.compact(run_dir, config).ok
    assert engine.expand(run_dir).ok
    after = snapshot(run_dir)
    for rel, digest in before.items():
        if rel.endswith(".1aln") or rel == "run_manifest.json":
            continue
        assert after[rel] == digest, rel


def test_check_recipes_reports_each_kind_and_changes_nothing(run_dir, config):
    before = snapshot(run_dir)
    outcome = engine.check_recipes(run_dir, config)
    assert outcome.ok, outcome.lines
    assert snapshot(run_dir) == before
    text = "\n".join(outcome.lines)
    assert "sat_micro" in text and "lengths" in text and "gff" in text
    assert "monomers" in text
    assert "checks passed" in text


def test_check_recipes_fails_when_a_recipe_stops_reproducing(run_dir, config):
    # Corrupt one derived view so the classifier's output no longer matches it.
    path = os.path.join(run_dir, f"{PREFIX}.micro.sat.gz")
    with gzip.open(path, "wb") as fh:
        fh.write(b"project\ttrf_id\n")
    outcome = engine.check_recipes(run_dir, config)
    assert not outcome.ok
    assert any("FAIL" in line and "sat_micro" in line for line in outcome.lines)


def test_no_verify_drops_still_records_a_recipe_and_a_digest(run_dir):
    config = CompactConfig(min_array_length=1000, level=3, verify_drops=False)
    engine.compact(run_dir, config)
    rec = record.load_record(run_dir)
    dropped = [e for e in rec["files"] if e["action"] == "drop"]
    assert dropped
    assert all(e["before"]["md5"] for e in dropped)
    # And it says plainly that it did not check them.
    unverified = [e for e in dropped
                  if e["recipe"].get("verified_at_compaction") is False]
    assert unverified
