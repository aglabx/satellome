"""Tests for re-running individual steps on a finished output directory.

The contract: a re-run amends an existing analysis. It must reuse the original
run's parameters (not the current command line), refuse clearly when the
directory is not a finished run, and leave the manifest an accurate record —
otherwise `--verify-run` starts lying about what the directory contains.
"""

import json
import os

import pytest

from satellome.core_functions.tools.rerun import (
    RERUNNABLE,
    RERUNNABLE_BY_NAME,
    RerunError,
    check_requirements,
    format_step_list,
    load_run_context,
    merge_steps,
    parse_step_names,
    record_rerun,
)
from satellome.core_functions.tools.run_manifest import MANIFEST_NAME


def write_manifest(directory, **overrides):
    manifest = {
        "manifest_version": 1,
        "satellome_version": "1.13.0",
        "project": "demo",
        "argv": ["-i", "/data/genome.fasta", "-o", str(directory), "-t", "8"],
        "input_fasta": "/data/genome.fasta",
        "taxon_name": "Homo sapiens",
        "taxid": "9606",
        "genome_size": 3_000_000_000,
        "steps": {"classification": "ok", "drawing": "ok"},
        "files": [],
        "complete": True,
        "host": "h", "platform": "p", "pid": 1,
        "leftover_partials": [],
    }
    manifest.update(overrides)
    (directory / MANIFEST_NAME).write_text(json.dumps(manifest))
    return manifest


class TestStepNames:
    def test_single_and_multiple(self):
        assert parse_step_names("ucsc_track") == ["ucsc_track"]
        assert set(parse_step_names("drawing,ucsc_track")) == {"drawing", "ucsc_track"}

    def test_order_follows_the_pipeline_not_the_command_line(self):
        """Two steps must apply in the order a full run would apply them."""
        pipeline_order = [s.name for s in RERUNNABLE]
        got = parse_step_names("ucsc_track,classification")
        assert got == [n for n in pipeline_order if n in got]
        assert got[0] == "classification"

    def test_whitespace_and_empties_tolerated(self):
        assert parse_step_names(" drawing , ,ucsc_track ") == ["drawing", "ucsc_track"]

    def test_unknown_step_is_named_with_the_alternatives(self):
        with pytest.raises(RerunError, match="unknown step"):
            parse_step_names("trf_search")
        try:
            parse_step_names("nope")
        except RerunError as error:
            assert "ucsc_track" in str(error), "must list what is accepted"

    def test_empty_is_rejected(self):
        with pytest.raises(RerunError, match="at least one step"):
            parse_step_names(" , ")

    def test_expensive_search_is_deliberately_not_rerunnable(self):
        """Hiding hours of recomputation behind --rerun would be a trap."""
        for name in ("trf_search", "fastan", "search"):
            assert name not in RERUNNABLE_BY_NAME

    def test_step_list_is_printable_and_mentions_force(self):
        text = "\n".join(format_step_list())
        for spec in RERUNNABLE:
            assert spec.name in text and spec.description in text
        assert "--force" in text


class TestContext:
    def test_reads_parameters_from_the_manifest(self, tmp_path):
        write_manifest(tmp_path)
        context = load_run_context(str(tmp_path))
        assert context["project"] == "demo"
        assert context["input_fasta"] == "/data/genome.fasta"
        assert context["genome_size"] == 3_000_000_000
        assert context["steps"] == {"classification": "ok", "drawing": "ok"}
        assert context["argv"][0] == "-i"

    def test_missing_directory(self, tmp_path):
        with pytest.raises(RerunError, match="not a directory"):
            load_run_context(str(tmp_path / "nope"))

    def test_directory_without_a_manifest_says_what_to_do(self, tmp_path):
        with pytest.raises(RerunError, match="run the full pipeline first"):
            load_run_context(str(tmp_path))

    def test_corrupt_manifest_is_not_guessed_around(self, tmp_path):
        (tmp_path / MANIFEST_NAME).write_text("{not json")
        with pytest.raises(RerunError):
            load_run_context(str(tmp_path))

    def test_manifest_without_argv_is_refused(self, tmp_path):
        """Guessing the options would amend the run with different settings."""
        write_manifest(tmp_path, argv=None)
        with pytest.raises(RerunError, match="original command line"):
            load_run_context(str(tmp_path))

    def test_manifest_without_project_is_refused(self, tmp_path):
        write_manifest(tmp_path, project="")
        with pytest.raises(RerunError, match="does not name a project"):
            load_run_context(str(tmp_path))


class TestRequirements:
    def test_missing_sat_file_is_reported_before_any_work(self, tmp_path):
        settings = {"trf_file": str(tmp_path / "absent.sat"),
                    "fasta_file": str(tmp_path / "absent.fa")}
        with pytest.raises(RerunError, match="tandem repeat file"):
            check_requirements(["ucsc_track"], settings)

    def test_ucsc_track_also_needs_the_genome(self, tmp_path):
        sat = tmp_path / "p.sat"
        sat.write_text("x\n")
        settings = {"trf_file": str(sat), "fasta_file": str(tmp_path / "absent.fa")}
        with pytest.raises(RerunError, match="genome FASTA"):
            check_requirements(["ucsc_track"], settings)

    def test_satisfied_requirements_pass(self, tmp_path):
        sat = tmp_path / "p.sat"
        fasta = tmp_path / "g.fa"
        sat.write_text("x\n")
        fasta.write_text(">c\nACGT\n")
        check_requirements(
            ["ucsc_track", "drawing"],
            {"trf_file": str(sat), "fasta_file": str(fasta)},
        )

    def test_every_step_declares_its_requirements(self, tmp_path):
        settings = {"trf_file": str(tmp_path / "a"), "fasta_file": str(tmp_path / "b")}
        for spec in RERUNNABLE:
            required = spec.requires(settings)
            assert required, f"{spec.name} declares no inputs"


class TestManifestUpdate:
    def test_rerun_steps_replace_and_others_are_kept(self):
        merged = merge_steps(
            {"classification": "ok", "drawing": "failed", "sat_family": "ok"},
            {"drawing": "ok", "ucsc_track": "ok"},
        )
        assert merged == {
            "classification": "ok",   # untouched step keeps its status
            "drawing": "ok",          # re-run replaces the old failure
            "sat_family": "ok",
            "ucsc_track": "ok",       # new step appears
        }

    def test_a_failed_rerun_overwrites_a_previous_success(self):
        """Otherwise the manifest would keep claiming a step that now fails."""
        merged = merge_steps({"drawing": "ok"}, {"drawing": "failed"})
        assert merged["drawing"] == "failed"

    def test_history_accumulates(self):
        manifest = {}
        first = record_rerun(manifest, ["ucsc_track"], "1.13.0")
        assert len(first) == 1
        manifest["reruns"] = first
        second = record_rerun(manifest, ["drawing"], "1.13.0")
        assert len(second) == 2
        assert second[0]["steps"] == ["ucsc_track"]
        assert second[1]["steps"] == ["drawing"]
        assert all("at" in entry for entry in second)
