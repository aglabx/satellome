"""Tests for the bundled Rust helper tool installer.

These tools had no install path at all before: absent from the wheel (they are
native per-platform binaries), absent from the sdist, and unknown to
--install-all. Two of them decide whether a run produces family-clustering and
telomere output, so "partially installed" must never read as success.
"""

import subprocess
from pathlib import Path

import pytest

from satellome.installers import rust_tools
from satellome.installers.rust_tools import (
    RUST_TOOLS,
    check_cargo,
    find_rust_sources,
    install_rust_tools,
)


class TestSpecCoverage:
    def test_every_rust_tool_has_a_crate_in_the_repo(self):
        sources = find_rust_sources()
        if sources is None:
            pytest.skip("no local rust/ sources (wheel install)")
        for name in RUST_TOOLS:
            assert (sources / name / "Cargo.toml").is_file(), f"no crate for {name}"

    def test_rust_tools_match_the_doctor_registry(self):
        """A tool the doctor reports must be one the installer can provide."""
        from satellome.core_functions.tools.env_check import TOOL_SPECS, RUST_TOOLS_INSTALL

        advertised = {s.name for s in TOOL_SPECS if s.install == RUST_TOOLS_INSTALL}
        assert advertised == set(RUST_TOOLS), (
            "the doctor promises --install-rust-tools for tools it cannot install"
        )


class TestCargoCheck:
    def test_missing_cargo_explains_how_to_get_it(self, monkeypatch):
        monkeypatch.setattr(rust_tools.shutil, "which", lambda name: None)
        ok, message = check_cargo()
        assert ok is False
        assert "rustup" in message
        assert "satellome --install-rust-tools" in message

    def test_install_aborts_with_a_reason_when_cargo_is_absent(self, tmp_path, monkeypatch, caplog):
        monkeypatch.setattr(rust_tools, "get_satellome_bin_dir", lambda: tmp_path)
        monkeypatch.setattr(rust_tools.shutil, "which", lambda name: None)

        with caplog.at_level("ERROR"):
            result = install_rust_tools(force=True, tools=["genome-size"])

        assert result is False
        assert "rustup" in caplog.text


class TestInstallOutcome:
    def _fake_sources(self, tmp_path, names):
        sources = tmp_path / "rust"
        for name in names:
            (sources / name).mkdir(parents=True)
            (sources / name / "Cargo.toml").write_text("[package]\n")
        return sources

    def test_partial_success_is_reported_as_failure(self, tmp_path, monkeypatch, caplog):
        bin_dir = tmp_path / "bin"
        bin_dir.mkdir()
        sources = self._fake_sources(tmp_path, ["sat-family", "find-gaps"])
        monkeypatch.setattr(rust_tools, "get_satellome_bin_dir", lambda: bin_dir)
        monkeypatch.setattr(rust_tools, "find_rust_sources", lambda: sources)
        monkeypatch.setattr(rust_tools, "check_cargo", lambda: (True, "cargo"))

        def build(crate_dir, name, target_dir):
            return name == "find-gaps"

        monkeypatch.setattr(rust_tools, "_build_one", build)

        with caplog.at_level("ERROR"):
            result = install_rust_tools(force=True, tools=["sat-family", "find-gaps"])

        assert result is False, "three-of-five must not be reported as success"
        assert "Still missing after installation: sat-family" in caplog.text

    def test_all_built_is_success(self, tmp_path, monkeypatch):
        bin_dir = tmp_path / "bin"
        bin_dir.mkdir()
        sources = self._fake_sources(tmp_path, ["find-gaps"])
        monkeypatch.setattr(rust_tools, "get_satellome_bin_dir", lambda: bin_dir)
        monkeypatch.setattr(rust_tools, "find_rust_sources", lambda: sources)
        monkeypatch.setattr(rust_tools, "check_cargo", lambda: (True, "cargo"))
        monkeypatch.setattr(rust_tools, "_build_one", lambda *a: True)

        assert install_rust_tools(force=True, tools=["find-gaps"]) is True

    def test_missing_crate_is_named_not_skipped(self, tmp_path, monkeypatch, caplog):
        bin_dir = tmp_path / "bin"
        bin_dir.mkdir()
        sources = self._fake_sources(tmp_path, ["find-gaps"])
        monkeypatch.setattr(rust_tools, "get_satellome_bin_dir", lambda: bin_dir)
        monkeypatch.setattr(rust_tools, "find_rust_sources", lambda: sources)
        monkeypatch.setattr(rust_tools, "check_cargo", lambda: (True, "cargo"))
        monkeypatch.setattr(rust_tools, "_build_one", lambda *a: True)

        with caplog.at_level("ERROR"):
            result = install_rust_tools(force=True, tools=["sat-family"])

        assert result is False
        assert "No crate for sat-family" in caplog.text

    def test_already_installed_short_circuits_without_building(self, tmp_path, monkeypatch):
        bin_dir = tmp_path / "bin"
        bin_dir.mkdir()
        existing = bin_dir / "find-gaps"
        existing.write_text("#!/bin/sh\n")
        existing.chmod(0o755)
        monkeypatch.setattr(rust_tools, "get_satellome_bin_dir", lambda: bin_dir)

        def fail(*args, **kwargs):
            raise AssertionError("must not rebuild an installed tool")

        monkeypatch.setattr(rust_tools, "check_cargo", fail)

        assert install_rust_tools(force=False, tools=["find-gaps"]) is True

    def test_build_failure_surfaces_the_compiler_output(self, tmp_path, monkeypatch, caplog):
        crate = tmp_path / "sat-family"
        crate.mkdir()
        (crate / "Cargo.toml").write_text("[package]\n")

        def fake_run(cmd, **kwargs):
            return subprocess.CompletedProcess(
                cmd, 1, stdout=b"", stderr=b"error[E0432]: unresolved import `foo`"
            )

        monkeypatch.setattr(rust_tools.subprocess, "run", fake_run)

        with caplog.at_level("ERROR"):
            ok = rust_tools._build_one(crate, "sat-family", tmp_path)

        assert ok is False
        assert "E0432" in caplog.text, "a bare 'build failed' leaves nothing to act on"

    def test_binary_missing_after_a_clean_build_is_an_error(self, tmp_path, monkeypatch, caplog):
        crate = tmp_path / "find-gaps"
        crate.mkdir()

        def fake_run(cmd, **kwargs):
            return subprocess.CompletedProcess(cmd, 0, stdout=b"", stderr=b"")

        monkeypatch.setattr(rust_tools.subprocess, "run", fake_run)

        with caplog.at_level("ERROR"):
            ok = rust_tools._build_one(crate, "find-gaps", tmp_path)

        assert ok is False
        assert "no binary at" in caplog.text
