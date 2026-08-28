"""Tests for the bundled Rust helper tool installer.

These tools had no install path at all before: absent from the wheel (they are
native per-platform binaries), absent from the sdist, and unknown to
--install-all. Two of them decide whether a run produces family-clustering and
telomere output, so "partially installed" must never read as success.
"""

import os
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


class TestPrebuiltDownload:
    """The fast path: binaries built by CI and attached to the release."""

    def _no_network(self, monkeypatch):
        monkeypatch.setattr(rust_tools, "_fetch", lambda url, timeout=None: None)

    def test_platform_suffix_maps_known_platforms(self, monkeypatch):
        import satellome.installers.base as base

        for (plat, arch), expected in [
            (("linux", "x86_64"), "linux-x86_64"),
            (("darwin", "arm64"), "macos-arm64"),
            (("darwin", "x86_64"), "macos-x86_64"),
        ]:
            monkeypatch.setattr(base, "detect_platform", lambda p=plat, a=arch: (p, a))
            assert rust_tools.platform_asset_suffix() == expected

    def test_unknown_platform_has_no_asset(self, monkeypatch):
        import satellome.installers.base as base

        monkeypatch.setattr(base, "detect_platform", lambda: ("linux", "aarch64"))
        assert rust_tools.platform_asset_suffix() is None

    def test_checksum_mismatch_refuses_to_install(self, tmp_path, monkeypatch, caplog):
        monkeypatch.setattr(rust_tools, "platform_asset_suffix", lambda: "linux-x86_64")
        monkeypatch.setattr(
            rust_tools, "_expected_checksums",
            lambda v, s: {"genome-size-linux-x86_64": "0" * 64},
        )
        monkeypatch.setattr(rust_tools, "_fetch", lambda url, timeout=None: b"tampered")

        with caplog.at_level("ERROR"):
            ok = rust_tools.download_prebuilt("genome-size", "1.11.0", tmp_path)

        assert ok is False
        assert not (tmp_path / "genome-size").exists(), (
            "a wrong binary is worse than a missing one - the pipeline would run it"
        )
        assert "Checksum mismatch" in caplog.text

    def test_matching_checksum_installs_and_is_executable(self, tmp_path, monkeypatch):
        import hashlib

        payload = b"\x7fELF fake binary"
        digest = hashlib.sha256(payload).hexdigest()
        monkeypatch.setattr(rust_tools, "platform_asset_suffix", lambda: "linux-x86_64")
        monkeypatch.setattr(
            rust_tools, "_expected_checksums",
            lambda v, s: {"genome-size-linux-x86_64": digest},
        )
        monkeypatch.setattr(rust_tools, "_fetch", lambda url, timeout=None: payload)

        assert rust_tools.download_prebuilt("genome-size", "1.11.0", tmp_path) is True

        installed = tmp_path / "genome-size"
        assert installed.read_bytes() == payload
        assert os.access(installed, os.X_OK)
        assert not (tmp_path / "genome-size.partial").exists(), "no leftover partial"

    def test_missing_asset_falls_back_rather_than_failing(self, tmp_path, monkeypatch):
        monkeypatch.setattr(rust_tools, "platform_asset_suffix", lambda: "linux-x86_64")
        self._no_network(monkeypatch)

        assert rust_tools.download_prebuilt("sat-family", "1.11.0", tmp_path) is False

    def test_checksum_file_is_parsed_from_sha256sum_format(self, monkeypatch):
        body = (
            b"aaaa  sat-family-linux-x86_64\n"
            b"bbbb *genome-size-linux-x86_64\n"
        )
        monkeypatch.setattr(rust_tools, "_fetch", lambda url, timeout=None: body)

        sums = rust_tools._expected_checksums("1.11.0", "linux-x86_64")

        assert sums["sat-family-linux-x86_64"] == "aaaa"
        assert sums["genome-size-linux-x86_64"] == "bbbb", "the '*' binary marker must be stripped"

    def test_install_prefers_download_and_skips_cargo_entirely(self, tmp_path, monkeypatch):
        bin_dir = tmp_path / "bin"
        bin_dir.mkdir()
        monkeypatch.delenv(rust_tools.ENV_NO_DOWNLOAD, raising=False)
        monkeypatch.setattr(rust_tools, "get_satellome_bin_dir", lambda: bin_dir)
        monkeypatch.setattr(rust_tools, "platform_asset_suffix", lambda: "linux-x86_64")
        monkeypatch.setattr(rust_tools, "download_prebuilt", lambda n, v, d: True)

        def fail():
            raise AssertionError("cargo must not be needed when downloads succeed")

        monkeypatch.setattr(rust_tools, "check_cargo", fail)

        assert rust_tools.install_rust_tools(force=True, tools=["find-gaps"]) is True

    def test_download_failure_falls_back_to_building(self, tmp_path, monkeypatch):
        bin_dir = tmp_path / "bin"
        bin_dir.mkdir()
        sources = tmp_path / "rust"
        (sources / "find-gaps").mkdir(parents=True)
        (sources / "find-gaps" / "Cargo.toml").write_text("[package]\n")
        monkeypatch.delenv(rust_tools.ENV_NO_DOWNLOAD, raising=False)
        monkeypatch.setattr(rust_tools, "get_satellome_bin_dir", lambda: bin_dir)
        monkeypatch.setattr(rust_tools, "platform_asset_suffix", lambda: "linux-x86_64")
        monkeypatch.setattr(rust_tools, "download_prebuilt", lambda n, v, d: False)
        monkeypatch.setattr(rust_tools, "check_cargo", lambda: (True, "cargo"))
        monkeypatch.setattr(rust_tools, "find_rust_sources", lambda: sources)
        built = []
        monkeypatch.setattr(
            rust_tools, "_build_one",
            lambda crate, name, target: built.append(name) or True,
        )

        assert rust_tools.install_rust_tools(force=True, tools=["find-gaps"]) is True
        assert built == ["find-gaps"], "must actually build what it could not download"

    def test_opt_out_skips_the_download_path(self, tmp_path, monkeypatch):
        bin_dir = tmp_path / "bin"
        bin_dir.mkdir()
        monkeypatch.setenv(rust_tools.ENV_NO_DOWNLOAD, "1")
        monkeypatch.setattr(rust_tools, "get_satellome_bin_dir", lambda: bin_dir)

        def fail(*a, **k):
            raise AssertionError("must not download when opted out")

        monkeypatch.setattr(rust_tools, "download_prebuilt", fail)
        monkeypatch.setattr(rust_tools, "check_cargo", lambda: (False, "no cargo"))

        assert rust_tools.install_rust_tools(force=True, tools=["find-gaps"]) is False
