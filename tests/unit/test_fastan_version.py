"""Tests for FasTAN thread forwarding (BUG2) and the minimum-version gate.

BUG2: satellome never forwarded ``-t/--threads`` to FasTAN, so FasTAN ran at its
own default (8) regardless of ``-t``. build_fastan_command now forwards it as
``-T<n>``; we assert the concrete command string.

Version gate: before running, satellome must ensure FasTAN is at least
FASTAN_MIN_VERSION, rebuilding from the canonical repo when the installed binary
is older. We assert the decision (update / abort / proceed) for each branch
rather than just that a call succeeded.
"""

import json
import logging

import pytest

from satellome.installers import fastan as fastan_mod
from satellome.installers.fastan import (
    FASTAN_MIN_VERSION,
    _version_tuple,
    get_installed_fastan_version,
    ensure_fastan_version,
)


# --------------------------------------------------------------------------- #
# BUG2 — thread forwarding
# --------------------------------------------------------------------------- #

def test_build_fastan_command_forwards_threads():
    from satellome.main import build_fastan_command

    cmd = build_fastan_command("/opt/fastan", "/data/g.fna", "/out/g.1aln", 8)
    assert cmd == "/opt/fastan -T8 /data/g.fna /out/g.1aln"
    assert "-T8" in cmd


def test_build_fastan_command_can_pin_single_thread():
    """`-t 1` -> `-T1`, the confirmed-stable mode while the FasTAN race is open."""
    from satellome.main import build_fastan_command

    cmd = build_fastan_command("fastan", "g.fna", "g.1aln", 1)
    assert "-T1" in cmd
    assert cmd.split()[1] == "-T1"


# --------------------------------------------------------------------------- #
# _version_tuple
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("text,expected", [
    ("FasTAN 0.8 [satellome-fork:ad3002 edb1e5a]", (0, 8)),
    ("0.5", (0, 5)),
    ("FasTAN 0.10 [x]", (0, 10)),     # 0.10 > 0.8 numerically
    ("1.2.3", (1, 2, 3)),
    ("no version here", None),
    ("", None),
    (None, None),
])
def test_version_tuple(text, expected):
    assert _version_tuple(text) == expected


def test_version_ordering_is_numeric():
    assert _version_tuple("0.5") < _version_tuple(FASTAN_MIN_VERSION)
    assert _version_tuple("0.10") > _version_tuple("0.8")  # not lexicographic


# --------------------------------------------------------------------------- #
# get_installed_fastan_version
# --------------------------------------------------------------------------- #

class _FakeProc:
    def __init__(self, stdout=b"", returncode=0):
        self.stdout = stdout
        self.returncode = returncode


def test_version_read_from_binary_V_output(monkeypatch, tmp_path):
    binary = tmp_path / "fastan"
    binary.write_bytes(b"\x7fELF")

    monkeypatch.setattr(
        fastan_mod.subprocess, "run",
        lambda *a, **k: _FakeProc(b"FasTAN 0.8 [satellome-fork:ad3002 abc1234]\n", 0),
    )
    assert get_installed_fastan_version(str(binary)) == "0.8"


def test_version_falls_back_to_manifest_when_binary_unrunnable(monkeypatch, tmp_path):
    binary = tmp_path / "fastan"
    binary.write_bytes(b"\x7fELF")
    (tmp_path / "fastan.manifest.json").write_text(json.dumps({"source_version": "0.7"}))

    # Binary cannot be executed (e.g. cross-arch) -> -V raises.
    def _boom(*a, **k):
        raise OSError("Exec format error")
    monkeypatch.setattr(fastan_mod.subprocess, "run", _boom)

    assert get_installed_fastan_version(str(binary)) == "0.7"


def test_version_none_when_undeterminable(monkeypatch, tmp_path):
    binary = tmp_path / "fastan"
    binary.write_bytes(b"\x7fELF")  # no manifest beside it

    monkeypatch.setattr(
        fastan_mod.subprocess, "run",
        lambda *a, **k: _FakeProc(b"some unrelated output", 0),
    )
    assert get_installed_fastan_version(str(binary)) is None


# --------------------------------------------------------------------------- #
# ensure_fastan_version
# --------------------------------------------------------------------------- #

def _patch_install(monkeypatch, result):
    calls = {"n": 0}

    def _fake_install(force=False):
        calls["n"] += 1
        return result
    monkeypatch.setattr(fastan_mod, "install_fastan", _fake_install)
    return calls


def test_current_version_ok_does_not_reinstall(monkeypatch):
    monkeypatch.setattr(fastan_mod, "get_installed_fastan_version", lambda b: "0.9")
    calls = _patch_install(monkeypatch, True)

    ok, used_bin, version = ensure_fastan_version("/bin/fastan")
    assert ok is True
    assert used_bin == "/bin/fastan"
    assert version == "0.9"
    assert calls["n"] == 0  # no rebuild when already up to date


def test_old_version_triggers_update_and_succeeds(monkeypatch, tmp_path):
    managed = tmp_path / "bin"
    managed.mkdir()
    (managed / "fastan").write_bytes(b"new")
    monkeypatch.setattr(fastan_mod, "get_satellome_bin_dir", lambda: managed)

    versions = iter(["0.5", "0.8"])  # before update, then after
    monkeypatch.setattr(fastan_mod, "get_installed_fastan_version",
                        lambda b: next(versions))
    calls = _patch_install(monkeypatch, True)

    ok, used_bin, version = ensure_fastan_version("/usr/bin/fastan")
    assert ok is True
    assert calls["n"] == 1                       # rebuilt once
    assert used_bin == str(managed / "fastan")   # re-bound to managed binary
    assert version == "0.8"


def test_old_version_update_fails_is_visible_error(monkeypatch, caplog):
    monkeypatch.setattr(fastan_mod, "get_installed_fastan_version", lambda b: "0.5")
    _patch_install(monkeypatch, False)  # build/clone fails

    with caplog.at_level(logging.ERROR):
        ok, _bin, _ver = ensure_fastan_version("/usr/bin/fastan")
    assert ok is False
    assert any("update failed" in r.message.lower() for r in caplog.records)


def test_still_old_after_update_aborts(monkeypatch, tmp_path):
    managed = tmp_path / "bin"
    managed.mkdir()
    (managed / "fastan").write_bytes(b"old")
    monkeypatch.setattr(fastan_mod, "get_satellome_bin_dir", lambda: managed)

    # Canonical repo still ships an old version: stays 0.5 even after rebuild.
    monkeypatch.setattr(fastan_mod, "get_installed_fastan_version", lambda b: "0.5")
    _patch_install(monkeypatch, True)

    ok, _bin, version = ensure_fastan_version("/usr/bin/fastan")
    assert ok is False
    assert version == "0.5"
