"""Tests for binary provenance/integrity manifests.

These guard the "silent stale/corrupt binary" failure mode: an installed binary
that no longer matches its recorded build must produce a *visible* error
(status "mismatch"), never a silent fallback. The assertions check the concrete
regression path (the message names the discrepancy), not merely truthiness.
"""

import json

from satellome.installers import base
from satellome.installers.base import (
    write_binary_manifest,
    verify_binary_manifest,
    manifest_path_for,
    resolve_binary,
    read_build_signature,
    compute_sha256,
)

MARKER = b"FASTAN_SATELLOME_SIGNATURE=ad3002:deadbee:0.5"


def _make_fake_binary(directory, name="fastan", body=b"\x7fELF fake binary "):
    """Create a non-executable file that embeds the greppable signature marker.

    The marker is NUL-terminated, mirroring a real C string constant in the
    binary image, so read_build_signature stops at the terminator.
    """
    path = directory / name
    path.write_bytes(body + MARKER + b"\x00 trailing bytes")
    return path


def test_manifest_roundtrip_verifies_ok(tmp_path):
    binary = _make_fake_binary(tmp_path)
    manifest = write_binary_manifest(
        binary, tool="fastan", repo="https://example/repo.git", git_sha="abc123"
    )
    assert manifest is not None
    assert manifest == manifest_path_for(binary)
    assert manifest.exists()

    data = json.loads(manifest.read_text())
    assert data["tool"] == "fastan"
    assert data["git_sha"] == "abc123"
    assert data["sha256"] == compute_sha256(binary)
    # Signature was captured from the embedded marker (binary is not runnable).
    assert "FASTAN_SATELLOME_SIGNATURE=ad3002:deadbee:0.5" in data["signature"]

    status, message = verify_binary_manifest(binary)
    assert status == "ok", message
    assert "abc123" in message


def test_corrupt_binary_is_a_visible_mismatch(tmp_path):
    """The core regression: a replaced/corrupted binary must NOT verify silently."""
    binary = _make_fake_binary(tmp_path)
    write_binary_manifest(binary, tool="fastan", git_sha="abc123")

    # Tamper with the binary after the manifest was written.
    binary.write_bytes(b"\x7fELF different binary " + MARKER)

    status, message = verify_binary_manifest(binary)
    assert status == "mismatch"
    # Names the concrete discrepancy, not a generic failure.
    assert "sha256" in message.lower()
    assert "mismatch" in message.lower()


def test_corrupt_manifest_is_a_visible_mismatch(tmp_path):
    binary = _make_fake_binary(tmp_path)
    write_binary_manifest(binary, tool="fastan", git_sha="abc123")

    # Corrupt the manifest JSON itself.
    manifest_path_for(binary).write_text("{ this is not valid json ")

    status, message = verify_binary_manifest(binary)
    assert status == "mismatch"
    assert "corrupt" in message.lower() or "unreadable" in message.lower()


def test_missing_manifest_is_unverified_not_silent_ok(tmp_path):
    """No manifest (e.g. user-supplied binary) -> visible 'unverified', not 'ok'."""
    binary = _make_fake_binary(tmp_path)
    # No manifest written.
    status, message = verify_binary_manifest(binary)
    assert status == "unverified"
    assert status != "ok"
    assert "manifest" in message.lower()


def test_manifest_without_sha256_is_mismatch(tmp_path):
    binary = _make_fake_binary(tmp_path)
    manifest_path_for(binary).write_text(json.dumps({"tool": "fastan"}))
    status, message = verify_binary_manifest(binary)
    assert status == "mismatch"
    assert "sha256" in message.lower()


def test_read_build_signature_from_embedded_marker(tmp_path):
    binary = _make_fake_binary(tmp_path)
    sig = read_build_signature(binary)
    assert sig == "FASTAN_SATELLOME_SIGNATURE=ad3002:deadbee:0.5"


def test_resolve_binary_prefers_managed_dir_over_path(tmp_path, monkeypatch):
    """The managed (manifest-backed) binary must win over an unrelated PATH one."""
    managed_dir = tmp_path / "managed"
    managed_dir.mkdir()
    managed = managed_dir / "fastan"
    managed.write_bytes(b"managed")
    managed.chmod(0o755)

    monkeypatch.setattr(base, "get_satellome_bin_dir", lambda: managed_dir)
    # Even if something named fastan exists on PATH, managed wins.
    monkeypatch.setattr(base.shutil, "which", lambda name: "/usr/bin/fastan")

    assert resolve_binary("fastan") == str(managed)


def test_resolve_binary_falls_back_to_path(tmp_path, monkeypatch):
    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()
    monkeypatch.setattr(base, "get_satellome_bin_dir", lambda: empty_dir)
    monkeypatch.setattr(base.shutil, "which", lambda name: "/usr/bin/fastan")
    assert resolve_binary("fastan") == "/usr/bin/fastan"


def test_resolve_binary_returns_none_when_absent(tmp_path, monkeypatch):
    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()
    monkeypatch.setattr(base, "get_satellome_bin_dir", lambda: empty_dir)
    monkeypatch.setattr(base.shutil, "which", lambda name: None)
    assert resolve_binary("fastan") is None
