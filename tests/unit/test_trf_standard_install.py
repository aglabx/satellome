"""Regression test for standard-TRF platform selection.

detect_platform() returns a (platform_name, arch) tuple. trf_standard used to
assign the whole tuple to `platform_name`, so `platform_name not in TRF_URLS`
was always True and the download was never attempted on ANY platform — a real
"always unavailable" bug. This test pins the URL-selection path.
"""

from unittest import mock

from satellome.installers import trf_standard


def test_install_trf_standard_selects_url_for_supported_platform(tmp_path, monkeypatch):
    monkeypatch.setattr(trf_standard, "detect_platform", lambda: ("darwin", "arm64"))
    monkeypatch.setattr(trf_standard, "get_satellome_bin_dir", lambda: tmp_path)
    # Pretend the binary is already-installed check fails so we proceed, and the
    # post-download verification passes.
    monkeypatch.setattr(trf_standard, "verify_installation", lambda name: True)

    fake_response = mock.MagicMock()
    fake_response.__enter__.return_value.read.return_value = b"fake-trf-binary-bytes"
    urlopen = mock.MagicMock(return_value=fake_response)
    monkeypatch.setattr(trf_standard.urllib.request, "urlopen", urlopen)

    result = trf_standard.install_trf_standard(force=True)

    assert result is True, "supported platform must reach the download + succeed"
    # The download was actually attempted with the darwin URL.
    urlopen.assert_called_once()
    called_url = urlopen.call_args.args[0]
    assert called_url == trf_standard.TRF_URLS["darwin"]
    # The binary landed on disk.
    assert (tmp_path / "trf").read_bytes() == b"fake-trf-binary-bytes"


def test_install_trf_standard_unsupported_platform_bails_cleanly(tmp_path, monkeypatch):
    monkeypatch.setattr(trf_standard, "detect_platform", lambda: ("unknown", "unknown"))
    monkeypatch.setattr(trf_standard, "get_satellome_bin_dir", lambda: tmp_path)
    monkeypatch.setattr(trf_standard, "verify_installation", lambda name: False)
    urlopen = mock.MagicMock()
    monkeypatch.setattr(trf_standard.urllib.request, "urlopen", urlopen)

    result = trf_standard.install_trf_standard(force=True)

    assert result is False
    urlopen.assert_not_called()
