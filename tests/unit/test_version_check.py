#!/usr/bin/env python
"""Unit tests for the GitHub update-check helper.

These tests verify the actual decision (update available / not / indeterminate)
and the failure surface — not just that functions run. In particular the
network-failure path must carry a visible error and must never raise.
"""

import json
import urllib.error

import pytest

from satellome.core_functions.tools import version_check as vc


class TestSelectLatest:
    def test_picks_highest_semver_ignoring_v_prefix(self):
        tags = ["v1.4.0", "v1.6.1", "v1.5.2", "v1.1.0"]
        assert vc._select_latest(tags) == "1.6.1"

    def test_ignores_non_version_tags(self):
        tags = ["nightly", "v1.2.0", "latest", "v1.10.0"]
        # 1.10.0 must beat 1.2.0 numerically (not lexically)
        assert vc._select_latest(tags) == "1.10.0"

    def test_returns_none_when_no_version_tags(self):
        assert vc._select_latest(["nightly", "latest"]) is None


class TestCheckForUpdates:
    def test_update_available_when_remote_is_newer(self, monkeypatch):
        monkeypatch.setattr(vc, "fetch_latest_version", lambda **kw: ("1.6.5", None))
        result = vc.check_for_updates("1.6.2", use_cache=False)
        assert result.update_available is True
        assert result.latest == "1.6.5"
        assert result.error is None

    def test_no_update_when_versions_equal(self, monkeypatch):
        monkeypatch.setattr(vc, "fetch_latest_version", lambda **kw: ("1.6.2", None))
        result = vc.check_for_updates("1.6.2", use_cache=False)
        assert result.update_available is False
        assert result.error is None

    def test_no_update_when_local_is_newer(self, monkeypatch):
        monkeypatch.setattr(vc, "fetch_latest_version", lambda **kw: ("1.6.1", None))
        result = vc.check_for_updates("1.6.2", use_cache=False)
        assert result.update_available is False

    def test_unknown_local_version_cannot_claim_update(self, monkeypatch):
        monkeypatch.setattr(vc, "fetch_latest_version", lambda **kw: ("1.6.2", None))
        result = vc.check_for_updates("unknown", use_cache=False)
        assert result.update_available is False
        assert result.latest == "1.6.2"

    def test_network_failure_is_surfaced_not_swallowed(self, monkeypatch):
        # Reproduces the silent-failure path: the check fails but must report a
        # visible reason and must not raise or claim an update.
        monkeypatch.setattr(
            vc, "fetch_latest_version", lambda **kw: (None, "GitHub unreachable: timed out")
        )
        result = vc.check_for_updates("1.6.2", use_cache=False)
        assert result.update_available is False
        assert result.latest is None
        assert result.error is not None
        assert "unreachable" in result.error


class TestCache:
    def test_fresh_cache_avoids_network(self, monkeypatch, tmp_path):
        cache_file = tmp_path / "version_check.json"
        monkeypatch.setattr(vc, "_cache_path", lambda: cache_file)

        def _boom(**kw):
            raise AssertionError("network must not be hit when cache is fresh")

        # Seed a fresh cache, then ensure no network call happens.
        vc._write_cache("1.6.9")
        monkeypatch.setattr(vc, "fetch_latest_version", _boom)
        result = vc.check_for_updates("1.6.2", use_cache=True)
        assert result.from_cache is True
        assert result.latest == "1.6.9"
        assert result.update_available is True

    def test_stale_cache_triggers_network(self, monkeypatch, tmp_path):
        cache_file = tmp_path / "version_check.json"
        cache_file.write_text(json.dumps({"checked_at": 0, "latest": "1.0.0"}))
        monkeypatch.setattr(vc, "_cache_path", lambda: cache_file)
        monkeypatch.setattr(vc, "fetch_latest_version", lambda **kw: ("1.6.3", None))
        result = vc.check_for_updates("1.6.2", use_cache=True)
        assert result.from_cache is False
        assert result.latest == "1.6.3"

    def test_recent_failure_is_cached_to_avoid_hammering(self, monkeypatch, tmp_path):
        cache_file = tmp_path / "version_check.json"
        # A failed check recorded a moment ago (ok=False) must be reused rather
        # than re-hitting the network on every offline run.
        cache_file.write_text(json.dumps({"checked_at": vc.time.time(), "latest": None, "ok": False}))
        monkeypatch.setattr(vc, "_cache_path", lambda: cache_file)

        def _boom(**kw):
            raise AssertionError("network must not be hit during the failure-cache window")

        monkeypatch.setattr(vc, "fetch_latest_version", _boom)
        result = vc.check_for_updates("1.6.2", use_cache=True)
        assert result.from_cache is True
        assert result.latest is None
        assert result.update_available is False

    def test_old_failure_is_retried(self, monkeypatch, tmp_path):
        cache_file = tmp_path / "version_check.json"
        # Failure older than FAILURE_CACHE_SECONDS but younger than the success
        # interval must still trigger a fresh attempt.
        stale = vc.time.time() - (vc.FAILURE_CACHE_SECONDS + 60)
        cache_file.write_text(json.dumps({"checked_at": stale, "latest": None, "ok": False}))
        monkeypatch.setattr(vc, "_cache_path", lambda: cache_file)
        monkeypatch.setattr(vc, "fetch_latest_version", lambda **kw: ("1.6.6", None))
        result = vc.check_for_updates("1.6.2", use_cache=True)
        assert result.from_cache is False
        assert result.latest == "1.6.6"

    def test_failed_fetch_is_persisted_as_not_ok(self, monkeypatch, tmp_path):
        cache_file = tmp_path / "version_check.json"
        monkeypatch.setattr(vc, "_cache_path", lambda: cache_file)
        monkeypatch.setattr(vc, "fetch_latest_version", lambda **kw: (None, "offline"))
        vc.check_for_updates("1.6.2", use_cache=True)
        saved = json.loads(cache_file.read_text())
        assert saved["ok"] is False
        assert saved["latest"] is None

    def test_corrupt_cache_is_ignored_not_fatal(self, monkeypatch, tmp_path):
        cache_file = tmp_path / "version_check.json"
        cache_file.write_text("{ this is not json ")
        monkeypatch.setattr(vc, "_cache_path", lambda: cache_file)
        monkeypatch.setattr(vc, "fetch_latest_version", lambda **kw: ("1.6.4", None))
        result = vc.check_for_updates("1.6.2", use_cache=True)
        # Falls through to the network instead of crashing on bad JSON.
        assert result.latest == "1.6.4"
        assert result.from_cache is False


class TestFetchLatestVersion:
    def _fake_urlopen(self, payload):
        class _Resp:
            def __enter__(self_inner):
                return self_inner

            def __exit__(self_inner, *a):
                return False

            def read(self_inner):
                return json.dumps(payload).encode()

        return lambda req, timeout=None: _Resp()

    def test_parses_tags_and_returns_highest(self, monkeypatch):
        tags = [{"name": "v1.4.0"}, {"name": "v1.6.1"}, {"name": "v1.5.0"}]
        monkeypatch.setattr(vc.urllib.request, "urlopen", self._fake_urlopen(tags))
        latest, error = vc.fetch_latest_version(timeout=1.0)
        assert latest == "1.6.1"
        assert error is None

    def test_http_error_returns_reason(self, monkeypatch):
        def _raise(req, timeout=None):
            raise urllib.error.HTTPError(vc.GITHUB_TAGS_URL, 403, "rate limited", {}, None)

        monkeypatch.setattr(vc.urllib.request, "urlopen", _raise)
        latest, error = vc.fetch_latest_version(timeout=1.0)
        assert latest is None
        assert "403" in error

    def test_url_error_returns_reason(self, monkeypatch):
        def _raise(req, timeout=None):
            raise urllib.error.URLError("name resolution failed")

        monkeypatch.setattr(vc.urllib.request, "urlopen", _raise)
        latest, error = vc.fetch_latest_version(timeout=1.0)
        assert latest is None
        assert "unreachable" in error

    def test_unexpected_shape_is_reported(self, monkeypatch):
        monkeypatch.setattr(vc.urllib.request, "urlopen", self._fake_urlopen({"not": "a list"}))
        latest, error = vc.fetch_latest_version(timeout=1.0)
        assert latest is None
        assert "shape" in error


class TestBanner:
    def test_banner_only_when_update_available(self):
        up = vc.VersionCheckResult("1.6.2", "1.7.0", True, None, False)
        banner = vc.format_update_banner(up)
        assert banner is not None
        assert "1.7.0" in banner
        assert "1.6.2" in banner
        assert "pip install --upgrade satellome" in banner
        assert vc.RELEASES_URL in banner

    def test_no_banner_when_up_to_date(self):
        same = vc.VersionCheckResult("1.6.2", "1.6.2", False, None, False)
        assert vc.format_update_banner(same) is None

    def test_no_banner_when_indeterminate(self):
        err = vc.VersionCheckResult("1.6.2", None, False, "offline", False)
        assert vc.format_update_banner(err) is None


class TestNotify:
    def test_disabled_via_env(self, monkeypatch):
        monkeypatch.setenv(vc.ENV_DISABLE, "1")
        result = vc.notify_if_update_available("1.6.2")
        assert result.update_available is False
        assert result.error == "disabled via environment"

    def test_notify_never_raises_on_internal_error(self, monkeypatch):
        monkeypatch.delenv(vc.ENV_DISABLE, raising=False)

        def _boom(*a, **kw):
            raise RuntimeError("kaboom")

        monkeypatch.setattr(vc, "check_for_updates", _boom)
        # Must swallow the *unexpected* error into a non-fatal, visible result.
        result = vc.notify_if_update_available("1.6.2")
        assert result.update_available is False
        assert "kaboom" in result.error


class TestIsDisabled:
    @pytest.mark.parametrize("value,expected", [
        ("1", True), ("true", True), ("yes", True), ("on", True),
        ("0", False), ("false", False), ("no", False), ("", False),
    ])
    def test_env_values(self, monkeypatch, value, expected):
        monkeypatch.setenv(vc.ENV_DISABLE, value)
        assert vc.is_disabled() is expected
