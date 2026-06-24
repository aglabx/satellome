"""GitHub update check for Satellome.

Compares the running version against the newest release tag published on the
project's GitHub repository so users learn when an update is available.

Silent-error policy (see project guidelines):
- An unreachable GitHub / network timeout is an *expected* absence of an
  optional source: it must never crash the pipeline. It is still surfaced as a
  short, visible note (the reason is carried in ``VersionCheckResult.error``),
  never swallowed into nothing.
- The result is cached for ``CHECK_INTERVAL_SECONDS`` so routine runs read a
  local file instead of hitting the API on every invocation, which keeps the
  pipeline fast and stays within GitHub's unauthenticated rate limit.

The check can be disabled entirely with the ``SATELLOME_NO_VERSION_CHECK``
environment variable (any non-empty/truthy value) or the ``--no-version-check``
CLI flag.
"""

import json
import logging
import os
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

logger = logging.getLogger(__name__)

GITHUB_REPO = "aglabx/satellome"
GITHUB_TAGS_URL = f"https://api.github.com/repos/{GITHUB_REPO}/tags"
RELEASES_URL = f"https://github.com/{GITHUB_REPO}/releases"
CHECK_INTERVAL_SECONDS = 24 * 3600  # at most one network check per day on success
FAILURE_CACHE_SECONDS = 3600  # after a failed check, wait an hour before retrying
DEFAULT_TIMEOUT = 3.0
ENV_DISABLE = "SATELLOME_NO_VERSION_CHECK"

# Prefer packaging.version for correct semver comparison; fall back to a small
# numeric-tuple parser so the check still works on a minimal install.
try:  # pragma: no cover - exercised indirectly
    from packaging.version import InvalidVersion, Version

    def _parse_version(value):
        try:
            return Version(str(value).lstrip("vV").strip())
        except (InvalidVersion, TypeError):
            return None

except Exception:  # pragma: no cover - fallback path
    import re

    def _parse_version(value):
        if not value:
            return None
        head = re.split(r"[-+]", str(value).lstrip("vV").strip(), maxsplit=1)[0]
        parts = []
        for chunk in head.split("."):
            match = re.match(r"\d+", chunk)
            if match is None:
                break
            parts.append(int(match.group()))
        return tuple(parts) if parts else None


@dataclass
class VersionCheckResult:
    """Outcome of an update check.

    Attributes:
        current: The running Satellome version (may be "unknown").
        latest: Newest tag found on GitHub, normalized without a leading "v"
            (None if the check could not determine it).
        update_available: True only when ``latest`` is parseable and strictly
            newer than ``current``.
        error: Human-readable reason the check could not complete, or None.
        from_cache: Whether ``latest`` came from the on-disk cache rather than
            a fresh network request.
    """

    current: str
    latest: Optional[str]
    update_available: bool
    error: Optional[str]
    from_cache: bool


def is_disabled() -> bool:
    """Return True when update checks are turned off via the environment."""
    value = os.environ.get(ENV_DISABLE, "").strip().lower()
    return value not in ("", "0", "false", "no", "off")


def _cache_path() -> Path:
    return Path.home() / ".satellome" / "version_check.json"


def _read_cache() -> Optional[dict]:
    path = _cache_path()
    try:
        with open(path, "r") as fh:
            data = json.load(fh)
        if isinstance(data, dict) and "checked_at" in data:
            return data
        return None
    except FileNotFoundError:
        return None
    except (OSError, ValueError) as exc:
        # A corrupt cache is not authoritative; log so it is not invisible,
        # then behave as if no cache exists.
        logger.debug("Ignoring unreadable version-check cache %s: %s", path, exc)
        return None


def _write_cache(latest: Optional[str], ok: bool = True) -> None:
    path = _cache_path()
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as fh:
            json.dump({"checked_at": time.time(), "latest": latest, "ok": ok}, fh)
    except OSError as exc:
        # Caching is best-effort; failing to persist must not break anything,
        # but record why so it is not a fully silent failure.
        logger.debug("Could not write version-check cache %s: %s", path, exc)


def _select_latest(tags: List[str]) -> Optional[str]:
    """Return the highest parseable tag, normalized without a leading 'v'."""
    best_raw = None
    best_key = None
    for tag in tags:
        key = _parse_version(tag)
        if key is None:
            continue
        if best_key is None or key > best_key:
            best_key = key
            best_raw = tag
    if best_raw is None:
        return None
    return str(best_raw).lstrip("vV").strip()


def fetch_latest_version(timeout: float = DEFAULT_TIMEOUT, url: str = GITHUB_TAGS_URL):
    """Fetch the newest release tag from GitHub.

    Returns:
        Tuple[Optional[str], Optional[str]]: ``(latest, error)`` where exactly
        one is non-None. ``latest`` is normalized without a leading "v".
    """
    request = urllib.request.Request(
        url,
        headers={
            "Accept": "application/vnd.github+json",
            "User-Agent": "satellome-version-check",
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            payload = response.read()
        tags = json.loads(payload)
    except urllib.error.HTTPError as exc:
        return None, f"GitHub returned HTTP {exc.code}"
    except urllib.error.URLError as exc:
        return None, f"GitHub unreachable: {exc.reason}"
    except (TimeoutError, OSError) as exc:
        return None, f"network error: {exc}"
    except ValueError as exc:
        return None, f"invalid response from GitHub: {exc}"

    if not isinstance(tags, list):
        return None, "unexpected GitHub response shape (expected a list of tags)"

    names = [t.get("name", "") for t in tags if isinstance(t, dict)]
    latest = _select_latest(names)
    if latest is None:
        return None, "no version-like tags found on GitHub"
    return latest, None


def check_for_updates(
    current_version: str,
    force: bool = False,
    use_cache: bool = True,
    timeout: float = DEFAULT_TIMEOUT,
) -> VersionCheckResult:
    """Determine whether a newer Satellome version is available.

    Args:
        current_version: The running version (e.g. ``satellome.__version__``).
        force: Bypass the daily cache and always query the network.
        use_cache: Read/write the on-disk daily cache.
        timeout: Per-request network timeout in seconds.
    """
    latest = None
    error = None
    from_cache = False

    if use_cache and not force:
        cache = _read_cache()
        if cache is not None:
            age = time.time() - cache.get("checked_at", 0)
            # A previously failed check is retried sooner than a successful one,
            # so a transient outage recovers within the hour while offline runs
            # stay fast and quiet in between.
            ttl = CHECK_INTERVAL_SECONDS if cache.get("ok", True) else FAILURE_CACHE_SECONDS
            if age < ttl:
                latest = cache.get("latest")
                from_cache = True

    if not from_cache:
        latest, error = fetch_latest_version(timeout=timeout)
        if use_cache:
            _write_cache(latest, ok=(error is None))

    update_available = False
    if latest is not None:
        cur_key = _parse_version(current_version)
        latest_key = _parse_version(latest)
        if cur_key is not None and latest_key is not None:
            update_available = latest_key > cur_key

    return VersionCheckResult(
        current=current_version,
        latest=latest,
        update_available=update_available,
        error=error,
        from_cache=from_cache,
    )


def format_update_banner(result: VersionCheckResult) -> Optional[str]:
    """Return a user-facing update notice, or None when there is nothing to say.

    A notice is produced only when a strictly newer version is available; an
    up-to-date or indeterminate result returns None so routine runs stay quiet.
    """
    if not result.update_available or not result.latest:
        return None
    return (
        f"A new Satellome version is available: {result.latest} "
        f"(installed: {result.current}).\n"
        f"  Update with: pip install --upgrade satellome\n"
        f"  Release notes: {RELEASES_URL}"
    )


def notify_if_update_available(current_version: str, log=logger, timeout: float = DEFAULT_TIMEOUT) -> VersionCheckResult:
    """Run a cached, best-effort check and emit a one-line notice.

    Never raises and never blocks beyond ``timeout``. When the check itself
    fails it logs one concise, visible line (not silent) and continues.
    """
    if is_disabled():
        return VersionCheckResult(current_version, None, False, "disabled via environment", False)

    try:
        result = check_for_updates(current_version, timeout=timeout)
    except Exception as exc:  # truly unexpected — stay non-fatal but visible
        log.warning("Update check skipped (unexpected error: %s)", exc)
        return VersionCheckResult(current_version, None, False, str(exc), False)

    banner = format_update_banner(result)
    if banner:
        for line in banner.split("\n"):
            log.info(line)
    elif result.error and not result.from_cache:
        # Optional source was unreachable: surface it as one visible line
        # rather than hiding it, but keep it low-key (it is not a failure of
        # the analysis itself).
        log.info("Update check skipped (%s)", result.error)
    return result
