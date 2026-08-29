"""Installer for the Rust helper tools whose sources live in this repository.

``rust/`` holds five small crates — sat-family, telomere-check, find-gaps,
bed-extract, genome-size — that the pipeline shells out to. Until now nothing
could install them: they are not in the wheel (they are per-platform native
binaries, and a ``py3-none-any`` wheel must not carry a macOS executable to a
Linux host), they were not in the sdist, and ``--install-all`` did not know
about them. A pip user therefore never had them, and two of them are the
difference between a run that clusters satellite families and one that quietly
does not.

Sources are taken from the local checkout when there is one, and cloned from
the project repository otherwise, so this works from a wheel install too.
"""

import logging
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

from satellome.installers.base import get_satellome_bin_dir

logger = logging.getLogger(__name__)

SATELLOME_REPO = "https://github.com/aglabx/satellome.git"

# Prebuilt binaries are attached to the GitHub Release for each version by
# .github/workflows/rust-tools.yml, together with a SHA256SUMS file per
# platform. Downloading one is the fast path; building with cargo is the
# fallback for a platform with no asset or an unreachable network.
RELEASE_ASSET_URL = (
    "https://github.com/aglabx/satellome/releases/download/v{version}/{asset}"
)
DOWNLOAD_TIMEOUT = 120
ENV_NO_DOWNLOAD = "SATELLOME_NO_PREBUILT"

# Crate directory name -> produced binary name (identical here, but keep the
# mapping explicit so a crate rename does not silently install nothing).
RUST_TOOLS = (
    "sat-family",
    "telomere-check",
    "find-gaps",
    "bed-extract",
    "genome-size",
)

BUILD_TIMEOUT = 1800
CLONE_TIMEOUT = 300


def find_rust_sources() -> Optional[Path]:
    """Locate the ``rust/`` source tree of an editable/source installation."""
    import satellome

    package_dir = Path(satellome.__file__).resolve().parent
    candidates = [
        package_dir / "rust",              # shipped as package data
        package_dir.parent.parent / "rust",  # src/satellome/.. -> repo root
        package_dir.parent / "rust",
    ]
    for candidate in candidates:
        if (candidate / RUST_TOOLS[0] / "Cargo.toml").is_file():
            return candidate
    return None


def source_tree_version() -> Optional[str]:
    """Version declared by the checkout we are running from, if any.

    ``satellome.__version__`` comes from installed package metadata, which in a
    source checkout can belong to an entirely different (older) install that
    happens to be on the path. Asking a release for assets under that version
    finds nothing and falls back to building, which works but for a reason the
    user has no way to guess. The checkout's own pyproject.toml is the truth
    when we are running from one.

    Returns None when this is a normal installation.
    """
    import satellome

    package_dir = Path(satellome.__file__).resolve().parent
    # src-layout checkout: <root>/src/satellome -> <root>; flat: <root>/satellome
    for root in (package_dir.parent.parent, package_dir.parent):
        pyproject = root / "pyproject.toml"
        if not pyproject.is_file():
            continue
        try:
            text = pyproject.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        # Only trust satellome's own pyproject, never an unrelated parent
        # project that happens to sit above the package directory.
        if 'name = "satellome"' not in text:
            continue
        for line in text.splitlines():
            stripped = line.strip()
            if stripped.startswith("version") and "=" in stripped:
                value = stripped.split("=", 1)[1].strip().strip('"').strip("'")
                if value:
                    return value
    return None


def asset_version(log=None) -> str:
    """The release whose assets match the code that is actually running."""
    log = log or logger
    from satellome import __version__ as installed

    from_source = source_tree_version()
    if from_source and from_source != installed:
        log.info(
            f"Running from a source checkout at version {from_source}, while the "
            f"installed metadata says {installed}. Fetching assets for "
            f"{from_source}."
        )
        return from_source
    return installed


def platform_asset_suffix() -> Optional[str]:
    """Asset suffix for this machine, or None if we publish nothing for it."""
    from satellome.installers.base import detect_platform

    platform_name, arch = detect_platform()
    mapping = {
        ("linux", "x86_64"): "linux-x86_64",
        ("darwin", "arm64"): "macos-arm64",
        ("darwin", "x86_64"): "macos-x86_64",
    }
    return mapping.get((platform_name, arch))


def _fetch(url: str, timeout: int = DOWNLOAD_TIMEOUT) -> Optional[bytes]:
    """GET a URL, returning None (with a reason logged) on any failure."""
    import urllib.error
    import urllib.request

    try:
        with urllib.request.urlopen(url, timeout=timeout) as response:
            return response.read()
    except urllib.error.HTTPError as error:
        # 404 is the expected "no asset for this version/platform", and is the
        # signal to fall back to building rather than an error to shout about.
        logger.debug(f"{url} -> HTTP {error.code}")
        return None
    except (urllib.error.URLError, OSError, ValueError) as error:
        logger.warning(f"Could not download {url}: {error}")
        return None


def _expected_checksums(version: str, suffix: str) -> dict:
    """filename -> sha256, from the SHA256SUMS file published with the release."""
    url = RELEASE_ASSET_URL.format(version=version, asset=f"SHA256SUMS-{suffix}.txt")
    data = _fetch(url, timeout=30)
    if data is None:
        return {}

    checksums = {}
    for line in data.decode("utf-8", "replace").splitlines():
        parts = line.split()
        if len(parts) >= 2:
            checksums[parts[-1].lstrip("*")] = parts[0]
    return checksums


def download_prebuilt(name: str, version: str, bin_dir: Path) -> bool:
    """Install one tool from the release assets, verifying its checksum.

    Returns False for "not available" as well as for a failed verification;
    the caller falls back to building. A binary whose checksum does not match
    is never installed - a wrong binary is worse than a missing one, because
    the pipeline would run it.
    """
    import hashlib

    suffix = platform_asset_suffix()
    if suffix is None:
        return False

    asset = f"{name}-{suffix}"
    expected = _expected_checksums(version, suffix).get(asset)
    if expected is None:
        logger.debug(f"No published checksum for {asset} in v{version}")
        return False

    data = _fetch(RELEASE_ASSET_URL.format(version=version, asset=asset))
    if data is None:
        return False

    actual = hashlib.sha256(data).hexdigest()
    if actual != expected:
        logger.error(
            f"Checksum mismatch for {asset}: expected {expected}, got {actual}. "
            "Refusing to install it."
        )
        return False

    target = bin_dir / name
    partial = bin_dir / f"{name}.partial"
    try:
        partial.write_bytes(data)
        partial.chmod(partial.stat().st_mode | 0o111)
        os.replace(str(partial), str(target))
    except OSError as error:
        logger.error(f"Could not install {name} into {bin_dir}: {error}")
        return False

    logger.info(f"✓ {name} downloaded ({asset}) to {target}")
    return True


def check_cargo() -> Tuple[bool, str]:
    """Is a Rust toolchain available to build with?"""
    cargo = shutil.which("cargo")
    if not cargo:
        return False, (
            "cargo (the Rust toolchain) was not found. These tools are compiled "
            "from source. Install Rust with:\n"
            "  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh\n"
            "then re-run: satellome --install-rust-tools"
        )
    return True, cargo


def _clone_sources(destination: Path) -> Optional[Path]:
    """Clone the project repository to obtain ``rust/`` for a wheel install."""
    if not shutil.which("git"):
        logger.error("git is required to fetch the Rust sources but was not found")
        return None

    repo_dir = destination / "satellome"
    logger.info(f"Fetching Rust sources from {SATELLOME_REPO}...")
    result = subprocess.run(
        ["git", "clone", "--depth", "1", SATELLOME_REPO, str(repo_dir)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=CLONE_TIMEOUT,
    )
    if result.returncode != 0:
        logger.error(
            "Could not clone the Rust sources: "
            f"{result.stderr.decode('utf-8', 'replace').strip()}"
        )
        return None

    rust_dir = repo_dir / "rust"
    if not (rust_dir / RUST_TOOLS[0] / "Cargo.toml").is_file():
        logger.error(f"Cloned repository has no usable rust/ directory at {rust_dir}")
        return None
    return rust_dir


def _build_one(crate_dir: Path, name: str, bin_dir: Path) -> bool:
    """cargo build --release, then publish the binary into the managed bin dir."""
    logger.info(f"Building {name}...")
    result = subprocess.run(
        ["cargo", "build", "--release"],
        cwd=str(crate_dir),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=BUILD_TIMEOUT,
    )
    if result.returncode != 0:
        # A compiler error must be shown: "build failed" alone leaves the user
        # with nothing to act on.
        logger.error(
            f"Building {name} failed:\n"
            + result.stderr.decode("utf-8", "replace").strip()[-2000:]
        )
        return False

    built = crate_dir / "target" / "release" / name
    if not built.is_file():
        logger.error(f"{name} built without error but no binary at {built}")
        return False

    target = bin_dir / name
    try:
        shutil.copy2(str(built), str(target))
        target.chmod(target.stat().st_mode | 0o111)
    except OSError as error:
        logger.error(f"Could not install {name} into {bin_dir}: {error}")
        return False

    logger.info(f"✓ {name} installed at {target}")
    return True


def install_rust_tools(force: bool = False, tools: Optional[Sequence[str]] = None) -> bool:
    """Build and install the bundled Rust helper tools.

    Returns True only when every requested tool ends up installed; a partial
    result is reported as a failure with the names that are still missing, so
    the caller cannot mistake "three of five" for success.
    """
    wanted: List[str] = list(tools) if tools else list(RUST_TOOLS)
    bin_dir = get_satellome_bin_dir()

    if not force:
        remaining = [
            name for name in wanted
            if not (bin_dir / name).is_file() or not os.access(bin_dir / name, os.X_OK)
        ]
        already = [name for name in wanted if name not in remaining]
        if already:
            logger.info(f"Already installed: {', '.join(already)}")
        if not remaining:
            return True
        wanted = remaining

    # Fast path: prebuilt binaries published with this version's release.
    if not os.environ.get(ENV_NO_DOWNLOAD):
        version = asset_version()

        suffix = platform_asset_suffix()
        if suffix is None:
            logger.info(
                "No prebuilt binaries are published for this platform; "
                "building from source."
            )
        else:
            logger.info(f"Fetching prebuilt binaries for {suffix} (v{version})...")
            still_needed = [
                name for name in wanted
                if not download_prebuilt(name, version, bin_dir)
            ]
            if not still_needed:
                logger.info(f"Installed: {', '.join(wanted)}")
                return True
            if len(still_needed) < len(wanted):
                logger.info(
                    f"Downloaded {len(wanted) - len(still_needed)} of {len(wanted)}; "
                    f"building the rest: {', '.join(still_needed)}"
                )
            else:
                logger.info(
                    "No prebuilt binaries available for this version; "
                    "building from source."
                )
            wanted = still_needed

    cargo_ok, message = check_cargo()
    if not cargo_ok:
        logger.error(message)
        return False

    sources = find_rust_sources()
    with tempfile.TemporaryDirectory(prefix="satellome_rust_") as tmp_dir:
        if sources is None:
            logger.info("No local rust/ sources (wheel install); fetching them.")
            sources = _clone_sources(Path(tmp_dir))
            if sources is None:
                return False
        else:
            logger.info(f"Building from local sources: {sources}")

        installed, failed = [], []
        for name in wanted:
            crate_dir = sources / name
            if not (crate_dir / "Cargo.toml").is_file():
                logger.error(f"No crate for {name} at {crate_dir}")
                failed.append(name)
                continue
            (installed if _build_one(crate_dir, name, bin_dir) else failed).append(name)

    if installed:
        logger.info(f"Installed: {', '.join(installed)}")
    if failed:
        logger.error(f"Still missing after installation: {', '.join(failed)}")
        return False
    return True
