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
