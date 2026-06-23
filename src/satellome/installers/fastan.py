"""
FasTAN installer
"""

import os
import json
import shutil
import subprocess
import tempfile
import logging
from pathlib import Path

import re

from .base import (
    get_satellome_bin_dir,
    check_build_dependencies,
    verify_installation,
    run_make_with_fallback,
    write_binary_manifest,
    get_git_commit,
    manifest_path_for,
)

logger = logging.getLogger(__name__)

FASTAN_REPO = "https://github.com/ad3002/FASTAN.git"

# Minimum FasTAN version satellome will run. Older binaries are auto-updated
# (rebuilt from the canonical repo) before the FasTAN step runs. Keep this in
# step with the FasTAN feature/format the pipeline depends on.
FASTAN_MIN_VERSION = "0.8"

_VERSION_RE = re.compile(r"(\d+(?:\.\d+)+)")


def _version_tuple(text):
    """Parse the first dotted-numeric token in *text* into a tuple of ints.

    e.g. ``"FasTAN 0.8 [satellome-fork:...]"`` -> ``(0, 8)``. Returns ``None``
    when no version-like token is present, so callers can tell "older" apart
    from "undeterminable".
    """
    if not text:
        return None
    match = _VERSION_RE.search(str(text))
    if not match:
        return None
    try:
        return tuple(int(part) for part in match.group(1).split("."))
    except ValueError:
        return None


def get_installed_fastan_version(fastan_bin):
    """Return the installed FasTAN version string (e.g. ``"0.8"``), or ``None``.

    The authoritative source is the binary's own ``-V`` output
    (``FasTAN 0.8 [satellome-fork:...]``), reflecting whatever will actually
    run. Falls back to the install-time provenance manifest's ``source_version``
    when the binary cannot be executed (e.g. cross-arch) but a manifest exists.
    """
    fastan_bin = str(fastan_bin)

    # 1) Ask the binary directly.
    try:
        result = subprocess.run(
            [fastan_bin, "-V"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=5,
            check=False,
        )
        out = result.stdout.decode(errors="replace").strip()
        if result.returncode == 0 and out:
            match = re.search(r"FasTAN\s+(\d+(?:\.\d+)+)", out)
            if match:
                return match.group(1)
    except (OSError, subprocess.SubprocessError) as e:
        logger.debug(f"Could not run '{fastan_bin} -V': {e}")

    # 2) Fall back to the provenance manifest recorded at install time.
    manifest = manifest_path_for(fastan_bin)
    try:
        if manifest.exists():
            with open(manifest) as f:
                data = json.load(f)
            source_version = data.get("source_version")
            if source_version:
                return str(source_version)
    except (OSError, ValueError) as e:
        logger.debug(f"Could not read source_version from {manifest}: {e}")

    return None


def ensure_fastan_version(fastan_bin, min_version=FASTAN_MIN_VERSION):
    """Ensure the resolved FasTAN binary is at least *min_version*.

    Reads the installed version; if it is older than required (or cannot be
    determined), forces a reinstall from the canonical repo and re-checks. Every
    decision is logged at a user-visible level — an out-of-date FasTAN is never
    silently tolerated.

    Returns ``(ok, fastan_bin, version)``:
      ok          -- ``False`` if the binary is still below *min_version* after
                     the update attempt (caller MUST abort the FasTAN step).
      fastan_bin  -- path to the binary to actually run; this changes after an
                     update, since the freshly built binary lives in the managed
                     bin dir (which resolve_binary ranks above a PATH copy).
      version     -- the resolved version string, or ``None`` if undeterminable.
    """
    required = _version_tuple(min_version)
    version = get_installed_fastan_version(fastan_bin)
    current = _version_tuple(version)

    if current is not None and (required is None or current >= required):
        logger.info(f"FasTAN version {version} (>= required {min_version})")
        return True, fastan_bin, version

    if current is None:
        logger.warning(
            f"Could not determine FasTAN version for {fastan_bin}; "
            f"updating to guarantee >= {min_version}..."
        )
    else:
        logger.warning(
            f"FasTAN {version} is older than the required {min_version}; updating..."
        )

    if not install_fastan(force=True):
        logger.error(
            f"FasTAN update failed; cannot guarantee version >= {min_version}. "
            f"Reinstall manually with: satellome --install-fastan"
        )
        return False, fastan_bin, version

    # The freshly built binary lives in the managed bin dir.
    new_bin = str(get_satellome_bin_dir() / "fastan")
    if not os.path.exists(new_bin):
        new_bin = fastan_bin
    new_version = get_installed_fastan_version(new_bin)
    new_current = _version_tuple(new_version)

    if new_current is None:
        logger.warning(
            f"FasTAN reinstalled but its version could not be read back "
            f"({new_bin}); proceeding with the freshly built binary."
        )
        return True, new_bin, new_version

    if required is not None and new_current < required:
        logger.error(
            f"FasTAN is still {new_version} after update (required >= {min_version}); "
            f"the canonical repo may not yet ship {min_version}. Aborting FasTAN step."
        )
        return False, new_bin, new_version

    logger.info(f"FasTAN updated to version {new_version}")
    return True, new_bin, new_version


def install_fastan(force: bool = False) -> bool:
    """
    Install FasTAN by cloning and compiling from source.

    Args:
        force: Force reinstallation even if binary already exists

    Returns:
        bool: True if installation successful, False otherwise
    """
    logger.info("Starting FasTAN installation...")

    # Check if already installed
    bin_dir = get_satellome_bin_dir()
    fastan_path = bin_dir / 'fastan'

    if fastan_path.exists() and not force:
        logger.info(f"FasTAN already installed at {fastan_path}")
        if verify_installation('fastan'):
            logger.info("FasTAN installation verified")
            return True
        else:
            logger.warning("Existing FasTAN binary failed verification, reinstalling...")

    # Check build dependencies
    deps_ok, error_msg = check_build_dependencies()
    if not deps_ok:
        logger.error(f"Build dependencies check failed:\n{error_msg}")
        return False

    # Create temporary directory for building
    with tempfile.TemporaryDirectory(prefix='fastan_build_') as tmp_dir:
        tmp_path = Path(tmp_dir)
        repo_dir = tmp_path / 'FASTAN'

        try:
            # Clone repository
            logger.info(f"Cloning FasTAN repository from {FASTAN_REPO}...")
            result = subprocess.run(
                ['git', 'clone', FASTAN_REPO, str(repo_dir)],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=300
            )

            if result.returncode != 0:
                logger.error(f"Failed to clone FasTAN repository:\n{result.stderr.decode()}")
                return False

            logger.info("Repository cloned successfully")

            # Capture provenance from the checkout before the temp dir is wiped.
            git_sha = get_git_commit(repo_dir)
            source_version = None
            fastan_source = repo_dir / 'FasTAN.c'
            if fastan_source.exists():
                try:
                    text = fastan_source.read_text(errors='replace')
                    match = re.search(r'#define\s+VERSION\s+"([^"]+)"', text)
                    if match:
                        source_version = match.group(1)
                except OSError as e:
                    logger.debug(f"Could not read source version from {fastan_source}: {e}")

            # Patch Makefile to add pthread linking flag
            makefile_path = repo_dir / 'Makefile'
            if makefile_path.exists():
                logger.info("Patching Makefile to add pthread support...")
                try:
                    with open(makefile_path, 'r') as f:
                        makefile_content = f.read()

                    # Add -lpthread to LDFLAGS if not already present
                    if '-lpthread' not in makefile_content:
                        # Try to find and modify LDFLAGS line
                        if 'LDFLAGS' in makefile_content:
                            makefile_content = makefile_content.replace(
                                'LDFLAGS =',
                                'LDFLAGS = -lpthread'
                            )
                        else:
                            # Add LDFLAGS at the beginning
                            makefile_content = 'LDFLAGS = -lpthread\n\n' + makefile_content

                        # Also ensure the linker command uses LDFLAGS
                        # Common pattern: gcc ... -o target
                        # Should be: gcc ... $(LDFLAGS) -o target
                        lines = makefile_content.split('\n')
                        modified_lines = []
                        for line in lines:
                            # If it's a gcc/cc link command without $(LDFLAGS)
                            if (('gcc' in line or 'cc' in line or '$(CC)' in line) and
                                '-o' in line and
                                '$(LDFLAGS)' not in line and
                                not line.strip().startswith('#')):
                                # Insert $(LDFLAGS) before -o
                                line = line.replace(' -o ', ' $(LDFLAGS) -o ')
                            modified_lines.append(line)
                        makefile_content = '\n'.join(modified_lines)

                        with open(makefile_path, 'w') as f:
                            f.write(makefile_content)
                        logger.info("Makefile patched successfully")
                except Exception as e:
                    logger.warning(f"Could not patch Makefile: {e}. Attempting build anyway...")

            # Build FasTAN
            logger.info("Compiling FasTAN...")

            success, error_msg = run_make_with_fallback(repo_dir)
            if not success:
                logger.error(f"Failed to compile FasTAN:\n{error_msg}")
                return False

            logger.info("FasTAN compiled successfully")

            # Find the binary (check common names in multiple locations)
            possible_names = ['FasTAN', 'fastan', 'FASTAN']
            search_dirs = [repo_dir, repo_dir / 'bin']
            binary_source = None

            for search_dir in search_dirs:
                if not search_dir.exists():
                    continue
                for name in possible_names:
                    candidate = search_dir / name
                    if candidate.exists() and os.access(candidate, os.X_OK):
                        binary_source = candidate
                        logger.info(f"Found FasTAN binary at {candidate}")
                        break
                if binary_source:
                    break

            if not binary_source:
                logger.error(f"Could not find FasTAN binary in {repo_dir} or {repo_dir / 'bin'}")
                logger.error(f"Repository root contents: {list(repo_dir.glob('*'))}")
                if (repo_dir / 'bin').exists():
                    logger.error(f"bin/ directory contents: {list((repo_dir / 'bin').glob('*'))}")
                return False

            # Copy binary to satellome bin directory
            logger.info(f"Installing FasTAN to {fastan_path}...")
            shutil.copy2(binary_source, fastan_path)
            os.chmod(fastan_path, 0o755)

            logger.info("FasTAN installed successfully!")

            # Record provenance/integrity manifest next to the binary.
            write_binary_manifest(
                fastan_path,
                tool='fastan',
                repo=FASTAN_REPO,
                git_sha=git_sha,
                source_version=source_version,
            )

            # Verify installation
            if verify_installation('fastan'):
                logger.info(f"FasTAN is ready to use at: {fastan_path}")
                return True
            else:
                logger.warning("FasTAN installed but verification failed")
                return False

        except subprocess.TimeoutExpired:
            logger.error("Installation timed out")
            return False
        except Exception as e:
            logger.error(f"Unexpected error during installation: {e}")
            return False


def uninstall_fastan() -> bool:
    """
    Uninstall FasTAN by removing the binary.

    Returns:
        bool: True if uninstallation successful, False otherwise
    """
    bin_dir = get_satellome_bin_dir()
    fastan_path = bin_dir / 'fastan'

    if not fastan_path.exists():
        logger.info("FasTAN is not installed")
        return True

    try:
        fastan_path.unlink()
        logger.info("FasTAN uninstalled successfully")
        return True
    except Exception as e:
        logger.error(f"Failed to uninstall FasTAN: {e}")
        return False
