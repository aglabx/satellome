"""
tanbed (from alntools) installer
"""

import os
import shutil
import subprocess
import tempfile
import logging
from pathlib import Path

from .base import (
    get_satellome_bin_dir,
    check_build_dependencies,
    verify_installation,
    run_make_with_fallback
)

logger = logging.getLogger(__name__)

ALNTOOLS_REPO = "https://github.com/richarddurbin/alntools.git"


def install_tanbed(force: bool = False) -> bool:
    """
    Install tanbed by cloning alntools and compiling from source.

    Args:
        force: Force reinstallation even if binary already exists

    Returns:
        bool: True if installation successful, False otherwise
    """
    logger.info("Starting tanbed installation...")

    # Check if already installed
    bin_dir = get_satellome_bin_dir()
    tanbed_path = bin_dir / 'tanbed'

    if tanbed_path.exists() and not force:
        logger.info(f"tanbed already installed at {tanbed_path}")
        if verify_installation('tanbed'):
            logger.info("tanbed installation verified")
            return True
        else:
            logger.warning("Existing tanbed binary failed verification, reinstalling...")

    # Check build dependencies
    deps_ok, error_msg = check_build_dependencies()
    if not deps_ok:
        logger.error(f"Build dependencies check failed:\n{error_msg}")
        return False

    # Create temporary directory for building
    with tempfile.TemporaryDirectory(prefix='tanbed_build_') as tmp_dir:
        tmp_path = Path(tmp_dir)
        repo_dir = tmp_path / 'alntools'

        try:
            # Clone repository
            logger.info(f"Cloning alntools repository from {ALNTOOLS_REPO}...")
            result = subprocess.run(
                ['git', 'clone', ALNTOOLS_REPO, str(repo_dir)],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=300
            )

            if result.returncode != 0:
                logger.error(f"Failed to clone alntools repository:\n{result.stderr.decode()}")
                return False

            logger.info("Repository cloned successfully")

            # Build alntools (which includes tanbed)
            logger.info("Compiling alntools (including tanbed)...")

            success, error_msg = run_make_with_fallback(repo_dir)
            if not success:
                logger.error(f"Failed to compile alntools:\n{error_msg}")
                return False

            logger.info("alntools compiled successfully")

            # Find tanbed binary (check multiple locations)
            search_dirs = [repo_dir, repo_dir / 'bin']
            binary_source = None

            for search_dir in search_dirs:
                if not search_dir.exists():
                    continue
                candidate = search_dir / 'tanbed'
                if candidate.exists() and os.access(candidate, os.X_OK):
                    binary_source = candidate
                    logger.info(f"Found tanbed binary at {candidate}")
                    break

            if not binary_source:
                logger.error(f"Could not find tanbed binary in {repo_dir} or {repo_dir / 'bin'}")
                logger.error(f"Repository root contents: {list(repo_dir.glob('*'))}")
                if (repo_dir / 'bin').exists():
                    logger.error(f"bin/ directory contents: {list((repo_dir / 'bin').glob('*'))}")
                return False

            # Copy binary to satellome bin directory
            logger.info(f"Installing tanbed to {tanbed_path}...")
            shutil.copy2(binary_source, tanbed_path)
            os.chmod(tanbed_path, 0o755)

            logger.info("tanbed installed successfully!")

            # Verify installation
            if verify_installation('tanbed'):
                logger.info(f"tanbed is ready to use at: {tanbed_path}")
                return True
            else:
                logger.warning("tanbed installed but verification failed")
                return False

        except subprocess.TimeoutExpired:
            logger.error("Installation timed out")
            return False
        except Exception as e:
            logger.error(f"Unexpected error during installation: {e}")
            return False


def uninstall_tanbed() -> bool:
    """
    Uninstall tanbed by removing the binary.

    Returns:
        bool: True if uninstallation successful, False otherwise
    """
    bin_dir = get_satellome_bin_dir()
    tanbed_path = bin_dir / 'tanbed'

    if not tanbed_path.exists():
        logger.info("tanbed is not installed")
        return True

    try:
        tanbed_path.unlink()
        logger.info("tanbed uninstalled successfully")
        return True
    except Exception as e:
        logger.error(f"Failed to uninstall tanbed: {e}")
        return False
