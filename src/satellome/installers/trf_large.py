"""
Modified TRF (for large genomes) installer
"""

import os
import shutil
import subprocess
import tempfile
import logging
from pathlib import Path
from typing import Optional

from .base import (
    get_satellome_bin_dir,
    check_build_dependencies,
    verify_installation
)

logger = logging.getLogger(__name__)

TRF_LARGE_REPO = "https://github.com/aglabx/trf.git"


def install_trf_large(force: bool = False) -> bool:
    """
    Install modified TRF (for large genomes) by cloning and compiling from source.

    Args:
        force: Force reinstallation even if binary already exists

    Returns:
        bool: True if installation successful, False otherwise
    """
    logger.info("Starting modified TRF (for large genomes) installation...")

    # Check platform - known issues on macOS
    from .base import detect_platform
    platform_name, _ = detect_platform()

    if platform_name == 'darwin':
        logger.info("Note: Compiling modified TRF on macOS using build.sh (requires Xcode Command Line Tools).")
        logger.info("If compilation fails, standard TRF will be used as fallback.")

    # Check if already installed
    bin_dir = get_satellome_bin_dir()
    trf_path = bin_dir / 'trf'

    if trf_path.exists() and not force:
        logger.info(f"Modified TRF already installed at {trf_path}")
        if verify_installation('trf'):
            logger.info("Modified TRF installation verified")
            return True
        else:
            logger.warning("Existing modified TRF binary failed verification, reinstalling...")

    # Check build dependencies
    deps_ok, error_msg = check_build_dependencies()
    if not deps_ok:
        logger.error(f"Build dependencies check failed:\n{error_msg}")
        return False

    # Create temporary directory for building
    with tempfile.TemporaryDirectory(prefix='trf_large_build_') as tmp_dir:
        tmp_path = Path(tmp_dir)
        repo_dir = tmp_path / 'trf'

        try:
            # Clone repository
            logger.info(f"Cloning modified TRF repository from {TRF_LARGE_REPO}...")
            result = subprocess.run(
                ['git', 'clone', TRF_LARGE_REPO, str(repo_dir)],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=300
            )

            if result.returncode != 0:
                logger.error(f"Failed to clone modified TRF repository:\n{result.stderr.decode()}")
                return False

            logger.info("Repository cloned successfully")

            # Make build.sh executable
            build_script = repo_dir / 'build.sh'
            if not build_script.exists():
                logger.error("build.sh not found in repository")
                return False

            os.chmod(build_script, 0o755)

            # Build modified TRF using build.sh
            logger.info("Compiling modified TRF using build.sh (this may take a minute)...")
            result = subprocess.run(
                ['./build.sh'],
                cwd=repo_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=300
            )

            if result.returncode != 0:
                logger.error(f"Failed to compile modified TRF:\n{result.stderr.decode()}")
                return False

            logger.info("Modified TRF compiled successfully")

            # Find the binary (should be in repo root)
            binary_source = repo_dir / 'trf'

            if not binary_source.exists():
                logger.error(f"Could not find TRF binary at {binary_source}")
                logger.error(f"Repository contents: {list(repo_dir.glob('*'))}")
                return False

            if not os.access(binary_source, os.X_OK):
                logger.error(f"TRF binary at {binary_source} is not executable")
                return False

            # Copy binary to satellome bin directory
            logger.info(f"Installing modified TRF to {trf_path}...")
            shutil.copy2(binary_source, trf_path)
            os.chmod(trf_path, 0o755)

            logger.info("Modified TRF installed successfully!")

            # Verify installation
            if verify_installation('trf'):
                logger.info(f"Modified TRF is ready to use at: {trf_path}")
                logger.info("Note: This version can handle chromosomes >2GB")
                return True
            else:
                logger.warning("Modified TRF installed but verification failed")
                return False

        except subprocess.TimeoutExpired:
            logger.error("Installation timed out")
            return False
        except Exception as e:
            logger.error(f"Unexpected error during installation: {e}")
            return False


def uninstall_trf_large() -> bool:
    """
    Uninstall modified TRF by removing the binary.

    Returns:
        bool: True if uninstallation successful, False otherwise
    """
    bin_dir = get_satellome_bin_dir()
    trf_path = bin_dir / 'trf'

    if not trf_path.exists():
        logger.info("Modified TRF is not installed")
        return True

    try:
        trf_path.unlink()
        logger.info("Modified TRF uninstalled successfully")
        return True
    except Exception as e:
        logger.error(f"Failed to uninstall modified TRF: {e}")
        return False
