"""
FasTAN installer
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

FASTAN_REPO = "https://github.com/ad3002/FASTAN.git"


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
