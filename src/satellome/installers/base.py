"""
Base utilities for installers
"""

import os
import sys
import shutil
import platform
import logging
from pathlib import Path
from typing import Optional, Tuple

logger = logging.getLogger(__name__)


def detect_platform() -> Tuple[str, str]:
    """
    Detect the current platform and architecture.

    Returns:
        Tuple[str, str]: (platform_name, architecture)
        platform_name: 'linux', 'darwin' (macOS), or 'unknown'
        architecture: 'x86_64', 'arm64', or 'unknown'
    """
    system = platform.system().lower()
    machine = platform.machine().lower()

    # Normalize platform names
    if system == 'linux':
        platform_name = 'linux'
    elif system == 'darwin':
        platform_name = 'darwin'
    else:
        platform_name = 'unknown'

    # Normalize architecture
    if machine in ['x86_64', 'amd64']:
        arch = 'x86_64'
    elif machine in ['arm64', 'aarch64']:
        arch = 'arm64'
    else:
        arch = 'unknown'

    return platform_name, arch


def get_satellome_bin_dir() -> Path:
    """
    Get or create the Satellome binary directory.

    Returns:
        Path: Path to ~/.satellome/bin/
    """
    bin_dir = Path.home() / '.satellome' / 'bin'
    bin_dir.mkdir(parents=True, exist_ok=True)
    return bin_dir


def check_command_exists(command: str) -> bool:
    """
    Check if a command is available in the system PATH.

    Args:
        command: Command name to check

    Returns:
        bool: True if command exists, False otherwise
    """
    return shutil.which(command) is not None


def check_binary_exists(binary_name: str, check_system_path: bool = True) -> Optional[str]:
    """
    Check if a binary exists in Satellome bin directory or system PATH.

    Args:
        binary_name: Name of the binary to check
        check_system_path: Whether to check system PATH in addition to ~/.satellome/bin/

    Returns:
        Optional[str]: Full path to binary if found, None otherwise
    """
    # Check ~/.satellome/bin/ first (higher priority)
    satellome_bin = get_satellome_bin_dir() / binary_name
    if satellome_bin.exists() and os.access(satellome_bin, os.X_OK):
        return str(satellome_bin)

    # Check system PATH if requested
    if check_system_path:
        system_path = shutil.which(binary_name)
        if system_path:
            return system_path

    return None


def verify_installation(binary_name: str, test_command: Optional[str] = None) -> bool:
    """
    Verify that a binary is properly installed and executable.

    Args:
        binary_name: Name of the binary to verify
        test_command: Optional command to test (e.g., "fastan --help")

    Returns:
        bool: True if verification successful, False otherwise
    """
    binary_path = check_binary_exists(binary_name)

    if not binary_path:
        logger.error(f"{binary_name} not found in PATH or ~/.satellome/bin/")
        return False

    logger.info(f"Found {binary_name} at: {binary_path}")

    # If test command provided, try to run it
    if test_command:
        import subprocess
        try:
            result = subprocess.run(
                test_command.split(),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=5
            )
            # Some programs return non-zero for --help, so we check if it runs without crash
            logger.info(f"{binary_name} test command executed successfully")
            return True
        except (subprocess.TimeoutExpired, subprocess.SubprocessError, FileNotFoundError) as e:
            logger.error(f"Failed to run test command for {binary_name}: {e}")
            return False

    return True


def check_build_dependencies() -> Tuple[bool, str]:
    """
    Check if required build tools are available.

    Returns:
        Tuple[bool, str]: (success, error_message)
    """
    missing = []

    # Check for git
    if not check_command_exists('git'):
        missing.append('git')

    # Check for make
    if not check_command_exists('make'):
        missing.append('make')

    # Check for C compiler (gcc or clang)
    has_compiler = check_command_exists('gcc') or check_command_exists('clang') or check_command_exists('cc')
    if not has_compiler:
        missing.append('gcc or clang')

    if missing:
        error_msg = f"Missing required build tools: {', '.join(missing)}\n"
        system_name, _ = detect_platform()

        if system_name == 'darwin':
            error_msg += "On macOS, install Xcode Command Line Tools: xcode-select --install"
        elif system_name == 'linux':
            error_msg += "On Ubuntu/Debian: sudo apt-get install build-essential git\n"
            error_msg += "On CentOS/RHEL: sudo yum groupinstall 'Development Tools' && sudo yum install git"

        return False, error_msg

    return True, ""
