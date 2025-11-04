#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Installer for standard TRF (Tandem Repeat Finder)
Downloads pre-compiled binary from official GitHub releases
"""

import logging
import urllib.request
import urllib.error
import stat
from pathlib import Path

from .base import (
    detect_platform,
    get_satellome_bin_dir,
    verify_installation,
)

logger = logging.getLogger(__name__)

# TRF download URLs
TRF_URLS = {
    'linux': 'https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.linux64',
    'darwin': 'https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.macos',
}


def install_trf_standard(force: bool = False) -> bool:
    """
    Install standard TRF by downloading pre-compiled binary.

    Args:
        force: If True, reinstall even if already exists

    Returns:
        True if installation successful, False otherwise
    """
    logger.info("Starting standard TRF installation...")

    platform_name = detect_platform()
    bin_dir = get_satellome_bin_dir()
    trf_binary = bin_dir / 'trf'

    # Check if already installed
    if not force and verify_installation('trf'):
        logger.info(f"Standard TRF is already installed at: {trf_binary}")
        return True

    # Get download URL for platform
    if platform_name not in TRF_URLS:
        logger.error(f"Standard TRF binary not available for {platform_name}")
        logger.info("Please download manually from: https://tandem.bu.edu/trf/trf.html")
        return False

    url = TRF_URLS[platform_name]
    logger.info(f"Downloading TRF from {url}...")

    try:
        # Download binary
        with urllib.request.urlopen(url, timeout=60) as response:
            binary_data = response.read()

        # Write to file
        trf_binary.write_bytes(binary_data)

        # Make executable
        trf_binary.chmod(trf_binary.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)

        logger.info(f"TRF binary downloaded to: {trf_binary}")

        # Verify installation
        if verify_installation('trf'):
            logger.info("Standard TRF installed successfully!")
            logger.info(f"TRF is ready to use at: {trf_binary}")
            return True
        else:
            logger.error("TRF binary downloaded but verification failed")
            return False

    except urllib.error.URLError as e:
        logger.error(f"Failed to download TRF: {e}")
        return False
    except Exception as e:
        logger.error(f"Unexpected error during TRF installation: {e}")
        return False
