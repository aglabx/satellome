"""
Installers module for external dependencies (FasTAN, tanbed, etc.)
"""

from .base import (
    detect_platform,
    check_binary_exists,
    check_command_exists,
    get_satellome_bin_dir,
    verify_installation
)
from .fastan import install_fastan
from .tanbed import install_tanbed

__all__ = [
    'detect_platform',
    'check_binary_exists',
    'check_command_exists',
    'get_satellome_bin_dir',
    'verify_installation',
    'install_fastan',
    'install_tanbed',
]
