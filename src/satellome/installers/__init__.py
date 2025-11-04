"""
Installers module for external dependencies (FasTAN, tanbed, TRF, etc.)
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
from .trf_large import install_trf_large
from .trf_standard import install_trf_standard

__all__ = [
    'detect_platform',
    'check_binary_exists',
    'check_command_exists',
    'get_satellome_bin_dir',
    'verify_installation',
    'install_fastan',
    'install_tanbed',
    'install_trf_large',
    'install_trf_standard',
]
