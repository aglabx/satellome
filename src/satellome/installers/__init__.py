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
from .fastan import (
    install_fastan,
    ensure_fastan_version,
    get_installed_fastan_version,
    FASTAN_MIN_VERSION,
)
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
    'ensure_fastan_version',
    'get_installed_fastan_version',
    'FASTAN_MIN_VERSION',
    'install_tanbed',
    'install_trf_large',
    'install_trf_standard',
]
