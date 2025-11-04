"""
Custom install command for setuptools to auto-install external tools
"""

import os
import sys
import logging
from setuptools.command.install import install

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(message)s'
)
logger = logging.getLogger('satellome.install')


class PostInstallCommand(install):
    """
    Post-installation command to automatically install external tools.

    Installs FasTAN, tanbed, and modified TRF during pip install.
    Can be skipped by setting SATELLOME_SKIP_AUTO_INSTALL environment variable.
    """

    def run(self):
        # Run standard install first
        install.run(self)

        # Check if auto-install is disabled
        if os.environ.get('SATELLOME_SKIP_AUTO_INSTALL'):
            logger.info("\n" + "="*60)
            logger.info("Satellome installed successfully!")
            logger.info("Skipping automatic tool installation (SATELLOME_SKIP_AUTO_INSTALL is set)")
            logger.info("You can install tools later with: satellome --install-all")
            logger.info("="*60 + "\n")
            return

        logger.info("\n" + "="*60)
        logger.info("Satellome post-installation: Installing external tools...")
        logger.info("="*60)
        logger.info("")
        logger.info("This will attempt to compile and install:")
        logger.info("  • FasTAN (alternative tandem repeat finder)")
        logger.info("  • tanbed (BED format converter)")
        logger.info("  • trf-large (modified TRF for large genomes)")
        logger.info("")
        logger.info("Note: Installation may fail on some platforms (especially macOS).")
        logger.info("      Satellome will work fine even if these tools fail to install.")
        logger.info("      You can install them later with: satellome --install-all")
        logger.info("")
        logger.info("To skip this step next time, set: SATELLOME_SKIP_AUTO_INSTALL=1")
        logger.info("="*60 + "\n")

        # Try to import and run installers
        try:
            from satellome.installers import install_fastan, install_tanbed, install_trf_large

            success_count = 0
            failed_tools = []

            # Install FasTAN
            logger.info("[1/3] Installing FasTAN...")
            try:
                if install_fastan(force=False):
                    logger.info("✓ FasTAN installed successfully\n")
                    success_count += 1
                else:
                    logger.warning("✗ FasTAN installation failed\n")
                    failed_tools.append("FasTAN")
            except Exception as e:
                logger.warning(f"✗ FasTAN installation failed: {e}\n")
                failed_tools.append("FasTAN")

            # Install tanbed
            logger.info("[2/3] Installing tanbed...")
            try:
                if install_tanbed(force=False):
                    logger.info("✓ tanbed installed successfully\n")
                    success_count += 1
                else:
                    logger.warning("✗ tanbed installation failed\n")
                    failed_tools.append("tanbed")
            except Exception as e:
                logger.warning(f"✗ tanbed installation failed: {e}\n")
                failed_tools.append("tanbed")

            # Install modified TRF
            logger.info("[3/3] Installing modified TRF (for large genomes)...")
            try:
                if install_trf_large(force=False):
                    logger.info("✓ Modified TRF installed successfully\n")
                    success_count += 1
                else:
                    logger.warning("✗ Modified TRF installation failed\n")
                    failed_tools.append("trf-large")
            except Exception as e:
                logger.warning(f"✗ Modified TRF installation failed: {e}\n")
                failed_tools.append("trf-large")

            # Print summary
            logger.info("="*60)
            logger.info("Satellome installation summary:")
            logger.info(f"  Successfully installed: {success_count}/3 tools")

            if failed_tools:
                logger.info(f"  Failed: {', '.join(failed_tools)}")
                logger.info("")
                logger.info("To retry installation:")
                logger.info("  satellome --install-all")
                logger.info("")
                logger.info("Or install individually:")
                for tool in failed_tools:
                    cmd = tool.lower().replace('-', '_')
                    logger.info(f"  satellome --install-{tool.lower()}")
            else:
                logger.info("  All tools installed successfully!")

            logger.info("="*60 + "\n")

        except ImportError as e:
            logger.warning(f"\nCould not import installers: {e}")
            logger.warning("External tools were not installed automatically.")
            logger.warning("You can install them later with: satellome --install-all\n")
        except Exception as e:
            logger.warning(f"\nUnexpected error during tool installation: {e}")
            logger.warning("External tools were not installed automatically.")
            logger.warning("You can install them later with: satellome --install-all\n")
