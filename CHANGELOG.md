# Changelog

All notable changes to Satellome will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.4.1] - 2025-01-05

### Fixed
- **Critical bugfix**: Removed broken import fallback pattern in main.py that caused ModuleNotFoundError when package was installed via pip
- Removed sys.path.append hack that interfered with proper package installation
- Fixed SyntaxWarning: invalid escape sequences in regex patterns (Python 3.12+ compatibility)
- Package now imports correctly when installed from PyPI

### Changed
- **Removed external dependencies**: Replaced `requests` with `urllib` from stdlib
  - Eliminates dependency conflicts with `urllib3` and `six`
  - Reduces total dependencies from 13 to 11 packages
  - No impact on functionality - NCBI taxon name fetching still works correctly
  - More reliable installation across different environments

## [1.4.0] - 2025-01-04

### Added
- **Automatic tool installation during pip install** (SAT-37)
  - External tools (FasTAN, tanbed, modified TRF) now install automatically with `pip install satellome`
  - Binaries installed to `<site-packages>/satellome/bin/` instead of `~/.satellome/bin/` (cleaner)
  - Graceful failure: Satellome installs successfully even if tool compilation fails
  - Can be skipped with `SATELLOME_SKIP_AUTO_INSTALL=1` environment variable
  - Fallback to `~/.satellome/bin/` if no write permissions to site-packages

- **Automatic installer for FasTAN and tanbed** (SAT-35)
  - CLI commands: `--install-fastan`, `--install-tanbed`, `--install-all`
  - Automatic compilation from source with dependency checking
  - Works on Linux and macOS
  - FasTAN: alternative tandem repeat finder by Gene Myers
  - tanbed: BED format converter from alntools

- **Automatic installer for modified TRF** (SAT-36)
  - CLI command: `--install-trf-large`
  - Modified TRF from aglabx for large genomes (chromosomes >2GB)
  - Supports large plant and amphibian genomes
  - Best support on Linux; manual installation recommended for macOS

### Changed
- Binary location priority changed: `<site-packages>/satellome/bin/` is now primary location
- Installation process significantly simplified - one command does everything
- No pollution of user's home directory by default

### Fixed
- External tool installation no longer blocks Satellome installation on failure

## [1.3.0] - Previous release

Earlier changes not documented in this format.
