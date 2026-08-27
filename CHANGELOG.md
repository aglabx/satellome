# Changelog

All notable changes to Satellome will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.10.0] - 2026-08-28

### Added
- **`satellome --fix-path`, and the same repair automatically**: when the
  launcher's directory is missing from `PATH`, satellome now adds it instead of
  only explaining how. One marked block (`# >>> satellome path >>>`) is appended
  to the startup file the user's shell already reads — bash, zsh (honouring
  `ZDOTDIR`), and fish (`fish_add_path`) — written in `$HOME/...` form so the
  file stays portable between machines. The edit is idempotent via its marker,
  does not duplicate an entry added by hand, and prefers a startup file that
  already exists: creating `~/.bash_profile` where none existed would stop a
  login shell from reading `~/.profile`.
  Since no process can change its parent shell's environment, the `source ...`
  command needed for the *current* shell is always printed.
  `SATELLOME_NO_PATH_FIX=1` keeps the diagnosis but never edits a file;
  `SATELLOME_NO_ENV_CHECK=1` disables the check entirely.

### Changed
- A shadowed launcher (a `satellome` on `PATH` belonging to a different install
  than the active interpreter) is reported but deliberately **not** auto-fixed:
  choosing which of two installations to remove is the user's decision, and no
  `PATH` edit can express it.

## [1.9.0] - 2026-08-28

### Added
- **`satellome --doctor`**: diagnoses the installation itself — where the
  package lives, which interpreter runs it, whether the shell can actually see
  the `satellome` launcher, which script directories are on `PATH`, and where
  every external tool resolves from. Exits 0 when healthy and 1 when it finds a
  problem, so a driver script can gate a setup on it.
- **Startup PATH check**: every run (and the no-argument info screen) now warns
  when the launcher is missing from `PATH` or is shadowed by a competing
  install running a different interpreter. `pip install --user` puts the
  launcher in `~/.local/bin`, which is frequently absent from `PATH`; the
  warning names the directory and gives the three fixes. Silence it with
  `SATELLOME_NO_ENV_CHECK=1` (`--doctor` still reports it).

### Fixed
- **Companion tools installed off `PATH` are no longer reported as missing.**
  `arraysplitter` was resolved with bare `shutil.which`, so on a machine whose
  `PATH` lacks the user scripts directory the tool was declared absent, an
  automatic `pip install` "succeeded", the second lookup failed again, and the
  step gave up with "arraysplitter installed but binary not found in PATH" —
  every run, while the executable sat in `~/.local/bin`. Resolution now falls
  back to the running interpreter's own scripts directories
  (`core_functions/tools/env_check.py`), so the tool is used by absolute path
  and the broken shell setup is reported instead of silently degrading the run.

## [1.8.0] - 2026-08-17

### Added
- **Run manifest and `satellome --verify-run <dir>`**: every run writes
  `run_manifest.json` **last**, listing each produced file with its byte size and
  the status of every step. `--verify-run` exits 0 only for a verifiably complete
  run; 1 for a missing or corrupt manifest, a failed step, a missing file, a
  leftover `*.partial`, or a file whose size no longer matches what the run
  wrote; 2 for a bad argument. Existence checks cannot detect truncation — this
  can. For `.gz` copies the gzip ISIZE trailer is compared to the recorded size,
  which catches an archive made from a partially read file (such an archive is
  structurally valid and passes `gzip -t`).
- **Output-directory lock**: a second run into the same `-o` is refused, naming
  the pid and host that holds the lock. Override explicitly with `--ignore-lock`.

### Changed
- **Atomic output publication**: every output file is written to
  `<path>.partial` and `os.replace`-d into place, so a final name never refers to
  a half-written file — `.sat` and FASTA extraction outputs, the gaps BED, the
  karyotype charts and the manifest itself. A driver script that compresses the
  output directory concurrently now gets either the complete file or no file.
- A failing drawing step no longer exits on the spot: the manifest is still
  written, with `drawing: failed`, and satellome exits 1. The failure is recorded
  in the run's own artifact instead of only in an exit code.
- The drawing step is invoked as an argument list instead of through a shell, so
  quotes in an organism name cannot mangle the command line.
- `atomic_io` reports when it has to create an output file's parent directory:
  the pipeline creates its directories up front, so a missing parent means a
  stray path separator in a file name.

### Fixed
- **Drawing crash on organism names containing a slash**: NCBI strain
  designations such as `MHOM/BR/75/M2904` or `CCAP 1055/1` went into the
  karyotype file name unchanged, so the slash was read as a path separator and
  drawing died with `FileNotFoundError` — after the whole pipeline had already
  written its data. Taxon names are now reduced to a single safe file-name
  component and the substitution is logged; plot titles keep the original name.
- A failure to write the gaps BED file propagates instead of being logged and
  ignored, and a gaps BED that exists but cannot be parsed is rewritten instead
  of being reused as cache by the next run.
- `redraw.py` no longer passes `-d` to `trf_draw.py`, which does not accept it.

## [1.5.2] - 2025-01-27

### Changed
- **FasTAN is now the default tool**: TRF no longer runs by default
  - Use `--run-trf` flag to enable TRF analysis
  - FasTAN runs by default (use `--nofastan` to disable)
  - `--notrf` flag is deprecated (TRF is already disabled by default)

## [1.5.1] - 2025-01-17

### Fixed
- **Fallback to system gcc when conda gcc fails**: When compiling tanbed or FasTAN in conda environments, if conda's gcc fails with "cannot find -lz" error, the installer now automatically retries with system gcc (/usr/bin/gcc)

## [1.5.0] - 2025-01-15

### Added
- **Comprehensive input validation** (SAT-5): Validation for FASTA, GFF, RepeatMasker, and TRF binary files
- **Memory-efficient streaming annotation** (SAT-7): Support for large genomes without loading entire files into memory
- **Comprehensive docstrings** (SAT-10): Added detailed documentation to all core modules:
  - I/O module (FASTA, GFF3, TRF, TAB handlers)
  - Tools module (clustering, statistics, HTML reports)
  - Analysis and utility tools
  - TRF model
- **Informative headers in TRF output files**: Files now include metadata headers with generation info
- **Size-filtered TRF files**: FasTAN output now includes filtered files (1kb, 3kb, 10kb thresholds)
- **Repeat length field** (SAT-49): Added to output TRF files for easier analysis
- **Sequence extraction**: FasTAN BED output now includes extracted sequences
- **Duplicate chromosome name validation** (SAT-49): Prevents issues with duplicate names in input

### Changed
- **Reduced TRF output fields from 37 to 18** (SAT-13): Streamlined output with automatic conversion for compatibility
- **FasTAN output now in TRF-compatible format**: With FASTA export support
- **Removed .pkl cache files** (SAT-47): Now uses only BED format for gaps data
- **Short chromosome names in output** (SAT-49): Cleaner output format
- **All file paths converted to absolute**: Prevents path resolution issues
- **Improved taxid resolution logging**: Clearer messages during taxonomy lookup
- **Removed outdated FASTA line length warning**

### Fixed
- **Atomic file operations** (SAT-54): `move_files_by_mask` now uses atomic operations to prevent data loss
- **Improved error handling** (SAT-22): Replaced generic exceptions and assertions with specific error types
- **TRF validation skip** when `--notrf` flag is set
- **Chromosome name parsing in BED files** (SAT-49): Fixed handling of names with spaces
- **Test fixes** (SAT-58): Fixed TRF file format in tests
- **Test compatibility**: Handle TRF header comments correctly

### Documentation
- **Decomposed TODO comments** (SAT-8): Converted inline TODOs into tracked issues
- Removed placeholder Zenodo DOI badge from README

## [1.4.3] - 2025-01-05

### Added
- **FasTAN integration** (SAT-24): Alternative tandem repeat finder by Gene Myers
  - FasTAN runs by default alongside TRF for comprehensive repeat analysis
  - Automatic installation with `--install-fastan` or during `pip install`
  - Use `--nofastan` to skip FasTAN and run TRF only
  - Use `--notrf` to skip TRF and run FasTAN only
  - Results saved to `fastan/` directory (`.1aln` and `.bed` formats)
  - Automatic tanbed installation for BED format conversion
- **Gaps annotation export to BED format** (SAT-41)
  - Genomic gaps automatically exported to `{project_name}.gaps.bed`
  - BED6 format with gap length as score
  - Includes detailed header with project, taxon, and total gap count
- **Professional README badges** (SAT-46)
  - CI/CD status (GitHub Actions)
  - Code coverage (Codecov)
  - Python version support (3.9, 3.10, 3.11)
  - MIT License badge
  - PyPI version badge
  - DOI badge for citations
- **GitHub Actions CI/CD** (SAT-42)
  - Automated testing on push and pull requests
  - Matrix testing across Python 3.9, 3.10, 3.11
  - Code coverage reporting to Codecov

### Fixed
- **Gaps BED file naming** (SAT-44): Now uses project name instead of taxon
  - Previous: `Unknown.gaps.bed` (using taxon)
  - Current: `{project_name}.gaps.bed` (using project name from input file)
- **Critical test fixes** (SAT-43): Fixed 12 failing unit tests
  - Overlap merging calculation fixed (coordinate-based, not array-based)
  - TRF format validation and parsing
  - Mock object handling in tests
- **Additional test fixes** (SAT-45): Fixed 10 more failing unit tests
  - TRF parse head now returns None for invalid input (proper handling in TRModel)
  - Fixed test fixtures to match real TRF format (15 fields, proper header format)
  - Fixed coordinate-length mismatches (1-based inclusive coordinates)
- **FasTAN output file naming**: Now uses genome filename instead of project name for consistency
- **FasTAN installer improvements**:
  - Better binary detection in PATH and satellome bin directories
  - Fixed compilation with pthread linker flag
  - Uses forked FasTAN repository with fixed Makefile

### Changed
- **Default behavior**: Both TRF and FasTAN run by default for comprehensive analysis
  - Use `--nofastan` to run TRF only
  - Use `--notrf` to run FasTAN only
  - At least one tool must run
- FasTAN and tanbed now auto-install on first use (like TRF)

## [1.4.1] - 2025-01-05

### Added
- **Automatic TRF installation on first run with smart fallback**:
  - Modified TRF tries to install first (for genomes with chromosomes >2GB)
    - Now uses simple `build.sh` - only needs git and gcc/clang (no automake/autoconf!)
  - Falls back to standard TRF if build tools are missing (downloads pre-compiled binary)
  - Works in any environment - containers, minimal installs, etc.
  - Works with both regular install and editable install (`pip install -e .`)
  - No manual installation needed - just run `satellome`
- **New installer**: `install_trf_standard()` - downloads pre-compiled TRF binary
  - Available via `--install-trf` flag
  - No build tools required
  - Works on Linux and macOS

### Fixed
- **Critical bugfix**: Removed broken import fallback pattern in main.py that caused ModuleNotFoundError when package was installed via pip
- Removed sys.path.append hack that interfered with proper package installation
- Fixed SyntaxWarning: invalid escape sequences in regex patterns (Python 3.12+ compatibility)
- Fixed subprocess calls to use correct Python interpreter (sys.executable instead of "python")
  - Fixes ModuleNotFoundError in subprocesses when system has multiple Python versions
- Removed interactive input() prompt in temp folder cleanup that blocked automatic pipeline execution
- Added explicit six>=1.16.0 dependency to fix pandas/python-dateutil compatibility in containers
  - Fixes "ModuleNotFoundError: No module named 'six.moves'" when system has outdated six package
- Package now imports correctly when installed from PyPI

### Changed
- **Removed external dependencies**: Replaced `requests` with `urllib` from stdlib
  - Eliminates dependency conflicts with `urllib3` and `six`
  - Reduces total dependencies from 13 to 11 packages
  - No impact on functionality - NCBI taxon name fetching still works correctly
  - More reliable installation across different environments
- **Eliminated kaleido/chromium dependency completely** (SAT-29): Now uses matplotlib for static images
  - Karyotype plots saved in BOTH formats: interactive HTML (plotly) + static PNG (matplotlib)
  - No external browser or chromium dependency required
  - Works in any environment: containers, minimal installs, headless servers
  - matplotlib backend (Agg) handles PNG/SVG export without X11/display
  - Significantly simpler installation and smaller footprint
- **Removed pandas dependency** (SAT-40): Replaced with stdlib csv module and native Python data structures
  - **Completely eliminates six.moves compatibility issues in containers**
  - Reduces dependencies by ~20MB
  - Faster module imports
  - DataFrame operations replaced with list of dicts (simpler, more Pythonic)
  - pd.read_csv() → csv.DictReader() for TRF file parsing
  - All data filtering/transformation uses native Python (list comprehensions, dict operations)
  - No datetime functionality was used, so no loss of features

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
