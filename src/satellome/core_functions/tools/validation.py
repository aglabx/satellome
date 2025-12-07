#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 2025-11-15
# @author: Claude Code
# @contact: noreply@anthropic.com

"""
Comprehensive input validation for Satellome pipeline.

Provides robust validation functions for all input files and system requirements
before running tandem repeat analysis. Catches common errors early with helpful
error messages and installation instructions.

Classes:
    ValidationError: Base exception for validation failures
    FastaValidationError: FASTA file validation failures
    GFFValidationError: GFF3 file validation failures
    BinaryValidationError: TRF binary validation failures
    OutputDirValidationError: Output directory validation failures

Functions:
    validate_fasta_file: Validate FASTA format, sequences, and content
    validate_gff_file: Validate GFF3 format and coordinates
    validate_repeatmasker_file: Validate RepeatMasker .out format
    validate_trf_binary: Validate TRF executable exists and is runnable
    validate_output_directory: Validate/create output directory with permissions check
    validate_input_files: Main validation function - validates all inputs at once

Key Features:
    - Early error detection before pipeline execution
    - Comprehensive format validation (headers, coordinates, fields)
    - File size and emptiness checks
    - Gzip-compressed FASTA support
    - Detailed error messages with fix suggestions
    - Installation instructions for missing TRF binary
    - Disk space warnings for output directory
    - IUPAC ambiguity code support in sequences
    - Malformed line detection and reporting

Validation Checks:
    **FASTA Files:**
    - File exists and is readable
    - Not empty (non-zero size)
    - Valid headers (start with '>')
    - Sequence content validation (IUPAC codes)
    - Zero-length sequence detection
    - Invalid character detection
    - Gzip format validation

    **GFF Files:**
    - 9 tab-separated fields per line
    - Valid integer coordinates
    - Start <= end positions
    - 1-based coordinate system check

    **RepeatMasker Files:**
    - At least 11 whitespace-separated fields
    - Valid integer coordinates (fields 5-6)
    - Start <= end positions

    **TRF Binary:**
    - File exists (full path or in $PATH)
    - Is executable (chmod +x)
    - Provides installation instructions if missing

    **Output Directory:**
    - Directory exists or can be created
    - Write permissions available
    - Sufficient disk space (warns if < 1GB)

Example:
    >>> from satellome.core_functions.tools.validation import validate_input_files
    >>> results = validate_input_files(
    ...     fasta_file="genome.fasta",
    ...     gff_file="annotations.gff3",
    ...     trf_binary="trf",
    ...     output_dir="output/"
    ... )
    INFO:...Validating FASTA file: genome.fasta
    INFO:...✓ FASTA validation passed: 24 sequences, 3000000000 bp total
    INFO:...Validating GFF file: annotations.gff3
    INFO:...✓ GFF validation passed: 12453 features
    INFO:...Validating TRF binary: trf
    INFO:...✓ TRF binary found: /usr/local/bin/trf
    INFO:...Validating output directory: output/
    INFO:...✓ Output directory ready: /path/to/output
    INFO:...All validations passed successfully!
    >>>
    >>> # Check statistics
    >>> print(results['fasta']['num_sequences'])
    24

Error Handling Example:
    >>> try:
    ...     validate_fasta_file("missing.fasta")
    ... except FastaValidationError as e:
    ...     print(f"Error: {e}")
    Error: FASTA file not found: missing.fasta

Typical Use Case:
    1. Call validate_input_files() at pipeline start
    2. Catch validation exceptions and display error to user
    3. If successful, proceed with analysis using validated paths
    4. Results dict contains statistics for logging/reporting

See Also:
    satellome.main: Pipeline entry point that uses these validators
    satellome.core_functions.io.fasta_file: FASTA parsing after validation
    satellome.core_functions.io.gff3_file: GFF3 parsing after validation
"""

import os
import sys
import shutil
import logging
import gzip
from pathlib import Path

logger = logging.getLogger(__name__)


class ValidationError(Exception):
    """Base exception for validation errors."""
    pass


class FastaValidationError(ValidationError):
    """Exception raised for FASTA validation errors."""
    pass


class GFFValidationError(ValidationError):
    """Exception raised for GFF validation errors."""
    pass


class BinaryValidationError(ValidationError):
    """Exception raised for binary validation errors."""
    pass


class OutputDirValidationError(ValidationError):
    """Exception raised for output directory validation errors."""
    pass


def validate_fasta_file(fasta_path, check_sequences=True):
    """
    Validate FASTA file format and content.

    Args:
        fasta_path: Path to FASTA file
        check_sequences: If True, validate sequence content (default: True)

    Raises:
        FastaValidationError: If validation fails

    Returns:
        dict: Validation statistics (num_sequences, total_length, warnings)
    """
    if not os.path.exists(fasta_path):
        raise FastaValidationError(f"FASTA file not found: {fasta_path}")

    if not os.path.isfile(fasta_path):
        raise FastaValidationError(f"FASTA path is not a file: {fasta_path}")

    file_size = os.path.getsize(fasta_path)
    if file_size == 0:
        raise FastaValidationError(f"FASTA file is empty: {fasta_path}")

    # Determine opener (gzip or regular)
    if fasta_path.endswith('.gz'):
        opener = gzip.open
        mode = 'rt'
    else:
        opener = open
        mode = 'r'

    num_sequences = 0
    total_length = 0
    warnings = []
    current_seq_name = None
    current_seq_length = 0
    line_num = 0
    has_valid_header = False

    # Valid nucleotide characters (including IUPAC ambiguity codes)
    valid_chars = set('ACGTURYSWKMBDHVNacgturyswkmbdhvn-.')

    try:
        with opener(fasta_path, mode) as fh:
            for line in fh:
                line_num += 1
                line = line.rstrip('\n\r')

                if not line:  # Skip empty lines
                    continue

                if line.startswith('>'):
                    has_valid_header = True

                    # Save previous sequence stats
                    if current_seq_name:
                        if current_seq_length == 0:
                            warnings.append(f"Sequence '{current_seq_name}' has zero length")
                        total_length += current_seq_length

                    # Start new sequence
                    num_sequences += 1
                    current_seq_name = line[1:].split()[0] if len(line) > 1 else f"seq_{num_sequences}"
                    current_seq_length = 0

                    if len(line) == 1:  # Just ">" with no name
                        warnings.append(f"Line {line_num}: Header has no sequence name")

                elif check_sequences:
                    # Validate sequence line
                    if not current_seq_name:
                        raise FastaValidationError(
                            f"Line {line_num}: Sequence data before first header"
                        )

                    # Check for invalid characters
                    invalid_chars = set(line) - valid_chars
                    if invalid_chars:
                        warnings.append(
                            f"Line {line_num}: Invalid characters in sequence '{current_seq_name}': "
                            f"{', '.join(sorted(invalid_chars))}"
                        )

                    current_seq_length += len(line)

        # Process last sequence
        if current_seq_name:
            if current_seq_length == 0:
                warnings.append(f"Sequence '{current_seq_name}' has zero length")
            total_length += current_seq_length

    except (IOError, OSError) as e:
        raise FastaValidationError(f"Error reading FASTA file: {e}")
    except gzip.BadGzipFile:
        raise FastaValidationError(f"File appears to be corrupted (bad gzip format): {fasta_path}")
    except UnicodeDecodeError as e:
        raise FastaValidationError(f"File contains invalid text encoding at line ~{line_num}: {e}")

    # Final validation checks
    if not has_valid_header:
        raise FastaValidationError(
            f"No valid FASTA headers found (headers must start with '>')"
        )

    if num_sequences == 0:
        raise FastaValidationError("No sequences found in FASTA file")

    if total_length == 0:
        raise FastaValidationError("Total sequence length is zero")

    return {
        'num_sequences': num_sequences,
        'total_length': total_length,
        'file_size': file_size,
        'warnings': warnings
    }


def validate_gff_file(gff_path, check_format=True):
    """
    Validate GFF3 file format.

    Args:
        gff_path: Path to GFF file
        check_format: If True, validate GFF format (default: True)

    Raises:
        GFFValidationError: If validation fails

    Returns:
        dict: Validation statistics (num_features, num_malformed, warnings)
    """
    if not os.path.exists(gff_path):
        raise GFFValidationError(f"GFF file not found: {gff_path}")

    if not os.path.isfile(gff_path):
        raise GFFValidationError(f"GFF path is not a file: {gff_path}")

    file_size = os.path.getsize(gff_path)
    if file_size == 0:
        raise GFFValidationError(f"GFF file is empty: {gff_path}")

    num_features = 0
    num_malformed = 0
    warnings = []
    line_num = 0

    try:
        with open(gff_path, 'r') as fh:
            for line in fh:
                line_num += 1
                line = line.strip()

                # Skip empty lines and comments
                if not line or line.startswith('#'):
                    continue

                if check_format:
                    fields = line.split('\t')

                    # GFF3 requires 9 tab-separated fields
                    if len(fields) < 9:
                        num_malformed += 1
                        if num_malformed <= 10:  # Only report first 10
                            warnings.append(
                                f"Line {line_num}: Expected 9 fields, got {len(fields)}"
                            )
                    else:
                        # Validate coordinate fields (4th and 5th fields)
                        try:
                            start = int(fields[3])
                            end = int(fields[4])
                            if start > end:
                                warnings.append(
                                    f"Line {line_num}: Start position ({start}) > end position ({end})"
                                )
                            if start < 1:
                                warnings.append(
                                    f"Line {line_num}: Start position ({start}) < 1 (GFF uses 1-based coords)"
                                )
                        except ValueError:
                            warnings.append(
                                f"Line {line_num}: Invalid coordinate values (must be integers)"
                            )

                        num_features += 1
                else:
                    num_features += 1

    except (IOError, OSError) as e:
        raise GFFValidationError(f"Error reading GFF file: {e}")
    except UnicodeDecodeError as e:
        raise GFFValidationError(f"File contains invalid text encoding at line ~{line_num}: {e}")

    if num_features == 0 and num_malformed == 0:
        raise GFFValidationError("No valid features found in GFF file")

    if num_malformed > num_features * 0.1:  # More than 10% malformed
        warnings.append(
            f"WARNING: High proportion of malformed lines ({num_malformed}/{num_features + num_malformed})"
        )

    return {
        'num_features': num_features,
        'num_malformed': num_malformed,
        'file_size': file_size,
        'warnings': warnings
    }


def validate_repeatmasker_file(rm_path, check_format=True):
    """
    Validate RepeatMasker .out file format.

    Args:
        rm_path: Path to RepeatMasker file
        check_format: If True, validate format (default: True)

    Raises:
        ValidationError: If validation fails

    Returns:
        dict: Validation statistics (num_features, num_malformed, warnings)
    """
    if not os.path.exists(rm_path):
        raise ValidationError(f"RepeatMasker file not found: {rm_path}")

    if not os.path.isfile(rm_path):
        raise ValidationError(f"RepeatMasker path is not a file: {rm_path}")

    file_size = os.path.getsize(rm_path)
    if file_size == 0:
        raise ValidationError(f"RepeatMasker file is empty: {rm_path}")

    num_features = 0
    num_malformed = 0
    warnings = []
    line_num = 0

    try:
        with open(rm_path, 'r') as fh:
            for line in fh:
                line_num += 1
                line = line.strip()

                # Skip empty lines and header lines
                if not line or line.startswith('SW') or line.startswith('score'):
                    continue

                if check_format:
                    fields = line.split()

                    # RepeatMasker .out requires at least 11 fields
                    # Fields we use: [4]=chrm, [5]=start, [6]=end, [10]=family
                    if len(fields) < 11:
                        num_malformed += 1
                        if num_malformed <= 10:  # Only report first 10
                            warnings.append(
                                f"Line {line_num}: Expected ≥11 fields, got {len(fields)}"
                            )
                    else:
                        # Validate coordinate fields
                        try:
                            start = int(fields[5])
                            end = int(fields[6])
                            if start > end:
                                warnings.append(
                                    f"Line {line_num}: Start position ({start}) > end position ({end})"
                                )
                        except (ValueError, IndexError):
                            num_malformed += 1
                            if num_malformed <= 10:
                                warnings.append(
                                    f"Line {line_num}: Invalid coordinate values"
                                )

                        num_features += 1
                else:
                    num_features += 1

    except (IOError, OSError) as e:
        raise ValidationError(f"Error reading RepeatMasker file: {e}")
    except UnicodeDecodeError as e:
        raise ValidationError(f"File contains invalid text encoding at line ~{line_num}: {e}")

    if num_features == 0 and num_malformed == 0:
        raise ValidationError("No valid features found in RepeatMasker file")

    return {
        'num_features': num_features,
        'num_malformed': num_malformed,
        'file_size': file_size,
        'warnings': warnings
    }


def validate_trf_binary(trf_path):
    """
    Validate that TRF binary exists and is executable.

    Args:
        trf_path: Path to TRF binary or command name

    Raises:
        BinaryValidationError: If validation fails

    Returns:
        str: Full path to TRF binary
    """
    # Check if it's a full path
    if os.path.sep in trf_path:
        if not os.path.exists(trf_path):
            raise BinaryValidationError(
                f"TRF binary not found at specified path: {trf_path}\n"
                f"Install TRF with:\n"
                f"  - satellome --install-trf-large (for large genomes, requires build tools)\n"
                f"  - satellome --install-trf (standard version, pre-compiled)\n"
                f"  - Download manually from: https://tandem.bu.edu/trf/trf.html"
            )

        if not os.path.isfile(trf_path):
            raise BinaryValidationError(f"TRF path is not a file: {trf_path}")

        if not os.access(trf_path, os.X_OK):
            raise BinaryValidationError(
                f"TRF binary is not executable: {trf_path}\n"
                f"Fix with: chmod +x {trf_path}"
            )

        return os.path.abspath(trf_path)

    # It's a command name - search in PATH
    trf_found = shutil.which(trf_path)
    if not trf_found:
        raise BinaryValidationError(
            f"TRF binary '{trf_path}' not found in PATH.\n"
            f"Install TRF with:\n"
            f"  - satellome --install-trf-large (for large genomes, requires build tools)\n"
            f"  - satellome --install-trf (standard version, pre-compiled)\n"
            f"  - Download manually from: https://tandem.bu.edu/trf/trf.html"
        )

    return trf_found


def validate_output_directory(output_path, create_if_missing=True):
    """
    Validate output directory and check write permissions.

    Args:
        output_path: Path to output directory
        create_if_missing: If True, create directory if it doesn't exist (default: True)

    Raises:
        OutputDirValidationError: If validation fails

    Returns:
        str: Absolute path to output directory
    """
    output_path = os.path.abspath(output_path)

    # Check if path exists
    if os.path.exists(output_path):
        # Check if it's a directory
        if not os.path.isdir(output_path):
            raise OutputDirValidationError(
                f"Output path exists but is not a directory: {output_path}"
            )

        # Check write permissions
        if not os.access(output_path, os.W_OK):
            raise OutputDirValidationError(
                f"No write permission for output directory: {output_path}\n"
                f"Fix with: chmod u+w {output_path}"
            )
    else:
        # Directory doesn't exist
        if not create_if_missing:
            raise OutputDirValidationError(
                f"Output directory does not exist: {output_path}"
            )

        # Try to create it
        try:
            os.makedirs(output_path, exist_ok=True)
            logger.info(f"Created output directory: {output_path}")
        except PermissionError:
            parent_dir = os.path.dirname(output_path)
            raise OutputDirValidationError(
                f"Permission denied: Cannot create output directory: {output_path}\n"
                f"Check permissions for parent directory: {parent_dir}"
            )
        except OSError as e:
            if "No space left on device" in str(e):
                raise OutputDirValidationError(
                    f"Disk full: Cannot create output directory: {output_path}\n"
                    f"Free up disk space and try again."
                )
            elif "File name too long" in str(e):
                raise OutputDirValidationError(
                    f"Path too long: Cannot create output directory: {output_path}\n"
                    f"Use a shorter path."
                )
            else:
                raise OutputDirValidationError(
                    f"Cannot create output directory: {output_path}\n"
                    f"Error: {e}"
                )

    # Check available disk space (warn if < 1GB)
    try:
        stat = os.statvfs(output_path)
        available_bytes = stat.f_bavail * stat.f_frsize
        available_gb = available_bytes / (1024 ** 3)

        if available_gb < 1:
            logger.warning(
                f"Low disk space: Only {available_gb:.2f} GB available in {output_path}\n"
                f"Large genome analysis may require several GB of disk space."
            )
    except (AttributeError, OSError):
        # statvfs not available on all platforms (e.g., Windows)
        pass

    return output_path


def validate_input_files(fasta_file, gff_file=None, rm_file=None,
                        trf_binary="trf", output_dir=None,
                        check_sequences=True):
    """
    Validate all input files and directories for Satellome pipeline.

    Args:
        fasta_file: Path to FASTA file (required)
        gff_file: Path to GFF file (optional)
        rm_file: Path to RepeatMasker file (optional)
        trf_binary: Path to TRF binary or command name (default: "trf")
        output_dir: Path to output directory (optional)
        check_sequences: If True, validate FASTA sequence content (default: True)

    Raises:
        ValidationError: If any validation fails

    Returns:
        dict: Validation results for all inputs
    """
    results = {}

    logger.info("Validating input files...")

    # Validate FASTA file
    try:
        logger.info(f"Validating FASTA file: {fasta_file}")
        fasta_stats = validate_fasta_file(fasta_file, check_sequences=check_sequences)
        results['fasta'] = fasta_stats

        logger.info(
            f"✓ FASTA validation passed: "
            f"{fasta_stats['num_sequences']} sequences, "
            f"{fasta_stats['total_length']:,} bp total"
        )

        if fasta_stats['warnings']:
            logger.warning(f"FASTA validation warnings ({len(fasta_stats['warnings'])}):")
            for warning in fasta_stats['warnings'][:5]:  # Show first 5
                logger.warning(f"  - {warning}")
            if len(fasta_stats['warnings']) > 5:
                logger.warning(f"  ... and {len(fasta_stats['warnings']) - 5} more warnings")

    except FastaValidationError as e:
        logger.error(f"✗ FASTA validation failed: {e}")
        raise

    # Validate GFF file (if provided)
    if gff_file:
        try:
            logger.info(f"Validating GFF file: {gff_file}")
            gff_stats = validate_gff_file(gff_file)
            results['gff'] = gff_stats

            logger.info(
                f"✓ GFF validation passed: "
                f"{gff_stats['num_features']} features"
            )

            if gff_stats['num_malformed'] > 0:
                logger.warning(
                    f"GFF file contains {gff_stats['num_malformed']} malformed lines "
                    f"(will be skipped during processing)"
                )

        except GFFValidationError as e:
            logger.error(f"✗ GFF validation failed: {e}")
            raise

    # Validate RepeatMasker file (if provided)
    if rm_file:
        try:
            logger.info(f"Validating RepeatMasker file: {rm_file}")
            rm_stats = validate_repeatmasker_file(rm_file)
            results['repeatmasker'] = rm_stats

            logger.info(
                f"✓ RepeatMasker validation passed: "
                f"{rm_stats['num_features']} features"
            )

            if rm_stats['num_malformed'] > 0:
                logger.warning(
                    f"RepeatMasker file contains {rm_stats['num_malformed']} malformed lines "
                    f"(will be skipped during processing)"
                )

        except ValidationError as e:
            logger.error(f"✗ RepeatMasker validation failed: {e}")
            raise

    # Validate TRF binary
    try:
        logger.info(f"Validating TRF binary: {trf_binary}")
        trf_path = validate_trf_binary(trf_binary)
        results['trf_binary'] = {'path': trf_path}
        logger.info(f"✓ TRF binary found: {trf_path}")

    except BinaryValidationError as e:
        logger.error(f"✗ TRF binary validation failed: {e}")
        raise

    # Validate output directory
    if output_dir:
        try:
            logger.info(f"Validating output directory: {output_dir}")
            output_path = validate_output_directory(output_dir)
            results['output_dir'] = {'path': output_path}
            logger.info(f"✓ Output directory ready: {output_path}")

        except OutputDirValidationError as e:
            logger.error(f"✗ Output directory validation failed: {e}")
            raise

    logger.info("All validations passed successfully!")
    return results
