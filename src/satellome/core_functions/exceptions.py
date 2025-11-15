#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 2025-11-15
# @author: Claude Code
# @contact: noreply@anthropic.com

"""
Custom exception classes for Satellome.

This module defines a hierarchy of exception classes for better error handling
and debugging throughout the Satellome pipeline. Using specific exception types
instead of generic Exception improves:

1. Error diagnosis - specific exception names indicate the error category
2. Error handling - allows catching specific error types
3. User experience - clearer error messages with context and suggestions
4. Testing - enables testing specific error conditions

Exception Hierarchy:
    SatellomeError (base)
    ├── StatisticsError - statistical computation failures
    ├── SequenceError - invalid sequence data or operations
    ├── FileFormatError - file format/parsing issues
    ├── ConfigurationError - configuration or setup issues
    └── ValidationError - input validation failures (aliased from validation.py)
"""


class SatellomeError(Exception):
    """
    Base exception for all Satellome-specific errors.

    All custom exceptions in Satellome inherit from this class, allowing
    code to catch all Satellome-specific errors with a single except clause
    if needed, while still allowing fine-grained error handling.

    Example:
        try:
            process_genome(fasta_file)
        except SatellomeError as e:
            logger.error(f"Satellome processing failed: {e}")
        except Exception as e:
            logger.error(f"Unexpected error: {e}")
            raise
    """
    pass


class StatisticsError(SatellomeError):
    """
    Raised when statistical computation fails.

    This includes:
    - Empty data arrays
    - Invalid variance/standard deviation values
    - Invalid sample sizes
    - Numerical computation issues (overflow, underflow)

    Example:
        raise StatisticsError(
            f"Cannot compute variance: empty data array. "
            f"Ensure the input contains at least one valid data point."
        )
    """
    pass


class SequenceError(SatellomeError):
    """
    Raised for invalid sequence data or sequence operations.

    This includes:
    - Invalid nucleotide sequences
    - Sequence length mismatches
    - Invalid monomer sequences
    - Empty consensus sequences
    - Array length != coordinate range

    Example:
        raise SequenceError(
            f"Invalid monomer sequence '{monomer}': contains non-ACGT characters. "
            f"Valid monomers must contain only A, C, G, T nucleotides. "
            f"Check TRF output for data corruption."
        )
    """
    pass


class FileFormatError(SatellomeError):
    """
    Raised for file format and parsing issues.

    This includes:
    - Malformed file formats (FASTA, GFF, TRF, etc.)
    - Missing required fields
    - Invalid file modes
    - Parsing errors

    Example:
        raise FileFormatError(
            f"Invalid file mode '{mode}': expected 'r', 'w', or 'a'. "
            f"Check file opening mode specification."
        )
    """
    pass


class ConfigurationError(SatellomeError):
    """
    Raised for configuration and setup issues.

    This includes:
    - Non-callable objects where functions expected
    - Invalid parameter values
    - Missing required configuration
    - Incompatible settings

    Example:
        raise ConfigurationError(
            f"format_func must be callable, got {type(obj).__name__}. "
            f"Provide a function that takes a model object and returns a formatted string."
        )
    """
    pass


# Import ValidationError from validation module to make it available
# as part of the exception hierarchy
try:
    from satellome.core_functions.tools.validation import (
        ValidationError,
        FastaValidationError,
        GFFValidationError,
        BinaryValidationError,
        OutputDirValidationError,
    )
except ImportError:
    # If validation module doesn't exist yet, define a basic ValidationError
    class ValidationError(SatellomeError):
        """
        Raised for input validation failures.

        This includes:
        - Invalid FASTA files
        - Invalid GFF files
        - Missing binaries
        - Invalid output directories
        """
        pass

    FastaValidationError = ValidationError
    GFFValidationError = ValidationError
    BinaryValidationError = ValidationError
    OutputDirValidationError = ValidationError


__all__ = [
    'SatellomeError',
    'StatisticsError',
    'SequenceError',
    'FileFormatError',
    'ConfigurationError',
    'ValidationError',
    'FastaValidationError',
    'GFFValidationError',
    'BinaryValidationError',
    'OutputDirValidationError',
]
