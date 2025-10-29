#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Constants and configuration parameters for Satellome pipeline.

This module contains all magic numbers and configuration constants
used throughout the Satellome pipeline to improve maintainability
and code readability.
"""

# ============================================================================
# TRF (Tandem Repeat Finder) Parameters
# ============================================================================

# TRF algorithm parameters
# Format: match, mismatch, indel, match_prob, indel_prob, min_score, max_period
TRF_MATCH_WEIGHT = 2          # Matching weight
TRF_MISMATCH_WEIGHT = 5       # Mismatching penalty
TRF_INDEL_WEIGHT = 7          # Insertion/deletion penalty
TRF_MATCH_PROBABILITY = 80    # Matching probability (%)
TRF_INDEL_PROBABILITY = 10    # Insertion/deletion probability (%)
TRF_MIN_SCORE = 50            # Minimum alignment score
TRF_MAX_PERIOD = 2000         # Maximum period size
TRF_MIN_LENGTH = 200          # Minimum tandem repeat length (-l parameter)

# Default TRF command parameters as a list
TRF_DEFAULT_PARAMS = [
    str(TRF_MATCH_WEIGHT),
    str(TRF_MISMATCH_WEIGHT),
    str(TRF_INDEL_WEIGHT),
    str(TRF_MATCH_PROBABILITY),
    str(TRF_INDEL_PROBABILITY),
    str(TRF_MIN_SCORE),
    str(TRF_MAX_PERIOD)
]

# TRF command flags
TRF_FLAGS = ["-l", str(TRF_MIN_LENGTH), "-d", "-h"]  # -d: data file, -h: suppress HTML

# ============================================================================
# Scaffold and Contig Size Thresholds
# ============================================================================

# Minimum scaffold sizes for different purposes
MIN_SCAFFOLD_LENGTH_DEFAULT = 10000      # Default minimum scaffold length (10 kb)
MIN_SCAFFOLD_LENGTH_LARGE = 1000000      # Large scaffold threshold (1 Mb)
MIN_SCAFFOLD_LENGTH_FILTER = 100000      # Filter threshold for drawing (100 kb)

# Large tandem repeat thresholds
TR_CUTOFF_DEFAULT = 1000                 # Default cutoff for large TRs (1 kb)
TR_CUTOFF_MEDIUM = 3000                  # Medium TR threshold (3 kb)
TR_CUTOFF_LARGE = 10000                  # Large TR threshold (10 kb)

# ============================================================================
# K-mer Filtering Parameters
# ============================================================================

KMER_THRESHOLD_DEFAULT = 90000           # Default unique k-mer threshold for repeat detection
KMER_SIZE = 21                            # Default k-mer size for analysis
MIN_FREE_DISK_GB = 500                    # Minimum free disk space required in GB

# ============================================================================
# Visualization Parameters
# ============================================================================

# Canvas dimensions
CANVAS_WIDTH_DEFAULT = 800               # Default canvas width for plots
CANVAS_HEIGHT_DEFAULT = 800              # Default canvas height for plots
CANVAS_WIDTH_KARYOTYPE = 1200           # Karyotype plot width
CANVAS_HEIGHT_MIN = 400                  # Minimum canvas height
CANVAS_HEIGHT_MAX = 8000                 # Maximum canvas height

# Chromosome visualization
CHROMOSOME_HEIGHT = 50                   # Height per chromosome in pixels
VERTICAL_SPACER = 20                     # Vertical spacing between chromosomes
BASE_HEIGHT = 150                        # Base height for margins, title, etc.
MARGIN_TOP = 100                         # Top margin
MARGIN_BOTTOM = 50                      # Bottom margin
MARGIN_LEFT = 100                       # Left margin
MARGIN_RIGHT = 50                       # Right margin

# Enhancement sizes for visualization
ENHANCE_DEFAULT = 1000000                # Default enhancement size (1 Mb)
ENHANCE_LARGE = 2000000                  # Large enhancement size (2 Mb)
ENHANCE_SMALL = 100000                   # Small enhancement size (100 kb)
DRAWING_ENHANCING_DEFAULT = 100000       # Default drawing enhancement

# Gap visualization
GAP_CUTOFF_DEFAULT = 1000                # Default gap cutoff (1 kb)
GAP_SEARCH_WINDOW = 10000               # Window for searching gaps around TRs

# ============================================================================
# Data Processing Thresholds
# ============================================================================

# Clustering and filtering
MAX_ITEMS_FOR_CLUSTERING = 2000         # Maximum items before sampling for clustering
SAMPLE_SIZE_FOR_CLUSTERING = 1999       # Sample size when too many items
MIN_CLUSTER_SIZE = 3                    # Minimum cluster size
START_CUTOFF_MAX = 50                   # Maximum start cutoff for distance calculation

# Batch processing
BATCH_SIZE_DEFAULT = 10000              # Default batch size for processing
PROGRESS_UPDATE_INTERVAL = 10000        # Update progress every N items

# ============================================================================
# Classification Categories
# ============================================================================

# Size ranges for tandem repeat classification
TR_SIZE_RANGES = [
    (1, 10, "Micro (1-10 bp)"),
    (11, 100, "Small (11-100 bp)"),
    (101, 1000, "Medium (101-1000 bp)"),
    (1001, 10000, "Medium-Large (1-10 kb)"),
    (10001, float('inf'), "Large (>10 kb)")
]

# Gap size ranges
GAP_SIZE_RANGES = [
    (1, 100, "Small gaps (1-100 bp)"),
    (101, 1000, "Medium gaps (101-1000 bp)"),
    (1001, 10000, "Large gaps (1-10 kb)"),
    (10001, 100000, "Very large gaps (10-100 kb)"),
    (100001, float('inf'), "Huge gaps (>100 kb)")
]

# ============================================================================
# File System and I/O
# ============================================================================

MAX_OUTPUT_CHARS = 30000                # Maximum characters in output before truncation
DEFAULT_TIMEOUT_MS = 120000             # Default command timeout (2 minutes)
MAX_TIMEOUT_MS = 600000                 # Maximum command timeout (10 minutes)

# ============================================================================
# Display and Formatting
# ============================================================================

SEPARATOR_LINE = "-" * 50               # Standard separator line for console output
SEPARATOR_LINE_DOUBLE = "=" * 50        # Double separator line for important sections

# ============================================================================
# Colors for Visualization (Hex colors)
# ============================================================================

COLORS = {
    "gaps": "rgba(0, 0, 0)",             # Black for gaps
    "aN_repeat": "#FF00FF",               # Magenta for aN repeats
    "aN_gap": "#663399",                  # Rebecca purple for aN gaps
    "Na_repeat": "#00CED1",               # Dark turquoise for Na repeats
    "Na_gap": "#00BFFF",                  # Deep sky blue for Na gaps
    "aNa_repeat": "#00FF7F",              # Spring green for aNa repeats
    "aNa_gap": "#228B22",                 # Forest green for aNa gaps
}

# ============================================================================
# Default Values
# ============================================================================

DEFAULT_PROJECT_NAME = "satellome"
DEFAULT_TAXON_NAME = "Unknown"
DEFAULT_THREADS = 10
DEFAULT_PARALLEL_JOBS = 10

RECURSION_LIMIT_DEFAULT = 20000000      # Default recursion limit for deep recursions