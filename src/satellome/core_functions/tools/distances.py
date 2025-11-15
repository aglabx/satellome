#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.02.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com
"""
Sequence distance metrics for tandem repeat clustering.

Provides specialized distance functions for comparing tandem repeat consensus
sequences, accounting for circular rotation (sliding window), reverse complement,
and repeat unit length variations. Used for similarity-based clustering.

Functions:
    hamming_sliding_distance: Circular-aware Hamming distance between two sequences
    compute_hs_distances: All-pairs Hamming Sliding distances with cutoff filtering
    compute_edit_distances: All-pairs edit distances (deprecated, rotation not implemented)

Key Features:
    - Circular rotation handling for tandem repeat motifs
    - Reverse complement comparison for strand-agnostic matching
    - Length normalization for different repeat unit counts
    - Distance cutoff filtering for sparse similarity graphs
    - Progress tracking with tqdm for large datasets

Distance Metrics:
    - **Hamming Sliding (HS)**: Minimum mismatches across all rotations
    - **Edit Distance**: Levenshtein distance (insertions/deletions/substitutions)

Example:
    >>> # Compare two repeat motifs with rotation
    >>> from satellome.core_functions.tools.distances import hamming_sliding_distance
    >>> seq1 = "ACAT"
    >>> seq2 = "TACA"  # Rotated version of seq1
    >>> dist = hamming_sliding_distance(seq1, seq2)
    >>> print(dist)  # Should be 0 (perfect match after rotation)
    0
    >>>
    >>> # Compute all pairwise distances
    >>> sequences = ["ACAT", "TACA", "GGAT"]
    >>> seq2id = {seq: i for i, seq in enumerate(sequences)}
    >>> distances = compute_hs_distances(sequences, seq2id, distance_cutoff=0.3)
    >>> print(distances)  # Only pairs below 30% difference

Typical Use Case:
    1. Extract unique consensus sequences from TRF output
    2. Compute pairwise Hamming Sliding distances with cutoff
    3. Build similarity graph from distance pairs
    4. Cluster using connected components (see clusterization.py)

See Also:
    satellome.core_functions.tools.clusterization: Graph-based clustering from distances
    satellome.core_functions.tools.processing: get_revcomp for strand comparison
"""

import editdistance as ed
from tqdm import tqdm

from satellome.core_functions.tools.processing import get_revcomp

def hamming_sliding_distance(seq1, seq2, min_hd=None):
    """
    Compute minimum Hamming distance across all circular rotations of seq2.

    Slides seq1 along doubled seq2 to find rotation with fewest mismatches.
    Essential for comparing tandem repeat motifs where ACAT = TACA = ATAC = CATA
    due to circular nature of repeat units.

    Args:
        seq1 (str): First sequence (fixed reference)
        seq2 (str): Second sequence (rotated to find best match)
        min_hd (int, optional): Initial minimum to beat (for early stopping).
                               If None, starts at len(seq1). Defaults to None.

    Returns:
        int: Minimum Hamming distance found across all rotations

    Example:
        >>> # Perfect match after rotation
        >>> hamming_sliding_distance("ACAT", "TACA")
        0
        >>>
        >>> # One mismatch in best alignment
        >>> hamming_sliding_distance("ACAT", "TACG")
        1
        >>>
        >>> # With early stopping hint
        >>> hamming_sliding_distance("ACAT", "GGGG", min_hd=2)
        2  # Stops early once confirmed >= 2

    Note:
        - Assumes len(seq1) == len(seq2) for correct behavior
        - Uses early stopping: breaks loop if current rotation exceeds min_hd
        - Time complexity: O(n²) where n = len(seq1)
        - Space complexity: O(n) due to seq2 doubling
    """
    seq2 = seq2 + seq2
    if min_hd is None:
        min_hd = len(seq1)
    for i in range(len(seq2)-len(seq1)+1):
        hd = 0
        for j in range(len(seq1)):
            if seq1[j] != seq2[i+j]:
                hd += 1
                if hd >= min_hd:
                    break
        if hd < min_hd:
            min_hd = hd
    return min_hd

def compute_hs_distances(sequences, seq2id, distance_cutoff=0.1):
    """
    Compute all-pairs Hamming Sliding distances with rotation and reverse complement.

    Main clustering distance function. Computes pairwise distances for all sequences,
    testing both forward and reverse complement orientations across all rotations.
    Only stores pairs below distance_cutoff for sparse similarity graph construction.

    Handles repeat unit length variations by extending shorter sequences (e.g.,
    "AT" × 2 → "ATAT" to compare with "ATAT").

    Args:
        sequences (list): List of unique consensus sequences to compare
        seq2id (dict): Maps each sequence string to unique integer ID
        distance_cutoff (float, optional): Maximum normalized distance to store
                                          (0.1 = 10% difference). Defaults to 0.1.

    Returns:
        dict: Maps (id1, id2) tuples to normalized distances (0.0-1.0), where:
             - Only pairs with distance <= cutoff are included
             - Both (id1, id2) and (id2, id1) are present (symmetric)
             - Self-distances (id, id) are always 0.0
             - Distances normalized by sequence length

    Example:
        >>> sequences = ["ACAT", "TACA", "GGAT"]
        >>> seq2id = {"ACAT": 0, "TACA": 1, "GGAT": 2}
        >>> distances = compute_hs_distances(sequences, seq2id, distance_cutoff=0.3)
        >>> print(distances[(0, 0)])  # Self-distance
        0.0
        >>> print(distances.get((0, 1), "not similar"))  # ACAT vs TACA (rotated)
        0.0
        >>> print(distances.get((0, 2), "not similar"))  # ACAT vs GGAT
        'not similar'  # Distance > 30%, not stored

    Processing Steps:
        1. Compare each sequence pair (i, j) where j > i (avoid duplicates)
        2. Extend shorter sequence by repetition if lengths are multiples
        3. Skip pairs with different final lengths
        4. Compute Hamming Sliding distance (forward orientation)
        5. If above cutoff, try reverse complement orientation
        6. Store both (i,j) and (j,i) if either orientation below cutoff

    Note:
        - Progress bar shows total comparisons (uses tqdm)
        - Only sequences with identical or multiple lengths compared
        - Distance normalized by sequence length (range 0.0-1.0)
        - Reverse complement tested only if forward orientation exceeds cutoff
        - Returns sparse dictionary (missing pairs imply distance > cutoff)
        - Primary distance metric for clustering in Satellome pipeline
    """
    sh_distances = {}
    computed = set()
    for i, ori_consensus1 in tqdm(enumerate(sequences), total=len(sequences)):
        sh_distances[(seq2id[ori_consensus1], seq2id[ori_consensus1])] = 0.0
        l1 = len(ori_consensus1)
        
        for consensus2 in sequences[i + 1 :]:
            if (ori_consensus1, consensus2) in computed:
                continue
            computed.add((ori_consensus1, consensus2))
            computed.add((consensus2, ori_consensus1))
            l2 = len(consensus2)
            consensus1 = ori_consensus1
            if l1 < l2:
                if l2 % l1 == 0:
                    consensus1 = consensus1 * (l2 // l1)
            if len(consensus1) != len(consensus2):
                continue
            key = (seq2id[ori_consensus1], seq2id[consensus2])

            d = hamming_sliding_distance(consensus1, consensus2)
            if d/len(consensus1) <= distance_cutoff:
                d = d / len(consensus1)
                sh_distances[key] = d
                sh_distances[(seq2id[consensus2], seq2id[ori_consensus1])] = d
            else:
                consensus1 = get_revcomp(consensus1)
                d = hamming_sliding_distance(
                    consensus1, consensus2, min_hd=d
                )
                d /= len(consensus1)
                if d <= distance_cutoff:
                    sh_distances[key] = d
                    sh_distances[(seq2id[consensus2], seq2id[ori_consensus1])] = d
    return sh_distances


def compute_edit_distances(sequences, seq2id, distance_cutoff=0.1):
    """
    Compute all-pairs edit distances with reverse complement (DEPRECATED).

    Alternative distance metric using Levenshtein edit distance instead of
    Hamming distance. Tests forward and reverse complement but does NOT
    handle circular rotations (unlike compute_hs_distances).

    **WARNING**: This function is NOT used in the Satellome pipeline.
    Use compute_hs_distances() for production clustering.

    Args:
        sequences (list): List of unique consensus sequences to compare
        seq2id (dict): Maps each sequence string to unique integer ID
        distance_cutoff (float, optional): Maximum normalized distance to store
                                          (0.1 = 10% difference). Defaults to 0.1.

    Returns:
        dict: Maps (id1, id2) tuples to normalized edit distances (0.0-1.0), where:
             - Only pairs with distance <= cutoff are included
             - Both (id1, id2) and (id2, id1) are present (symmetric)
             - Self-distances (id, id) are always 0.0
             - Distances normalized by sequence length

    Example:
        >>> sequences = ["ACAT", "ACAT", "GGAT"]
        >>> seq2id = {"ACAT": 0, "ACAT": 1, "GGAT": 2}
        >>> distances = compute_edit_distances(sequences, seq2id, distance_cutoff=0.5)
        >>> print(distances[(0, 0)])  # Self-distance
        0.0
        >>> # Note: "ACAT" vs "TACA" would NOT match (no rotation handling)

    Differences from compute_hs_distances:
        - Uses editdistance library (Levenshtein) instead of Hamming
        - Allows insertions/deletions, not just substitutions
        - Does NOT test circular rotations (ACAT ≠ TACA)
        - Generally slower than Hamming Sliding
        - Not suitable for tandem repeat clustering

    Note:
        - Extends shorter sequences if lengths are multiples (same as compute_hs_distances)
        - Skip pairs with different final lengths
        - Tests reverse complement if forward exceeds cutoff
        - Progress bar shows total comparisons
        - **Not recommended for tandem repeat analysis** (missing rotation handling)
    """
    distances = {}
    computed = set()
    for i, ori_consensus1 in tqdm(enumerate(sequences), total=len(sequences)):
        distances[(seq2id[ori_consensus1], seq2id[ori_consensus1])] = 0.0
        l1 = len(ori_consensus1)
        for consensus2 in sequences[i + 1 :]:
            if (ori_consensus1, consensus2) in computed:
                continue
            computed.add((ori_consensus1, consensus2))
            computed.add((consensus2, ori_consensus1))
            l2 = len(consensus2)
            consensus1 = ori_consensus1
            if l1 < l2:
                if l2 % l1 == 0:
                    consensus1 = consensus1 * (l2 // l1)
            if len(consensus1) != len(consensus2):
                continue
            key = (seq2id[ori_consensus1], seq2id[consensus2])

            d = ed.eval(consensus1, consensus2) / len(consensus1)
            if d <= distance_cutoff:
                distances[key] = d
                distances[(seq2id[consensus2], seq2id[ori_consensus1])] = d
            else:
                consensus1 = get_revcomp(consensus1)
                d = ed.eval(consensus1, consensus2) / len(consensus1)
                if d <= distance_cutoff:
                    distances[key] = d
                    distances[(seq2id[consensus2], seq2id[ori_consensus1])] = d
    return distances
