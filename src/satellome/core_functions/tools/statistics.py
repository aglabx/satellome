#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 26.09.2010
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com
"""
Basic statistical functions for genomic data analysis.

Provides simple descriptive statistics (mean, variance, standard deviation),
frequency counting, and t-test computations. Designed for analyzing numerical
properties of tandem repeats and genomic features.

Functions:
    get_mean: Calculate arithmetic mean
    get_variance: Calculate population variance
    get_sigma: Calculate sum of absolute deviations from mean
    get_standard_deviation: Calculate standard deviation from variance
    t_test: Perform Student's t-test
    get_element_frequencies: Count element occurrences
    get_simple_statistics: Compute all basic statistics at once

Key Features:
    - Input validation with informative error messages
    - Efficient frequency counting with defaultdict
    - Comprehensive statistics dict output
    - Handles empty data gracefully (where appropriate)

Example:
    >>> # Basic statistics
    >>> data = [100, 150, 200, 180, 120]
    >>> stats = get_simple_statistics(data)
    >>> print(f"Mean: {stats['mean']:.2f}")
    >>> print(f"SD: {stats['standard_deviation']:.2f}")
    >>>
    >>> # Frequency counting
    >>> categories = ['A', 'B', 'A', 'C', 'B', 'A']
    >>> freqs = get_element_frequencies(categories)
    >>> print(dict(freqs))  # {'A': 3, 'B': 2, 'C': 1}
    >>>
    >>> # T-test
    >>> sample_mean = 150.0
    >>> dist_mean = 120.0
    >>> variance = 625.0
    >>> N = 30
    >>> t = t_test(sample_mean, dist_mean, variance, N)

See Also:
    satellome.core_functions.exceptions: StatisticsError for validation failures
"""

import math
from collections import defaultdict
from satellome.core_functions.exceptions import StatisticsError


def get_variance(data):
    """
    Calculate population variance for numeric data.

    Args:
        data (list): List of numeric values

    Returns:
        float: Population variance (mean of squared deviations)

    Raises:
        StatisticsError: If data is empty

    Example:
        >>> get_variance([100, 150, 200])
        1666.67
    """
    if not data:
        raise StatisticsError(
            "Cannot compute variance: empty data array. "
            "Ensure the input contains at least one valid data point."
        )
    mean = get_mean(data)
    N = float(len(data))
    return sum([(x - mean) ** 2 for x in data]) / N


def get_sigma(data):
    """
    Calculate sum of absolute deviations from mean.

    Args:
        data (list): List of numeric values

    Returns:
        float: Sum of |xi - mean| for all values (0 if single element)

    Raises:
        StatisticsError: If data is empty

    Example:
        >>> get_sigma([100, 150, 200])
        100.0
    """
    if not data:
        raise StatisticsError(
            "Cannot compute sigma: empty data array. "
            "Ensure the input contains at least one valid data point."
        )
    n = len(data)
    if n == 1:
        return 0
    mean = get_mean(data)
    return sum([abs(x - mean) for x in data])


def get_mean(data):
    """
    Calculate arithmetic mean of numeric data.

    Args:
        data (list): List of numeric values

    Returns:
        float: Arithmetic mean (sum / count)

    Raises:
        StatisticsError: If data is empty

    Example:
        >>> get_mean([100, 150, 200])
        150.0
    """
    if not data:
        raise StatisticsError(
            "Cannot compute mean: empty data array. "
            "Ensure the input contains at least one valid data point."
        )
    sum_x = sum(data)
    mean = float(sum_x) / len(data)
    return mean


def get_standard_deviation(variance):
    """
    Calculate standard deviation from variance.

    Args:
        variance (float): Population or sample variance (must be non-negative)

    Returns:
        float: Standard deviation (square root of variance)

    Raises:
        StatisticsError: If variance is negative

    Example:
        >>> get_standard_deviation(625.0)
        25.0
    """
    if variance < 0:
        raise StatisticsError(
            f"Invalid variance value: {variance}. "
            f"Variance must be non-negative. "
            f"Check input data for outliers or computation errors."
        )
    return math.sqrt(variance)


def t_test(sample_mean, dist_mean, variance, N):
    """
    Perform Student's t-test for mean comparison.

    Args:
        sample_mean (float): Observed sample mean
        dist_mean (float): Expected/population mean
        variance (float): Sample variance (must be positive)
        N (int): Sample size (must be positive)

    Returns:
        float: T-statistic value

    Raises:
        StatisticsError: If N <= 0 or variance <= 0

    Example:
        >>> t_test(sample_mean=150, dist_mean=120, variance=625, N=30)
        6.57...

    Note:
        Formula: t = (sample_mean - dist_mean) / sqrt(variance / N)
    """
    if N <= 0:
        raise StatisticsError(
            f"Invalid sample size: N={N}. "
            f"Sample size must be positive (N > 0). "
            f"Check that the data array is not empty and N is correctly computed."
        )
    if variance <= 0:
        raise StatisticsError(
            f"Invalid variance value: {variance}. "
            f"Variance must be positive for t-test computation. "
            f"Check input data variance calculation or use alternative statistical test."
        )
    return (sample_mean - dist_mean) / float(math.sqrt(variance / N))


def get_element_frequencies(data):
    """
    Count occurrences of each element in data.

    Args:
        data (iterable): Sequence of hashable elements

    Returns:
        defaultdict(int): Maps element to occurrence count

    Example:
        >>> freqs = get_element_frequencies(['A', 'B', 'A', 'C', 'B', 'A'])
        >>> dict(freqs)
        {'A': 3, 'B': 2, 'C': 1}
    """
    d = defaultdict(int)
    for element in data:
        d[element] += 1
    return d


def get_simple_statistics(data):
    """
    Compute all basic statistics in one call.

    Convenience function that calculates mean, variance, sigma,
    and standard deviation. Returns zeros for empty data.

    Args:
        data (list): List of numeric values (empty list returns zeros)

    Returns:
        dict: Statistics dictionary with keys:
            - 'mean' (float): Arithmetic mean
            - 'variance' (float): Population variance
            - 'sigma' (float): Sum of absolute deviations
            - 'standard_deviation' (float): Standard deviation

    Example:
        >>> data = [100, 150, 200, 180, 120]
        >>> stats = get_simple_statistics(data)
        >>> print(f"Mean: {stats['mean']:.1f}")
        Mean: 150.0
        >>> print(f"SD: {stats['standard_deviation']:.1f}")
        SD: 35.8

    Note:
        - Empty data returns all zeros (no exception raised)
        - Uses population variance (divide by N, not N-1)
    """
    if not data:
        result = {
            "mean": 0,
            "variance": 0,
            "sigma": 0,
            "standard_deviation": 0,
        }
        return result
    variance = get_variance(data)
    mean = get_mean(data)
    result = {
        "mean": mean,
        "variance": variance,
        "sigma": get_sigma(data),
        "standard_deviation": get_standard_deviation(variance),
    }
    return result
