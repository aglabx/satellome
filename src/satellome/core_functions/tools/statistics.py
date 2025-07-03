#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 26.09.2010
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com
"""
Collection of functions related to simple statistics.
"""

import math
from collections import defaultdict


def get_variance(data):
    """
    Calculated variance for given list.
    @param data: list of numbers
    @return: variance
    """
    if not data:
        raise Exception("Empty arrays for variance computation.")
    mean = get_mean(data)
    N = float(len(data))
    return sum([(x - mean) ** 2 for x in data]) / N


def get_sigma(data):
    """
    Calculate sigma, return sum(module(xi - mean).
    @param data: list of numbers
    @return: sigma
    """
    if not data:
        raise Exception("Empty arrays for variance computation.")
    n = len(data)
    if n == 1:
        return 0
    mean = get_mean(data)
    return sum([abs(x - mean) for x in data])


def get_mean(data):
    """
    Calculated mean for given list.
    @param data: list of numbers
    @return: arithmetic mean
    """
    if not data:
        raise Exception("Empty data.")
    sum_x = sum(data)
    mean = float(sum_x) / len(data)
    return mean


def get_standard_deviation(variance):
    """
    Get sample deviation for variance.
    @param variance:
    @return: standard deviation
    """
    if variance < 0:
        raise Exception("Wrong variance value %s" % variance)
    return math.sqrt(variance)


def t_test(sample_mean, dist_mean, variance, N):
    """
    Compute t-test.
    @param sample_mean: sample mean
    @param dist_mean: distribution mean
    @param variance: sample variance
    @param N: sample size
    @return: t_test
    """
    if N <= 0:
        raise Exception("Wrong N value %s" % N)
    if variance <= 0:
        raise Exception("Wrong variance value %s" % variance)
    return (sample_mean - dist_mean) / float(math.sqrt(variance / N))


def get_element_frequencies(data):
    """
    Get default dictionary of element frequencies in given list
    @param data:
    @return: default dictionary element to tf
    """
    d = defaultdict(int)
    for element in data:
        d[element] += 1
    return d


def get_simple_statistics(data):
    """
    Compute simple statistics
    @param data: list of numbers
    @return: dictionary with keys (mean, variance, sigma, standard_deviation)
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
