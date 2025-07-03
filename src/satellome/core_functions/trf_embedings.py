#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.02.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import math
from tqdm import tqdm
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity

from satellome.core_functions.tools.processing import get_revcomp


def get_pentatokens():
    token2id = {}
    token2revtoken = {}
    i = 0
    for n1 in "ACGT":
        for n2 in "ACGT":
            for n3 in "ACGT":
                for n4 in "ACGT":
                    for n5 in "ACGT":
                        token = "".join([n1, n2, n3, n4, n5])
                        token2id[token] = i
                        i += 1
                        token2revtoken[token] = get_revcomp(token)
    return token2id, token2revtoken


def create_vector(token2id, token2revtoken, seq, k=5):
    seq = seq.upper()
    vector = np.zeros((1, len(token2id)))
    N = len(seq) - k + 1
    for i in range(N):
        token = seq[i : i + k]
        if "N" in token:
            continue
        vector[0, token2id[token]] += 1
        vector[0, token2id[token2revtoken[token]]] += 1
    vector /= 2 * N
    return vector


def fill_vectors(df_trs, token2id, token2revtoken, k=5):
    tr2vector = {}
    for id1, x in df_trs.iterrows():
        tr2vector[id1] = create_vector(token2id, token2revtoken, x["seq"], k)
    return tr2vector


def fill_vectors_arrays(arrays, token2id, token2revtoken, k=5):
    tr2vector = {}
    for iid, array in enumerate(arrays):
        vector = np.zeros((1, len(token2id)))
        seq = array.upper()
        N = len(seq) - k + 1
        for i in range(N):
            token = seq[i : i + k]
            if "N" in token:
                continue
            vector[0, token2id[token]] += 1
            vector[0, token2id[token2revtoken[token]]] += 1
        vector /= 2 * N
        tr2vector[iid] = vector
    return tr2vector


def compute_distances(tr2vector):
    return compute_distances_cosine(tr2vector)


def compute_distances_cosine(tr2vector):
    distances = {}
    keys = list(tr2vector.keys())
    for i, id1 in tqdm(enumerate(keys), total=len(keys), desc="Compute distances"):
        for id2 in keys[i:]:
            distances[(id1, id2)] = get_cosine_distance(tr2vector[id1], tr2vector[id2])
            distances[(id2, id1)] = distances[(id1, id2)]
    return distances


def get_cosine_distance(vector1, vector2):
    return 100 * (1 - cosine_similarity(vector1, vector2)[0][0])


def compute_distances_euclidean(tr2vector):
    distances = {}
    keys = list(tr2vector.keys())
    for i, id1 in tqdm(enumerate(keys), total=len(keys), desc="Compute distances"):
        for id2 in keys[i:]:
            distances[(id1, id2)] = 100 * math.dis(tr2vector[id1], tr2vector[id2])[0][0]
            distances[(id2, id1)] = distances[(id1, id2)]
    return distances


def get_disances(df_trs):
    print("Computing distances...")
    token2id, token2revtoken = get_pentatokens()
    tr2vector = fill_vectors(df_trs, token2id, token2revtoken, k=5)
    distances = compute_distances(tr2vector)
    return distances, tr2vector


token2id, token2revtoken = get_pentatokens()
