#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.02.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import editdistance as ed
from suffix_trees import STree
from tqdm import tqdm

from satelome.core_functions.tools.processing import get_revcomp


def hamming_sliding_distance(seq1, seq2, lcs_cutoff=0.2):
    """Compute Hamming Sliding distance between two sequences."""
    st = STree.STree([seq1, seq2])
    suff = st.lcs()
    if len(suff) / len(seq1) < lcs_cutoff:
        return len(seq1)
    pos1 = seq1.find(suff)
    pos2 = seq2.find(suff)
    if pos1 > 0:
        seq1 = seq1[pos1:] + seq1[:pos1]
    if pos2 > 0:
        seq2 = seq2[pos2:] + seq2[:pos2]
    hd = sum([seq1[j] != seq2[j] for j in range(len(seq2))])
    return hd


def compute_hs_distances(sequences, seq2id, distance_cutoff=0.1, lcs_cutoff=0.2):
    """Compute Hamming Sliding distance between all sequences."""
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

            d = hamming_sliding_distance(consensus1, consensus2, lcs_cutoff=lcs_cutoff)
            d /= len(consensus1)
            if d <= distance_cutoff:
                sh_distances[key] = d
                sh_distances[(seq2id[consensus2], seq2id[ori_consensus1])] = d
            else:
                consensus1 = get_revcomp(consensus1)
                d = hamming_sliding_distance(
                    consensus1, consensus2, lcs_cutoff=lcs_cutoff
                )
                d /= len(consensus1)
                if d <= distance_cutoff:
                    sh_distances[key] = d
                    sh_distances[(seq2id[consensus2], seq2id[ori_consensus1])] = d
    return sh_distances


def compute_edit_distances(sequences, seq2id, distance_cutoff=0.1):
    """Compute edit distance between all sequences."""
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
