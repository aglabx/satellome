#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.02.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

from collections import Counter
from dataclasses import dataclass

import networkx as nx


@dataclass
class AnnotatedComponent:
    cid: int
    comp: set
    n_uniq_consensus: int
    n_arrays: int
    comp_fraq: float
    taxons: list
    taxons_string: str
    motif_sizes: dict


def make_taxon_flags(taxon_dict, flags, tx):
    for t in tx:
        flags[taxon_dict[t]] = 1
    return flags


def get_flag_string(taxon_dict, flags):
    result = []
    for tx, i in taxon_dict.items():
        flag = flags[i]
        result.append(f"{tx}{flag}")
    return " ".join(result)


def get_connected_components(distances, id2size):
    """Get connected components sorted by real component size."""

    G = nx.Graph()
    for pair in distances.keys():
        G.add_edge(*pair)
    connected_components = list(nx.connected_components(G))

    connected_components.sort(key=lambda x: -sum([id2size[iid] for iid in x]))

    return connected_components


def annotate_components(
    connected_components, id2size, taxon_dict, id2seq, consensuses_taxons
):
    annotated_components = []
    total = 0
    dataset_size = sum(id2size.values())
    for cid, comp in enumerate(connected_components):
        t = sum([id2size[iid] for iid in comp])
        comp_fraq = 100.0 * t / dataset_size
        total += comp_fraq

        flags = [0 for x in range(len(taxon_dict))]

        sizes = Counter()
        for iid in comp:
            s = id2seq[iid]
            tx = consensuses_taxons[s]
            make_taxon_flags(taxon_dict, flags, tx)
            sizes[len(s)] += 1

        sizes = dict(sizes)

        ac = AnnotatedComponent(
            cid,
            comp,
            len(comp),
            t,
            comp_fraq,
            flags,
            get_flag_string(taxon_dict, flags),
            sizes,
        )
        annotated_components.append(ac)

    return annotated_components
