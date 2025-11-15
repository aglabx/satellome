#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.02.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com
"""
Graph-based clustering utilities for tandem repeat analysis.

Provides network-based clustering using sequence similarity distances,
with annotation by taxon presence and motif size distribution. Uses
NetworkX for connected component detection.

Classes:
    AnnotatedComponent: Dataclass for cluster metadata

Functions:
    get_connected_components: Find clusters from similarity graph
    annotate_components: Add taxon and size annotations to clusters
    make_taxon_flags: Set taxon presence flags
    get_flag_string: Format taxon flags as string

Key Features:
    - NetworkX-based graph clustering
    - Size-weighted component sorting
    - Multi-taxon presence tracking
    - Motif size distribution analysis
    - Fractional abundance calculation

Example:
    >>> # Build similarity graph and cluster
    >>> distances = {(1, 2): 0.05, (2, 3): 0.03}  # ID pairs
    >>> id2size = {1: 100, 2: 150, 3: 200}  # Array lengths
    >>> components = get_connected_components(distances, id2size)
    >>>
    >>> # Annotate with taxon info
    >>> id2seq = {1: "ACT", 2: "ACT", 3: "GAT"}
    >>> taxon_dict = {"human": 0, "mouse": 1}
    >>> consensuses_taxons = {"ACT": ["human"], "GAT": ["mouse"]}
    >>> annotated = annotate_components(
    ...     components, id2size, taxon_dict, id2seq, consensuses_taxons
    ... )

See Also:
    networkx.connected_components: Underlying clustering algorithm
"""

from collections import Counter
from dataclasses import dataclass

import networkx as nx


@dataclass
class AnnotatedComponent:
    """
    Cluster metadata with taxon and size annotations.

    Attributes:
        cid (int): Cluster ID (0-indexed)
        comp (set): Set of sequence IDs in this cluster
        n_uniq_consensus (int): Number of unique consensus sequences
        n_arrays (int): Total array count (sum of id2size for all IDs)
        comp_fraq (float): Cluster fraction of dataset (percentage)
        taxons (list): Taxon presence flags (binary list)
        taxons_string (str): Formatted taxon flags (e.g., "human1 mouse0")
        motif_sizes (dict): Maps motif length to occurrence count
    """
    cid: int
    comp: set
    n_uniq_consensus: int
    n_arrays: int
    comp_fraq: float
    taxons: list
    taxons_string: str
    motif_sizes: dict


def make_taxon_flags(taxon_dict, flags, tx):
    """
    Set presence flags for given taxons.

    Args:
        taxon_dict (dict): Maps taxon name to flag index
        flags (list): Binary flag list to modify in-place
        tx (list): Taxon names to mark as present

    Returns:
        list: Modified flags (same as input, modified in-place)

    Example:
        >>> taxon_dict = {"human": 0, "mouse": 1, "dog": 2}
        >>> flags = [0, 0, 0]
        >>> make_taxon_flags(taxon_dict, flags, ["human", "dog"])
        [1, 0, 1]
    """
    for t in tx:
        flags[taxon_dict[t]] = 1
    return flags


def get_flag_string(taxon_dict, flags):
    """
    Format taxon flags as human-readable string.

    Args:
        taxon_dict (dict): Maps taxon name to flag index
        flags (list): Binary flag list

    Returns:
        str: Space-separated "taxonN" strings (e.g., "human1 mouse0")

    Example:
        >>> taxon_dict = {"human": 0, "mouse": 1}
        >>> flags = [1, 0]
        >>> get_flag_string(taxon_dict, flags)
        'human1 mouse0'
    """
    result = []
    for tx, i in taxon_dict.items():
        flag = flags[i]
        result.append(f"{tx}{flag}")
    return " ".join(result)


def get_connected_components(distances, id2size):
    """
    Find clusters from pairwise similarity graph, sorted by total size.

    Builds undirected graph from distance pairs and extracts connected
    components. Components sorted by sum of array sizes (descending).

    Args:
        distances (dict): Maps (id1, id2) tuple to distance value
        id2size (dict): Maps sequence ID to array size/count

    Returns:
        list: Connected components (each is set of sequence IDs),
             sorted by total size (largest first)

    Example:
        >>> distances = {(1, 2): 0.05, (2, 3): 0.03, (4, 5): 0.02}
        >>> id2size = {1: 100, 2: 150, 3: 200, 4: 50, 5: 60}
        >>> components = get_connected_components(distances, id2size)
        >>> len(components)
        2
        >>> sum(id2size[i] for i in components[0])  # Largest cluster
        450

    Note:
        - Uses NetworkX Graph for efficiency
        - Distance values not used for clustering (presence/absence only)
        - Sorting ensures largest clusters processed first
    """

    G = nx.Graph()
    for pair in distances.keys():
        G.add_edge(*pair)
    connected_components = list(nx.connected_components(G))

    connected_components.sort(key=lambda x: -sum([id2size[iid] for iid in x]))

    return connected_components


def annotate_components(
    connected_components, id2size, taxon_dict, id2seq, consensuses_taxons
):
    """
    Add taxon and motif size annotations to connected components.

    Enriches raw clusters with metadata: taxon presence flags, motif size
    distributions, fractional abundances, and unique consensus counts.
    Each component gets formatted as AnnotatedComponent dataclass.

    Args:
        connected_components (list): List of components (each is set of sequence IDs)
        id2size (dict): Maps sequence ID to array size/count
        taxon_dict (dict): Maps taxon name to flag index (e.g., {"human": 0, "mouse": 1})
        id2seq (dict): Maps sequence ID to consensus sequence string
        consensuses_taxons (dict): Maps consensus sequence to list of taxon names

    Returns:
        list: AnnotatedComponent objects with full metadata, one per input component

    Example:
        >>> components = [{1, 2}, {3, 4, 5}]
        >>> id2size = {1: 100, 2: 150, 3: 200, 4: 80, 5: 120}
        >>> taxon_dict = {"human": 0, "mouse": 1}
        >>> id2seq = {1: "ACT", 2: "ACT", 3: "GAT", 4: "GAT", 5: "GAT"}
        >>> consensuses_taxons = {"ACT": ["human"], "GAT": ["mouse"]}
        >>> annotated = annotate_components(
        ...     components, id2size, taxon_dict, id2seq, consensuses_taxons
        ... )
        >>> print(annotated[0].n_arrays)  # First cluster total size
        250
        >>> print(annotated[0].taxons_string)  # Taxon presence
        'human1 mouse0'
        >>> print(annotated[1].motif_sizes)  # Motif length distribution
        {3: 3}

    Note:
        - comp_fraq is percentage of total dataset (not fraction)
        - Motif sizes count unique sequences, not array counts
        - Taxon flags are cumulative across all sequences in component
        - Components processed in input order (typically largest-first from get_connected_components)
    """
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
