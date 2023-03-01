#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 31.08.2022
# @author: Marina Popova
# @contact: marinaalexpopova@yandex.ru


from itertools import product
from sklearn.cluster import AgglomerativeClustering
from collections import defaultdict
from sklearn.metrics import pairwise
from datetime import datetime
import re
import numpy as np
import click


def read_fasta(fasta):
    sequences = {}
    with open(fasta) as file:
        line = file.readline().strip()
        sequence = []
        while True:
            if not line:
                sequences[name] = "".join(sequence)
                break
            if line.startswith(">"):
                if sequence:
                    sequences[name] = "".join(sequence)
                sequence = []
                name = line[1:]
            else:
                sequence.append(line.lower())
            line = file.readline().strip()
    return sequences


def make_matrix(k, sequences, sequences_names):
    k_mers = ["".join(list(i)) for i in product("actg", repeat=k)]
    matrix = np.zeros((len(sequences_names), len(k_mers)))
    for i in range(len(sequences_names)):
        for j in range(len(k_mers)):
            if k_mers[j] in sequences[sequences_names[i]]:
                matrix[i][j] = len(
                    re.findall(rf"(?={k_mers[j]})", sequences[sequences_names[i]])
                )
    distance_matrix = pairwise.cosine_distances(matrix)
    return distance_matrix


def final_matrix(fasta):
    sequences = read_fasta(fasta)
    sequences_names = [i for i in sequences.keys()]
    final_matrix = np.zeros((len(sequences_names), len(sequences_names)))
    for i in range(2, 6):
        distance_matrix = make_matrix(i, sequences, sequences_names)
        final_matrix = final_matrix + distance_matrix
    final_matrix = final_matrix / 4
    return final_matrix, sequences_names


def write_output(out, clusters_dict, time, mode_value):
    with open(out, "w") as file:
        file.write(f"Run time: {time}\n")
        file.write(f"Parametrs: threshold - {mode_value}\n")
        file.write("Clusters:\n")
        for i in sorted(clusters_dict.keys()):
            file.write(f"Cluster #{i}: {','.join(clusters_dict[i])}\n")


@click.command()
@click.option("--fasta", type=click.Path(), required=True, help="fasta file")
@click.option(
    "--mode",
    type=click.Choice(["clust", "threshold"]),
    default="clust",
    help="clust mode if you know how many clusters should be or threshold mode if you do not know how many clusters should be",
)
@click.option(
    "--mode_value",
    default=100.0,
    help="n clusters for clust mode/ float for threshold mode: the linkage distance threshold above which, clusters will not be merged",
)
@click.option(
    "--out", type=click.Path(), help="output file name", default="clusters.txt"
)
def cluster(fasta, mode, mode_value, out):
    start = datetime.now()
    matrix, sequences_names = final_matrix(fasta)
    if mode == "clust":
        clustering = AgglomerativeClustering(
            n_clusters=int(mode_value), affinity="precomputed", linkage="average"
        ).fit(matrix)
    else:
        clustering = AgglomerativeClustering(
            distance_threshold=float(mode_value),
            affinity="precomputed",
            linkage="average",
            n_clusters=None,
        ).fit(matrix)
    clusters_dict = defaultdict(list)
    for i in range(len(sequences_names)):
        clusters_dict[clustering.labels_[i]].append(sequences_names[i])
    finish = datetime.now()
    time = finish - start
    print(f"spent time {time}")
    write_output(out, clusters_dict, time, mode_value)


if __name__ == "__main__":
    cluster()
