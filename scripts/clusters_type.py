#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 04.09.2022
# @author: Marina Popova
# @contact: marinaalexpopova@yandex.ru

from collections import defaultdict
import click

# import yaml


def checking(species, cluster):
    cluster_species = [name.split("_")[:-1] for name in cluster]
    cluster_species = set([" ".join(name).lower() for name in cluster_species])
    intersection_species = cluster_species.intersection(species)
    return intersection_species


@click.command()
@click.option("--file_path", type=click.Path(), required=True)
@click.option("--species", type=click.Path(), required=True)
@click.option("--tandems_path", type=click.Path(), required=True)
@click.option("--all", type=click.Path(), default="clusters_with_all_species.txt")
@click.option("--unique", type=click.Path(), default="clusters_with_unique_species.txt")
def clusters(file_path, species, tandems_path, all, unique):
    with open(species) as file:
        species = set()
        line = file.readline().strip()
        while True:
            if not line:
                break
            species.add(line.lower().replace("_", " "))
            line = file.readline().strip()
    tandems = {}
    with open(tandems_path) as file:
        line = file.readline().strip()
        sequence = []
        while True:
            if not line:
                tandems[name] = "".join(sequence)
                break
            if line.startswith(">"):
                if sequence:
                    tandems[name] = "".join(sequence)
                sequence = []
                name = line[1:]
            else:
                sequence.append(line.lower())
            line = file.readline().strip()
    with open(file_path) as file:
        clusters = {}
        line = file.readline().strip()
        while True:
            if not line:
                break
            if line.startswith("Cluster #"):
                name = line.split(":")[0].split("#")[1]
                sequences = line.split(": ")[1].split(",")
                clusters[name] = sequences
            line = file.readline().strip()

    multi_species_clusters = defaultdict()
    for cluster_number in clusters.keys():
        intersection = checking(species, clusters[cluster_number])
        if len(intersection) / len(species) > 0.8:
            for seq in clusters[cluster_number]:
                if frozenset(intersection) not in multi_species_clusters.keys():
                    multi_species_clusters[frozenset(intersection)] = defaultdict(list)
                multi_species_clusters[frozenset(intersection)][cluster_number].append(
                    [seq, tandems[seq]]
                )
    uniq_species_clusters = defaultdict()
    for cluster_number in clusters.keys():
        intersection = checking(species, clusters[cluster_number])
        if len(intersection) == 1:
            for seq in clusters[cluster_number]:
                if frozenset(intersection) not in uniq_species_clusters.keys():
                    uniq_species_clusters[frozenset(intersection)] = defaultdict(list)
                uniq_species_clusters[frozenset(intersection)][cluster_number].append(
                    [seq, tandems[seq]]
                )

    with open(all, "w") as file:
        file.write("Clusters with multi named species\n")
        for species in multi_species_clusters.keys():
            file.write(f"Clusters with {species}\n")
            for cluster, seqs in multi_species_clusters[species].items():
                file.write(f"\tCluster #{cluster}\n")
                tandem_seqs = "\n".join("\t".join(str(i) for i in seq) for seq in seqs)
                file.write(f"\t\t{tandem_seqs}\n")

    with open(unique, "w") as file:
        file.write("Clusters with unique species\n")
        for species in uniq_species_clusters.keys():
            file.write(f"Clusters with {species}\n")
            for cluster, seqs in uniq_species_clusters[species].items():
                file.write(f"\tCluster #{cluster}\n")
                tandem_seqs = "\n".join("\t".join(str(i) for i in seq) for seq in seqs)
                file.write(f"\t\t{tandem_seqs}\n")


if __name__ == "__main__":
    clusters()
