#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 31.08.2022
# @author: Marina Popova
# @contact: marinaalexpopova@yandex.ru

import click


@click.command()
@click.option("--file_path", type=click.Path(), required=True)
@click.option("--out_path", type=click.Path(), default="clusters.txt")
def clusters(file_path, out_path):
    with open(file_path) as file:
        clusters = []
        line = file.readline().strip()
        while True:
            if not line:
                break
            if line.startswith("Cluster #"):
                name = line.split(":")[0]
                sequences = [
                    " ".join(name.split("_")[:-1])
                    for name in line.split(":")[1][1:].split(",")
                ]
                if set(sequences) not in clusters:
                    clusters.append(set(sequences))
            line = file.readline().strip()

    with open(out_path, "w") as file:
        for i in range(len(clusters)):
            file.write(f"{i}: {clusters[i]}\n")


if __name__ == "__main__":
    clusters()
