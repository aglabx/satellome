#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 19.08.2022
# @author: Marina Popova
# @contact: marinaalexpopova@yandex.ru


import click
import re


@click.command()
@click.argument("file_path", type=click.Path())
@click.argument("out_path", type=click.Path())
def clusters(file_path, out_path):
    with open(file_path) as file:
        line = file.readline().strip()
        representative = re.split(r",|\t|;", line)[0]
        clusters = {}
        cluster = set()
        i = 1
        while True:
            if not line:
                if cluster not in clusters.values():
                    clusters[i] = cluster
                    i += 1
                break
            if representative == re.split(r",|\t|;", line)[0]:
                cluster.add(representative)
                cluster.add(re.split(r",|\t|;", line)[1])
            else:
                if cluster not in clusters.values():
                    clusters[i] = cluster
                    i += 1
                representative = re.split(r",|\t|;", line)[0]
                cluster = set()
                cluster.add(representative)
                cluster.add(re.split(r",|\t|;", line)[1])
            line = file.readline().strip()

    with open(out_path, "w") as file:
        for k, clu in clusters.items():
            line = f"Cluster # {k}: {','.join([i for i in clu])}\n"
            file.write(line)


if __name__ == "__main__":
    clusters()
