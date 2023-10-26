#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 29.09.2022
# @author: Marina Popova
# @contact: marinaalexpopova@yandex.ru


from collections import defaultdict

import click


@click.command()
@click.option("--file_path", type=click.Path(), required=True)
@click.option(
    "--out_path", type=click.Path(), help="output file name", default="clusters.txt"
)
def repeat_clusters(file_path, out_path):
    with open(file_path) as file:
        line = file.readline()
        repeats = defaultdict(int)
        repeat = [0, "Unknown", 0]
        total = 0
        while True:
            if not line:
                break
            if line.startswith("#=GF ID"):
                repeat[0] = line.split("    ")[1]
            if line.startswith("#=GF TP"):
                repeat[1] = line.split("    ")[1]
            if line.startswith("#=GF SQ"):
                repeat[2] = int(line.split("    ")[1])
                total += int(line.split("    ")[1])
            if repeat[2]:
                repeats[repeat[1]] += repeat[2]
                repeat = [0, "Unknown", 0]
            line = file.readline()
    print(total)
    with open(out_path, "w") as file:
        for key in repeats.keys():
            file.write(f"{key}: {repeats[key]}\n")


if __name__ == "__main__":
    repeat_clusters()
