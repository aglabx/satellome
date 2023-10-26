#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 31.08.2022
# @author: Marina Popova
# @contact: marinaalexpopova@yandex.ru

import click
import matplotlib
import seaborn as sns


def read_fasta(fasta):
    length = []
    gc_content = []
    with open(fasta) as file:
        line = file.readline().strip()
        sequence = []
        while True:
            if not line:
                sequence_length = len("".join(sequence))
                length.append(sequence_length)
                gc = "".join(sequence).lower().count("g") + "".join(
                    sequence
                ).lower().count("c")
                gc_content.append(round(gc / sequence_length * 100, 2))
                break
            if line.startswith(">"):
                if sequence:
                    sequence_length = len("".join(sequence))
                    length.append(sequence_length)
                    gc = "".join(sequence).lower().count("g") + "".join(
                        sequence
                    ).lower().count("c")
                    gc_content.append(round(gc / sequence_length * 100, 2))
                    sequence = []
                name = line[1:]
            else:
                sequence.append(line.lower())
            line = file.readline().strip()
    return length, gc_content


@click.command()
@click.option("--file_path", type=click.Path(), required=True)
@click.option("--gc_path", type=click.Path(), default="gc.png")
@click.option("--length_path", type=click.Path(), default="distribution.png")
def distribution(file_path, gc_path, length_path):
    matplotlib.use("Agg")
    length, gc_content = read_fasta(file_path)
    fig = sns.displot(length)
    fig.set(xlabel="length of tandem repeat")
    fig.savefig(length_path)

    fig = sns.displot(gc_content)
    fig.set(xlabel="% GC")
    fig.savefig(gc_path)


if __name__ == "__main__":
    distribution()
