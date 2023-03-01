#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 10.07.2022
# @author: Marina Popova
# @contact: marinaalexpopova@yandex.ru


import click


def read_file(file_path):
    with open(file_path) as file:
        file_names = [line.strip() for line in file]
    return file_names


def make_fasta(file_path, path):
    file_names = read_file(file_path)
    fasta_table = {}
    for file_name in file_names:
        with open(f"{path}/{file_name}.1kb.trf") as file:
            line = file.readline()
            while True:
                if not line:
                    break
                name = "_".join([file_name, line.split("\t")[2]])
                monomer = line.split("\t")[13]
                fasta_table[name] = monomer * 3
                line = file.readline()
    return fasta_table


@click.command()
@click.option("--file_path", type=click.Path())
@click.option("--path", type=click.Path())
def save_fasta(file_path, path):
    fasta_table = make_fasta(file_path, path)
    with open("tandem_1kb.fasta", "w") as file:
        for name, monomer in fasta_table.items():
            file.write(f">{name}\n{monomer}\n")


if __name__ == "__main__":
    save_fasta()
