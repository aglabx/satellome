#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 13.10.2022
# @author: Marina Popova
# @contact: marinaalexpopova@yandex.ru


import pathlib

import click
import pandas as pd
import yaml


def read_yaml(file_name):
    with open(file_name) as f:
        yaml_file = yaml.load(f, Loader=yaml.FullLoader)
    column_names = ["10kb", "3kb", "1kb", "compex", "fSSR", "tSSR", "micro", "pmicro"]
    row = [yaml_file["species"]]
    row.append(yaml_file["work_files"]["assembly_stats"]["dataset"]["total_length"])
    for i in column_names:
        if yaml_file["work_files"]["repeats"]["dataset"]["trevis"][i]:
            row.append(yaml_file["work_files"]["repeats"]["dataset"]["trevis"][i]["n"])
            row.append(
                yaml_file["work_files"]["repeats"]["dataset"]["trevis"][i]["pgenome"]
            )
        else:
            continue
    return row


def make_table(directory_path):
    table = []
    for file_name in pathlib.Path(directory_path).glob("*.yaml"):
        table.append(read_yaml(file_name))
    return table


@click.command()
@click.option("--directory_path", type=click.Path(), required=True)
@click.option(
    "--output_file",
    type=click.Path(),
    help="output file name",
    default="repeat_stat.tsv",
)
def make_csv(directory_path, output_file):
    table = make_table(directory_path)
    df = pd.DataFrame(
        table,
        columns=[
            "Species",
            "total length of assembly",
            "# of 10kb repeat",
            "% in genome of 10kb repeat",
            "# of 3kb repeat",
            "% in genome of 3kb repeat",
            "# of 1kb repeat",
            "% in genome of 1kb repeat",
            "# of complex",
            "% of complex",
            "# of fSSR",
            "% of fSSR",
            "# of tSSR",
            "% of tSSR",
            "# of micro",
            "% of micro",
            "# of pmicro",
            "% of pmicro",
        ],
    )
    df.to_csv(output_file, sep="\t")


if __name__ == "__main__":
    make_csv()
