#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 10.07.2022
# @author: Marina Popova
# @contact: marinaalexpopova@yandex.ru


import click
import re
from Bio import Entrez


def get_id_list(txid: str) -> list:
    Entrez.email = "user@gmail.com"
    term = f'txid{txid}[Organism] AND "latest refseq"[filter]'
    id_list = []
    retstart = 0
    retmax = 50000
    while True:
        handle = Entrez.esearch(
            db="assembly", term=term, retstart=retstart, retmax=retmax, idtype="acc"
        )
        record = Entrez.read(handle)
        if len(record["IdList"]) == 0:
            break
        id_list.extend(record["IdList"])
        handle.close()
        retstart += retmax
    return id_list


def get_statistics(statistics_name: list, meta: str) -> list:
    statistics = []
    for name in statistics_name:
        query = re.compile(rf"(?<=Stat category=\"{name}\" sequence_tag=\"all\">)[\d]+")
        number = query.findall(meta)[0]
        statistics.append(number)
    return statistics


def write(table: list, file_name: str) -> None:
    with open(file_name, "w") as file:
        file.write(
            f"Organism\tAssembly status\tRefSeq path\tNumber of contigs\tContig L50\t Contig N50\t Number of scaffolds\tScaffold L50\tScaffold N50\tTotal sequence length\n"
        )
        for row in table:
            file.write("\t".join(row))
            file.write("\n")


@click.command()
@click.option("--txid", default="7742", help="Taxon ID in NCBI")
@click.option(
    "--output_file",
    default="table.tsv",
    help="Output file name. Output table will be in tsv format",
)
def get_table(txid: str, output_file: str) -> None:
    STATISTICS_NAME = [
        "contig_count",
        "contig_l50",
        "contig_n50",
        "scaffold_count",
        "scaffold_l50",
        "scaffold_n50",
        "total_length",
    ]

    id_list = get_id_list(txid)
    table = []

    for uid in id_list:
        handle = Entrez.esummary(db="assembly", id=uid, report="full")
        record = Entrez.read(handle)
        handle.close()
        name = record["DocumentSummarySet"]["DocumentSummary"][0]["Organism"]
        status = record["DocumentSummarySet"]["DocumentSummary"][0]["AssemblyStatus"]
        refseq_path = record["DocumentSummarySet"]["DocumentSummary"][0][
            "FtpPath_RefSeq"
        ]
        row = [name, status, refseq_path]
        meta = record["DocumentSummarySet"]["DocumentSummary"][0]["Meta"]
        row.extend(get_statistics(STATISTICS_NAME, meta))
        table.append(row)

    write(table, output_file)


if __name__ == "__main__":
    get_table()
