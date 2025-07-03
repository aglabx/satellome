#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 02.12.2022
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import math
import re

import pandas as pd
from tqdm import tqdm

from satellome.core_functions.io.fasta_file import sc_iter_fasta_brute

CENPB_REGEXP = re.compile(r".ttcg....a..cggg.")
TELOMERE_REGEXP = re.compile(r"ttagggttagggttagggttagggttaggg")
CHRM_REGEXP = re.compile("chromosome\: (.*)")

chm2name = {
    "NC_060925.1": "Chr1",
    "NC_060926.1": "Chr2",
    "NC_060927.1": "Chr3",
    "NC_060928.1": "Chr4",
    "NC_060929.1": "Chr5",
    "NC_060930.1": "Chr6",
    "NC_060931.1": "Chr7",
    "NC_060947.1": "ChrX",
    "NC_060932.1": "Chr8",
    "NC_060933.1": "Chr9",
    "NC_060934.1": "Chr10",
    "NC_060935.1": "Chr11",
    "NC_060936.1": "Chr12",
    "NC_060937.1": "Chr13",
    "NC_060938.1": "Chr14",
    "NC_060939.1": "Chr15",
    "NC_060940.1": "Chr16",
    "NC_060941.1": "Chr17",
    "NC_060942.1": "Chr18",
    "NC_060943.1": "Chr19",
    "NC_060944.1": "Chr20",
    "NC_060945.1": "Chr21",
    "NC_060946.1": "Chr22",
    "NC_060948.1": "ChrY",
}


def sort_chrm(name):
    v = name.replace("Chr", "")
    print(v)
    if v == "Y":
        return 24
    if v == "X":
        return 24
    return int(v)


def scaffold_length_sort_dict(
    fasta_file, lenght_cutoff=100000, name_regexp=None, chm2name=None
):
    """Function that calculates length of scaffolds
    and return table with scaffold data from fasta file
    """
    scaffold_length = []
    for header, seq in sc_iter_fasta_brute(fasta_file):
        name = header[1:].split()[0]
        if len(seq) < lenght_cutoff:
            continue
        if name_regexp:
            new_name = re.findall(name_regexp, header)
            if new_name:
                name = new_name[0]
        if chm2name:
            name = chm2name[name]
        scaffold_length.append((name, 1, len(seq)))

    scaffold_length.sort(key=lambda x: sort_chrm(x[0]))

    scaffold_df = pd.DataFrame(scaffold_length, columns=["scaffold", "start", "end"])
    return scaffold_df


def scaffold_length_sort_length(
    fasta_file, lenght_cutoff=100000, name_regexp=None, chm2name=None
):
    """Function that calculates length of scaffolds
    and return table with scaffold data from fasta file
    """
    scaffold_length = []
    for header, seq in sc_iter_fasta_brute(fasta_file):
        name = header[1:].split()[0]
        if len(seq) < lenght_cutoff:
            continue
        if name_regexp:
            new_name = re.findall(name_regexp, header)
            if new_name:
                name = new_name[0]
        if chm2name:
            name = chm2name[name]
        scaffold_length.append((name, 1, len(seq)))

    scaffold_df = pd.DataFrame(scaffold_length, columns=["scaffold", "start", "end"])
    scaffold_df.sort_values(by=["end"], inplace=True, ascending=False)
    return scaffold_df


def read_trf_file(trf_file):
    """Function that convert Aleksey script's trf table to csv."""

    
    data = pd.read_csv(
        trf_file,
        sep="\t",
        low_memory=False,
    )
    data["start"] = data["trf_l_ind"]
    data["end"] = data["trf_r_ind"]
    data["period"] = data["trf_period"]
    data["pmatch"] = data["trf_pmatch"]
    data["mono"] = data["trf_consensus"]
    data["array"] = data["trf_array"]
    data["gc"] = data["trf_array_gc"]
    data["scaffold"] = data["trf_head"]
    data["length"] = data["trf_array_length"]
    data["seq"] = data["array"]
    data["mono*3"] = data["mono"] * 3
    data["centromere"] = [1 if CENPB_REGEXP.findall(i) else 0 for i in data["array"]]
    data["telomere"] = [1 if TELOMERE_REGEXP.findall(i) else 0 for i in data["array"]]
    data["final_id"] = [f"{x['scaffold']}_{x['id']}" for i, x in data.iterrows()]
    data["class_name"] = [
        "CENPB" if x["centromere"] else "UNK" for i, x in data.iterrows()
    ]
    data["class_name"] = [
        "TEL" if x["telomere"] else x["class_name"] for i, x in data.iterrows()
    ]
    data["family_name"] = None
    data["locus_name"] = None
    data["log_length"] = [math.log(x["length"]) for i, x in data.iterrows()]
    data["scaffold"] = [x["scaffold"].split()[0] for i, x in data.iterrows()]
    return data


def check_patterns(data):
    """ """
    centromers = data.loc[data["centromere"] == 1]
    telomers = data.loc[data["telomere"] == 1]
    return (centromers, telomers)


def get_gaps_annotation(fasta_file, genome_size, lenght_cutoff=100000):
    """Function that finding all gaps."""
    gaps = []

    with tqdm(total=genome_size, desc="Find gaps") as pbar:
        for header, seq in sc_iter_fasta_brute(fasta_file):
            name = header[1:].split()[0]
            if len(seq) < lenght_cutoff:
                continue
            in_gap = False
            gap_start = None
            # print(name)
            for i in range(len(seq)):
                
                if seq[i] == "N":
                    if not in_gap:
                        in_gap = True
                        gap_start = i
                    continue
                if in_gap:
                    in_gap = False
                    gaps.append([name, gap_start, i, abs(gap_start - i)])
            if in_gap:
                in_gap = False
                gaps.append([name, gap_start, i, abs(gap_start - i)])
            pbar.update(len(seq))
    return gaps


def get_gaps_annotation_re(fasta_file, genome_size, lenght_cutoff=100000):
    """Function that finding all gaps."""
    gaps = []
    with tqdm(total=genome_size, desc="Find gaps") as pbar:
        for header, seq in sc_iter_fasta_brute(fasta_file):
            pbar.update(len(seq))
            name = header[1:].split()[0]
            if len(seq) < lenght_cutoff:
                continue
            print(name)
            hits = re.findall("N+", seq)
            print(hits)
            for pos, item in hits:
                gaps.append((name, pos, pos + len(item)))
    return gaps
