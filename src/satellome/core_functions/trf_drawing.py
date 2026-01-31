#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 02.12.2022
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import logging
import math
import re
import csv
import sys

from tqdm import tqdm

from satellome.core_functions.io.fasta_file import sc_iter_fasta_brute

logger = logging.getLogger(__name__)

CENPB_REGEXP = re.compile(r".ttcg....a..cggg.")
TELOMERE_REGEXP = re.compile(r"ttagggttagggttagggttagggttaggg")
CHRM_REGEXP = re.compile(r"chromosome\: (.*)")

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
    logger.debug(v)
    if v == "Y":
        return 24
    if v == "X":
        return 24
    return int(v)


def scaffold_length_sort_dict(
    fasta_file, lenght_cutoff=100000, name_regexp=None, chm2name=None
):
    """
    Calculate scaffold lengths and return sorted scaffold data from FASTA file.

    Reads a FASTA file, filters scaffolds by minimum length, optionally
    applies name transformations, and returns scaffold coordinates sorted
    by chromosome name.

    Args:
        fasta_file (str): Path to FASTA file containing genomic scaffolds
        lenght_cutoff (int, optional): Minimum scaffold length in bp.
                                       Scaffolds shorter than this are excluded.
                                       Defaults to 100000 (100kb).
        name_regexp (str, optional): Regular expression pattern to extract
                                     scaffold name from FASTA header.
                                     If None, uses first space-delimited word.
                                     Defaults to None.
        chm2name (dict, optional): Dictionary mapping scaffold names to
                                   chromosome names for renaming.
                                   Defaults to None.

    Returns:
        dict: Dictionary with three keys, each containing a list:
            - 'scaffold' (list of str): Scaffold/chromosome names
            - 'start' (list of int): Start coordinates (all set to 1)
            - 'end' (list of int): End coordinates (scaffold lengths)
            All lists are sorted by chromosome name using sort_chrm().

    Example:
        >>> data = scaffold_length_sort_dict("genome.fasta", lenght_cutoff=50000)
        >>> data['scaffold']
        ['Chr1', 'Chr2', 'Chr3']
        >>> data['end']
        [150000, 200000, 180000]
    """
    scaffolds = []
    starts = []
    ends = []

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
        scaffolds.append(name)
        starts.append(1)
        ends.append(len(seq))

    # Sort by chromosome name
    sorted_data = sorted(zip(scaffolds, starts, ends), key=lambda x: sort_chrm(x[0]))

    return {
        "scaffold": [x[0] for x in sorted_data],
        "start": [x[1] for x in sorted_data],
        "end": [x[2] for x in sorted_data]
    }


def scaffold_length_sort_length(
    fasta_file, lenght_cutoff=100000, name_regexp=None, chm2name=None
):
    """Function that calculates length of scaffolds
    and return dict with scaffold data from fasta file, sorted by length

    Returns:
        dict with keys: 'scaffold', 'start', 'end' (lists), sorted by end descending
    """
    scaffolds = []
    starts = []
    ends = []

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
        scaffolds.append(name)
        starts.append(1)
        ends.append(len(seq))

    # Sort by length (descending)
    sorted_data = sorted(zip(scaffolds, starts, ends), key=lambda x: x[2], reverse=True)

    return {
        "scaffold": [x[0] for x in sorted_data],
        "start": [x[1] for x in sorted_data],
        "end": [x[2] for x in sorted_data]
    }


def read_trf_file(trf_file):
    """Function that convert Aleksey script's trf table to list of dicts.

    Returns:
        list of dicts, each representing one TRF record
    """
    # Increase CSV field size limit for large satellite arrays (can be several megabases)
    csv.field_size_limit(sys.maxsize)

    # Define expected field names for backward compatibility with files without header
    TRF_FIELDNAMES = ["project", "trf_id", "trf_head", "trf_l_ind", "trf_r_ind", "trf_period", "trf_n_copy",
                      "trf_pmatch", "trf_pvar", "trf_entropy", "trf_consensus", "trf_array",
                      "trf_array_gc", "trf_consensus_gc", "trf_array_length", "trf_joined", "trf_family", "trf_ref_annotation"]

    data = []

    with open(trf_file, 'r') as f:
        # Filter out comment lines starting with '#'
        lines = list(line for line in f if not line.startswith('#'))

        # Check if first line is a header row
        if lines and lines[0].strip().startswith("project"):
            # File has header row, use standard DictReader
            reader = csv.DictReader(iter(lines), delimiter='\t')
        else:
            # File has no header row, provide field names
            reader = csv.DictReader(iter(lines), fieldnames=TRF_FIELDNAMES, delimiter='\t')

        for row in reader:
            # Create computed fields
            record = dict(row)  # Copy all original fields

            # Add renamed/computed fields
            record["start"] = row.get("trf_l_ind")
            record["end"] = row.get("trf_r_ind")
            record["period"] = row.get("trf_period")
            record["pmatch"] = row.get("trf_pmatch")
            record["mono"] = row.get("trf_consensus")
            record["array"] = row.get("trf_array")
            record["gc"] = row.get("trf_array_gc")
            record["scaffold"] = row.get("trf_head")
            record["length"] = row.get("trf_array_length")
            record["seq"] = record["array"]
            record["mono*3"] = record["mono"] * 3 if record.get("mono") else None

            # Pattern matching
            array_val = record.get("array") or ""
            record["centromere"] = 1 if array_val and CENPB_REGEXP.findall(array_val) else 0
            record["telomere"] = 1 if array_val and TELOMERE_REGEXP.findall(array_val) else 0

            # Computed fields that need other fields first
            record["final_id"] = f"{record['scaffold']}_{record.get('id', '')}"
            record["class_name"] = "CENPB" if record["centromere"] else "UNK"
            record["class_name"] = "TEL" if record["telomere"] else record["class_name"]
            record["family_name"] = None
            record["locus_name"] = None

            # Numeric computations
            length_val = record.get("length")
            if length_val:
                try:
                    record["log_length"] = math.log(float(length_val))
                except (ValueError, TypeError):
                    record["log_length"] = None
            else:
                record["log_length"] = None

            # Clean scaffold name
            if record.get("scaffold"):
                record["scaffold"] = record["scaffold"].split()[0]

            data.append(record)

    return data


def check_patterns(data):
    """Filter data for centromere and telomere patterns.

    Args:
        data: list of dicts from read_trf_file()

    Returns:
        tuple of (centromers, telomers) - both are lists of dicts
    """
    centromers = [record for record in data if record.get("centromere") == 1]
    telomers = [record for record in data if record.get("telomere") == 1]
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
            logger.debug(name)
            hits = re.findall("N+", seq)
            logger.debug(hits)
            for pos, item in hits:
                gaps.append((name, pos, pos + len(item)))
    return gaps
