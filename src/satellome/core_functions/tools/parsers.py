#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2007-2009 Aleksey Komissarov ( ad3002@gmail.com )
# All rights reserved.
#
# This software is licensed as described in the file COPYING, which
# you should have received as part of this distribution.
"""
Function for various parsing tasks.

- parse_fasta_head(fa_head) -> (P1,P2,P3)
- parse_chromosome_name(head) -> string
- trf_parse_line(line) -> [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, '0', '0']
- trf_parse_param(line) -> string
- trf_parse_head(line) -> string  
"""
import re


def parse_fasta_head(fa_head):
    """Parsing fasta head.

    - head -> (...) -> head parameters
    - Head parsing <= fasta subtypes

    #TODO: Write other variants

    Valid examples of head format:

    - >gi|20928815|ref|NW_003237.1|MmUn_WIFeb01_12612 Mus musculus chromosome Un genomic contig
    - Short variant >(\d*)\t(\S*)
    - >134124-14124\n
    - >134124\n
    - >probe|misat|ref|CAAA01154094|start|991|end|1019
    - Repbase: >TC5A    Mariner/Tc1    Caenorhabditis elegans (split by \t <=== settings)
    - >lcl|HmaUn_WGA106_1 Hydra magnipapillata genomic contig, re...
    - plants: >gi|293886233|dbj|BABO01423189.1|
    """

    head_regexp = re.compile("^>?gi\|(\d+)\|(?:ref|dbj)\|([\w.]+)\|(.*)")
    regexp_obj = head_regexp.search(fa_head)

    head_regexp_lcl = re.compile("^>?lcl\|(.*?) (.*)")
    regexp_obj_lcl = head_regexp_lcl.search(fa_head)

    head_regexp_psu = re.compile("^>?psu\|(.*?) (.*)")
    regexp_obj_psu = head_regexp_psu.search(fa_head)

    head_regexp_short = re.compile("^>(\d+)\t(\S+)")
    regexp_obj_short = head_regexp_short.search(fa_head)

    head_regexp_comp = re.compile("^>(\d+)-(\d+)")
    regexp_obj_comp = head_regexp_comp.search(fa_head)

    head_regexp_number = re.compile("^>(\d+)")
    regexp_obj_number = head_regexp_number.search(fa_head)

    head_regexp_probe_index = re.compile(
        "^>probe\|(.*?)\|ref\|(.*?)\|start\|(\d*?)\|end\|(\d*)"
    )
    regexp_obj_probe_index = head_regexp_probe_index.search(fa_head)

    head_wgs = re.compile("^>?gi\|(\d+)\|\w+?\|([\w.]+)\|(.*)")
    regexp_obj_wgs = head_wgs.search(fa_head)

    head_trace = re.compile("^>gnl\|ti\|(\d+) (.*)")
    regexp_obj_trace = head_trace.search(fa_head)

    if regexp_obj:
        match = regexp_obj.groups()
        return list(match)
    elif regexp_obj_probe_index:
        match_l = []

        match = regexp_obj_probe_index.groups()

        gi = "%s_%s_%s" % (match[1], match[2], match[3])
        desc = "%s_%s_%s_%s" % (match[1], match[2], match[3], match[0])

        match_l.append(gi)
        match_l.append("None")
        match_l.append(desc)
        return list(match_l)
    elif regexp_obj_lcl:
        match_l = []
        match = regexp_obj_lcl.groups()
        match_l.append(match[0])
        match_l.append(match[0])
        match_l.append(match[1])
        return list(match_l)
    elif regexp_obj_psu:
        match_l = []
        match = regexp_obj_psu.groups()
        match_l.append(match[0])
        match_l.append(match[0])
        match_l.append(match[1])
        return list(match_l)
    elif regexp_obj_short:
        match = regexp_obj_short.groups()
        match = list(match)
        match.append("Unknown")
        return list(match)
    elif regexp_obj_comp:
        match = regexp_obj_comp.groups()
        match = list(match)
        match.append("Unknown")
        return list(match)
    elif regexp_obj_number:
        match = regexp_obj_number.groups()
        match = list(match)
        match.append("Unknown")
        match.append("Unknown")
        return list(match)
    elif regexp_obj_wgs:
        match = regexp_obj_wgs.groups()
        return list(match)
    elif regexp_obj_trace:
        match = list(regexp_obj_trace.groups())
        match.append(match[-1])
        return list(match)
    else:
        match = ("Unknown", "Unknown", "Unknown")
        # print "Failed parse sequence head: %s" % (fa_head)
        return list(match)


def parse_chromosome_name(head):
    """Parse chromosome name in ncbi fasta head"""
    # Head -> (...) -> Chromosome name or ""
    # TODO: write parse_chromosome_name function
    try:
        chr0 = re.compile("chromosome ([^, ]+)").findall(head)
        chr1 = re.compile("chromosome (\S+?),").findall(head)
        chr2 = re.compile("chromosome (\S+?),?").findall(head)
        chr3 = re.compile("chr(\S+?) ").findall(head)
        mit = re.compile(" (mitochon\S+?) ").findall(head)
        if chr0:
            return chr0[0]
        if chr1:
            return chr1[0]
        if chr2:
            return chr2[0]
        if chr3:
            return chr3[0]
        if mit:
            return "MT"
        return "?"
    except:
        return "?"


def trf_parse_line(line):
    """Parse TRF data line"""
    line = line.strip()
    groups = re.split("\s", line)

    if groups and len(groups) == 15:
        return list(groups)
    else:
        print(("Failed parse ta: %s" % (line)))
        return [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "0", "0"]


def trf_parse_param(line):
    """TRF parameters line"""
    try:
        res = re.compile("Parameters: ([\d ]*)", re.S).findall(line)[0]
        return res
    except:
        res = "Unknown"
        print(("Failed parse param: %s" % (line)))
        return res


def trf_parse_head(line):
    """Parse TRF head"""
    try:
        res = re.compile("Sequence: (.*?)\n", re.S).findall(line)
        res2 = re.compile("Sequence: (.*)", re.S).findall(line)
        if res:
            return res[0]
        if res2:
            return res2[0]
    except:
        res = "Unknown"
        print(("Failed parse head: %s" % (line)))
        return res


def get_wgs_prefix_from_ref(ref):
    """Function parse WGS prefix from Genbank ref."""
    reg_exp = "([A-Z]+)"
    res = re.search(reg_exp, ref)
    if res:
        return res.group(0)
    else:
        return "UNKN"


def get_wgs_prefix_from_head(head):
    """Function parse WGS prefix from fasta head."""
    reg_exp = "ref.([A-Z]{4,})"
    res = re.search(reg_exp, head)
    if res:
        return res.group(1)
    else:
        reg_exp = "gb.([A-Z]{4,})"
        res = re.search(reg_exp, head)
        if res:
            return res.group(1)
        else:
            return None


def refine_name(i, trf_obj):
    """Refine TRF name
    :param trf_obj: TRF object
    :return: TRF object with refined name
    """
    name = trf_obj.trf_head.split()
    if len(name):
        name = name[0]
    else:
        name = name
    trf_obj.trf_id = f"{name}_{trf_obj.trf_l_ind}_{trf_obj.trf_r_ind}"
    trf_obj.id = f"AGT{(i+1) * 100:013d}"
    trf_obj.trf_consensus = trf_obj.trf_consensus.upper()
    trf_obj.trf_array = trf_obj.trf_array.upper()
    return trf_obj
