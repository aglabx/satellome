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
import logging
import re

logger = logging.getLogger(__name__)


def parse_fasta_head(fa_head):
    """
    Parse FASTA header into structured components (GI, accession, description).

    Supports multiple FASTA header formats from NCBI, custom databases, and other
    sources. Tries multiple regular expression patterns to extract GI number,
    accession, and description from the header.

    Args:
        fa_head (str): FASTA header line (with or without leading '>')

    Returns:
        list: [gi, accession, description] where:
            - gi (str): GenBank identifier or custom ID
            - accession (str): Sequence accession number or reference
            - description (str): Sequence description or name
            Returns ["Unknown", "Unknown", "Unknown"] if no pattern matches.

    Supported formats:
        - NCBI RefSeq: >gi|20928815|ref|NW_003237.1|MmUn_WIFeb01_12612...
        - Local sequences: >lcl|HmaUn_WGA106_1 Hydra magnipapillata...
        - Short variant: >134124\\tSEQNAME
        - Numeric range: >134124-14124
        - Numeric only: >134124
        - Probe format: >probe|misat|ref|CAAA01154094|start|991|end|1019
        - WGS entries: >gi|293886233|dbj|BABO01423189.1|
        - Trace archives: >gnl|ti|123456 description

    Example:
        >>> parse_fasta_head(">gi|123|ref|NW_001|description text")
        ['123', 'NW_001', 'description text']
        >>> parse_fasta_head(">12345")
        ['12345', 'Unknown', 'Unknown']
    """

    head_regexp = re.compile(r"^>?gi\|(\d+)\|(?:ref|dbj)\|([\w.]+)\|(.*)")
    regexp_obj = head_regexp.search(fa_head)

    head_regexp_lcl = re.compile(r"^>?lcl\|(.*?) (.*)")
    regexp_obj_lcl = head_regexp_lcl.search(fa_head)

    head_regexp_psu = re.compile(r"^>?psu\|(.*?) (.*)")
    regexp_obj_psu = head_regexp_psu.search(fa_head)

    head_regexp_short = re.compile(r"^>(\d+)\t(\S+)")
    regexp_obj_short = head_regexp_short.search(fa_head)

    head_regexp_comp = re.compile(r"^>(\d+)-(\d+)")
    regexp_obj_comp = head_regexp_comp.search(fa_head)

    head_regexp_number = re.compile(r"^>(\d+)")
    regexp_obj_number = head_regexp_number.search(fa_head)

    head_regexp_probe_index = re.compile(
        r"^>probe\|(.*?)\|ref\|(.*?)\|start\|(\d*?)\|end\|(\d*)"
    )
    regexp_obj_probe_index = head_regexp_probe_index.search(fa_head)

    head_wgs = re.compile(r"^>?gi\|(\d+)\|\w+?\|([\w.]+)\|(.*)")
    regexp_obj_wgs = head_wgs.search(fa_head)

    head_trace = re.compile(r"^>gnl\|ti\|(\d+) (.*)")
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
        return list(match)


def parse_chromosome_name(head):
    """
    Extract chromosome name from NCBI FASTA header.

    Tries multiple regular expression patterns to identify chromosome names,
    including special handling for mitochondrial DNA.

    Args:
        head (str): FASTA header line containing chromosome information

    Returns:
        str: Chromosome name (e.g., "1", "X", "Y", "MT") or "?" if not found

    Supported patterns:
        - "chromosome 1", "chromosome X" → "1", "X"
        - "chr1 ", "chrX " → "1", "X"
        - "mitochondrion", "mitochondrial" → "MT"

    Example:
        >>> parse_chromosome_name("Mus musculus chromosome 1, GRCm38")
        '1'
        >>> parse_chromosome_name("chrX complete sequence")
        'X'
        >>> parse_chromosome_name("mitochondrion genome")
        'MT'
    """
    # Head -> (...) -> Chromosome name or ""
    try:
        chr0 = re.compile(r"chromosome ([^, ]+)").findall(head)
        chr1 = re.compile(r"chromosome (\S+?),").findall(head)
        chr2 = re.compile(r"chromosome (\S+?),?").findall(head)
        chr3 = re.compile(r"chr(\S+?) ").findall(head)
        mit = re.compile(r" (mitochon\S+?) ").findall(head)
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
    except (AttributeError, IndexError, TypeError) as e:
        logger.error(f"Error parsing chromosome name from {head}: {e}")
        return "?"


def trf_parse_line(line):
    """
    Parse TRF (Tandem Repeat Finder) data line into fields.

    TRF output lines contain 15 space-separated fields describing a tandem repeat.

    Args:
        line (str): Single line from TRF output containing repeat data

    Returns:
        list: List of 15 string values representing TRF fields:
            [start, end, period, copy_number, consensus_size, percent_matches,
             percent_indels, score, A_pct, C_pct, G_pct, T_pct, entropy,
             consensus_sequence, repeat_sequence]
            Returns list of zeros and "0" if parsing fails.

    Example:
        >>> trf_parse_line("100 200 5 20.0 5 95 5 190 25 25 25 25 1.50 AAAAT AAAAT...")
        ['100', '200', '5', '20.0', '5', '95', '5', '190', '25', '25', '25', '25', '1.50', 'AAAAT', 'AAAAT...']
    """
    line = line.strip()
    groups = re.split(r"\s", line)

    if groups and len(groups) == 15:
        return list(groups)
    else:
        logger.error(f"Failed parse ta: {line}")
        return [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "0", "0"]


def trf_parse_param(line):
    """
    Extract TRF run parameters from parameters line.

    TRF output includes a "Parameters:" line with algorithm settings.

    Args:
        line (str): TRF parameters line (e.g., "Parameters: 2 7 7 80 10 50 500")

    Returns:
        str: Parameter values as space-separated string, or "Unknown" if parsing fails

    Example:
        >>> trf_parse_param("Parameters: 2 7 7 80 10 50 500")
        '2 7 7 80 10 50 500'
    """
    try:
        res = re.compile(r"Parameters: ([\d ]*)", re.S).findall(line)[0]
        return res
    except (IndexError, AttributeError) as e:
        logger.error(f"Failed parse param: {line}, error: {e}")
        return "Unknown"


def trf_parse_head(line):
    """
    Extract sequence name from TRF sequence header line.

    TRF output includes "Sequence: <name>" lines before each sequence's repeats.

    Args:
        line (str): TRF sequence header line (e.g., "Sequence: chr1")

    Returns:
        str: Sequence name, or None if parsing fails

    Example:
        >>> trf_parse_head("Sequence: chr1\\n")
        'chr1'
        >>> trf_parse_head("Sequence: scaffold_123")
        'scaffold_123'
    """
    try:
        res = re.compile("Sequence: (.*?)\n", re.S).findall(line)
        res2 = re.compile("Sequence: (.*)", re.S).findall(line)
        if res:
            return res[0]
        if res2:
            return res2[0]
        return None  # Return None if no pattern matches
    except (IndexError, AttributeError) as e:
        logger.error(f"Failed parse head: {line}, error: {e}")
        return None


def get_wgs_prefix_from_ref(ref):
    """
    Extract WGS (Whole Genome Shotgun) prefix from GenBank reference.

    Extracts the uppercase letter prefix from GenBank WGS accessions,
    which identify the sequencing project.

    Args:
        ref (str): GenBank reference/accession (e.g., "AAAA01000001.1")

    Returns:
        str: WGS prefix (uppercase letters), or "UNKN" if not found

    Example:
        >>> get_wgs_prefix_from_ref("AAAA01000001.1")
        'AAAA'
        >>> get_wgs_prefix_from_ref("NW_003237.1")
        'NW'
    """
    reg_exp = "([A-Z]+)"
    res = re.search(reg_exp, ref)
    if res:
        return res.group(0)
    else:
        return "UNKN"


def get_wgs_prefix_from_head(head):
    """
    Extract WGS prefix from FASTA header.

    Searches for WGS project prefixes in RefSeq (ref|) or GenBank (gb|)
    formatted FASTA headers.

    Args:
        head (str): FASTA header line

    Returns:
        str: WGS prefix (4+ uppercase letters), or None if not found

    Example:
        >>> get_wgs_prefix_from_head(">gi|123|ref|AAAA01000001.1|...")
        'AAAA'
        >>> get_wgs_prefix_from_head(">gi|456|gb|BBBB02000001.1|...")
        'BBBB'
    """
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
    """
    Generate standardized IDs and normalize sequences for TRF object.

    Creates unique identifiers based on genomic location and ensures
    consensus/array sequences are uppercase.

    Args:
        i (int): Index number for generating unique ID (0-based)
        trf_obj: TRF model object with attributes trf_head, trf_l_ind,
                trf_r_ind, trf_consensus, trf_array

    Returns:
        TRF object: Modified TRF object with updated fields:
            - trf_id: "{chromosome}_{start}_{end}" format
            - id: "AGT" + zero-padded number (e.g., "AGT0000000000100")
            - trf_consensus: Uppercase consensus sequence
            - trf_array: Uppercase tandem repeat array

    Example:
        >>> trf = TRModel()
        >>> trf.trf_head, trf.trf_l_ind, trf.trf_r_ind = "chr1 description", 1000, 1100
        >>> trf.trf_consensus, trf.trf_array = "acgt", "acgtacgt"
        >>> refined = refine_name(0, trf)
        >>> refined.trf_id
        'chr1_1000_1100'
        >>> refined.id
        'AGT0000000000100'
        >>> refined.trf_consensus
        'ACGT'
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
