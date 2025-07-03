#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 10.10.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import argparse
import os
from satellome.core_functions.io.fasta_file import sc_iter_fasta_brute
from intervaltree import IntervalTree
import re
from tqdm import tqdm

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Check mammalian telomeric repeats in assembly")
    parser.add_argument("-i", "--input", help="Input fasta file", required=True)
    parser.add_argument("-o", "--output", help="Output report file", required=True)
    args = vars(parser.parse_args())

    input_fasta_file = args["input"]
    output_report_file = args["output"]

    # forward_telomeric_repeats = "[TAACCC]{30,}"
    # reverse_telomeric_repeats = "[TTAGGG]{30,}"

    flank_size = 50

    forward_telomeric_repeats = "TAACCC" * 30
    reverse_telomeric_repeats = "TTAGGG" * 30

    forward_regexp = re.compile(forward_telomeric_repeats)
    reverse_regexp = re.compile(reverse_telomeric_repeats)
    
    dataset = {}
    header2length = {}

    fasta_length = os.path.getsize(input_fasta_file)

    
    with tqdm(total=fasta_length, desc="Searching") as pbar:
        for header, seq in sc_iter_fasta_brute(input_fasta_file):
            header2length[header] = len(seq)

            seq = seq.upper()

            forward_intervals = []
            for match in re.finditer(forward_regexp, seq):
                start = max(0, match.start() - flank_size)
                end = min(match.end() + flank_size, header2length[header])
                forward_intervals.append((start, end))
            forward_itree = IntervalTree.from_tuples(forward_intervals)
            forward_itree.merge_overlaps()

            reverse_intervals = []
            for match in re.finditer(reverse_regexp, seq):
                start = max(0, match.start() - flank_size)
                end = min(match.end() + flank_size, header2length[header])
                reverse_intervals.append((start, end))
            reverse_itree = IntervalTree.from_tuples(reverse_intervals)
            reverse_itree.merge_overlaps()

            dataset[header] = {
                "forward": sorted(forward_itree),
                "reverse": sorted(reverse_itree),
            }

            pbar.update(len(seq)+1)
            pbar.update(len(header)+1)
        pbar.total = pbar.n



    with open(output_report_file, "w") as fh:
        fh.write(f"#header\ttotal_length\direction\tbegin\tend\tlength\n")
        for header, data in dataset.items():
            forward_intervals = data["forward"]
            reverse_intervals = data["reverse"]

            fh.write(f"{header}\t{header2length[header]}\t\t\t\n")
            
            for interval in forward_intervals:
                length = interval.end - interval.begin
                fh.write(f"\t\tforward\t{interval.begin}\t{interval.end}\t{length}\n")
            
            for interval in reverse_intervals:
                length = interval.end - interval.begin
                fh.write(f"\t\treverse\t\t{interval.begin}\t{interval.end}\t{length}\n")

            
        
