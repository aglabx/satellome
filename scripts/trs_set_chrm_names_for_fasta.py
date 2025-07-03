#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.02.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import argparse
import re

from satellome.core_functions.io.fasta_file import sc_iter_fasta_brute


def main(args):

    header2name = {}
    header2size = {}
    for i, (header, seq) in enumerate(sc_iter_fasta_brute(args.fasta)):
        if len(seq) < args.cutoff:
            continue
        print(header)
        name = "chr%s" % i
        suggest = re.findall(r"chromosome (\d+)", header)
        if suggest:
            name = "chr%s" % suggest[0]
            if args.manual:
                name = input(f"Set name ({name})?") or name
        else:
            name = input(f"Set name ({name})?") or name
        header2name[header] = name
        header2size[header] = len(seq)
        print(f"Set name {name} with length ({len(seq)} bp)")

    with open(args.output, "w") as fh:
        for header, name in header2name.items():
            fh.write("%s\t%s\t%s\n" % (header, name, header2size[header]))


def get_args():
    parser = argparse.ArgumentParser(
        description="Create names for chrmosomes for fasta file"
    )
    parser.add_argument("-f", "--fasta", help="Fasta file", required=True)
    parser.add_argument("-o", "--output", help="Output file", required=True)
    parser.add_argument("--manual", help="Ask each contig", default=False)
    parser.add_argument(
        "-c", "--cutoff", type=int, default=1000000, help="Minimal length of scaffold"
    )
    args = parser.parse_args()
    return args


if __name__ == "__main__":

    args = get_args()
    main(args)
