#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.02.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import argparse
from trf_clusters import *

def main():
    args = get_args()
    trf_file = args.input
    output_folder = args.output
    lenght_cutoff = args.cutoff
    fasta_file = args.fasta
    enhance = args.enhance
    taxon = args.taxon

    chm2name = None
    
    draw_all(trf_file, fasta_file, chm2name, output_folder, taxon, lenght_cutoff=lenght_cutoff, level=1, enhance=enhance)
    

def get_args():
    parser = argparse.ArgumentParser(description='Tandem repeats classification')
    parser.add_argument('-f', '--fasta', type=str, help='Fasta file')
    parser.add_argument('-i', '--input', type=str, help='TRF file')
    parser.add_argument('-o', '--output', type=str, help='Output folder')
    parser.add_argument('-c', '--cutoff', type=int, default=1000000, help='Minimal length of scaffold')
    parser.add_argument('-e', '--enhance', type=int, default=1000000, help='Enhance level for TRs visualization')
    parser.add_argument('-t', '--taxon', type=str, default="Unk", help='Taxon name')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()