#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 02.12.2022
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import re
import pandas as pd
from trseeker.seqio.fasta_file import sc_iter_fasta_brute
import math

CENPB_REGEXP = re.compile(r'.ttcg....a..cggg.')
TELOMERE_REGEXP = re.compile(r'ttagggttagggttagggttagggttaggg')
ANURA_REGEXP = re.compile("chromosome\: (.*)")

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

def scaffold_length_sort_dict(fasta_file, chm2name, lenght_cutoff=100000, name_regexp=ANURA_REGEXP):
    ''' Function that calculates length of scaffolds 
        and return table with scaffold data from fasta file
    '''
    scaffold_length = []
    for header, seq in sc_iter_fasta_brute(fasta_file):
        name = header[1:].split()[0]
        if len(seq) < lenght_cutoff:
            continue
        new_name = re.findall(ANURA_REGEXP, header)
        if new_name:
            name = new_name[0]
            
        name = chm2name[name]
        scaffold_length.append((name, 1, len(seq)))
        
    scaffold_length.sort(key=lambda x: sort_chrm(x[0]))
    
    scaffold_df = pd.DataFrame(scaffold_length, columns=['scaffold', 'start', 'end'])
    return scaffold_df


def scaffold_length_sort_length(fasta_file, lenght_cutoff=100000, name_regexp=ANURA_REGEXP):
    ''' Function that calculates length of scaffolds 
        and return table with scaffold data from fasta file
    '''
    scaffold_length = []
    for header, seq in sc_iter_fasta_brute(fasta_file):
        name = header[1:].split()[0]
        if len(seq) < lenght_cutoff:
            continue
        new_name = re.findall(ANURA_REGEXP, header)
        if new_name:
            name = new_name[0]
            
        name = chm2name[name]
        scaffold_length.append((name, 1, len(seq)))
        
    
    scaffold_df = pd.DataFrame(scaffold_length, columns=['scaffold', 'start', 'end'])
    scaffold_df.sort_values(by=['end'], inplace=True, ascending=False)
    return scaffold_df


def read_table(trf_file):
    ''' Function that convert Aleksey script's trf table to csv.
    '''
    data = pd.read_csv(trf_file, sep='\t', names=['project', '1', 'id', '2', '3', '4', 'start', 'end', 'period', '5', 'pmatch','6', '7', 'mono', 'array', 'gc', '8', '9', 'scaffold', '10', 'length', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26'], low_memory=False)
    data.drop(columns=['project', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26'], inplace=True)
    data['mono*3'] = data['mono']*3
    data['centromere'] = [1 if CENPB_REGEXP.findall(i) else 0 for i in data['array']]
    data['telomere'] = [1 if TELOMERE_REGEXP.findall(i) else 0 for i in data['array']]
    data['final_id'] = [f"{x['scaffold']}_{x['id']}" for i, x in data.iterrows()]
    data['class_name'] = ["CENPB" if x['centromere'] else "UNK" for i, x in data.iterrows()]
    data['class_name'] = ["TEL" if x['telomere'] else x["class_name"] for i, x in data.iterrows()]
    data['family_name'] = None
    data['locus_name'] = None
    data['log_length'] = [math.log(x['length']) for i, x in data.iterrows()]
    data['scaffold'] = [x['scaffold'].split()[0] for i, x in data.iterrows()]
    return data


def check_patterns(data):
    '''
    '''
    centromers = data.loc[data['centromere'] == 1]
    telomers = data.loc[data['telomere'] == 1]
    return (centromers, telomers)


def draw_genome(fasta_file, trf_file):
    '''
    '''
    trf_df = read_table(trf_file)
    scaffold_df = scaffold_length(fasta_file, lenght_cutoff=100000, name_regexp=ANURA_REGEXP)
