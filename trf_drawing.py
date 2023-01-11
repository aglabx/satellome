#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 02.12.2022
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import re
import pandas as pd
from trseeker.seqio.fasta_file import sc_iter_fasta_brute


CENPB_REGEXP = re.compile(r'.ttcg....a..cggg.')
TELOMERE_REGEXP = re.compile(r'ttagggttagggttagggttagggttaggg')
ANURA_REGEXP = re.compile("chromosome\: (.*)")


def scaffold_length(fasta_file, lenght_cutoff=100000, name_regexp=ANURA_REGEXP):
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
        scaffold_length.append((name, 1, len(seq)))
    scaffold_df = pd.DataFrame(scaffold_length,columns=['scaffold', 'start', 'end'])
    scaffold_df.sort_values(by=['end'], inplace=True, ascending=False)
    return scaffold_df


def read_table(trf_file):
    ''' Function that convert Aleksey script's trf table to csv.
    '''
    data = pd.read_csv(trf_file, sep='\t', names=['project', '1', 'id', '2', '3', '4', 'start', 'end', 'period', '5', 'pmatch','6', '7', 'mono', 'array', 'gc', '8', '9', 'scaffold', '10', 'length', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26'])
    data.drop(columns=['project', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26'], inplace=True)
    data['mono*3'] = data['mono']*3
    data['centromere'] = [1 if CENPB_REGEXP.findall(i) else 0 for i in data['array']]
    data['telomere'] = [1 if TELOMERE_REGEXP.findall(i) else 0 for i in data['array']]
    data['final_id'] = [f"{x['scaffold']}_{x['id']}" for i, x in data.iterrows()]
    
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
