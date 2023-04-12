#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2007-2009 Aleksey Komissarov ( ad3002@gmail.com )
# All rights reserved.
#
# This software is licensed as described in the file COPYING, which
# you should have received as part of this distribution.
"""
Function for simple sequence analysis.
- get_revcomp(sequence)
- get_gc(sequence)
- check_gapped(sequence)
- clear_sequence(sequence)
- restriction(sequence, pattern, end="")
- restriction_fragments_n(sequence, pattern)
- get_subseq(seq, start, end)
- check_synonims(seq_a, seq_b)
- random_mutation(seq, n, alphabet="actgn ")
    
TODO: Not implemented:
- get_bin_sequence(sequence)
    
"""
import re
import random
from collections import defaultdict

REVCOMP_DICTIONARY = dict(zip('ATCGNatcgn~[]', 'TAGCNtagcn~]['))

def get_revcomp(sequence):
    '''Return reverse complementary sequence.
    >>> complementary('AT CG')
    'CGAT'
    '''
    return ''.join(REVCOMP_DICTIONARY.get(nucleotide, '') for nucleotide in reversed(sequence))

def clear_sequence(sequence, lower=True):
    """ Clear sequence (full alphabet):
    
    - lower case
    - \s -> ""
    - [^actgn] -> ""
    """
    sequence = sequence.strip()
    if lower:
        sequence = sequence.lower().strip()
    sequence = re.sub("\s+", "", sequence)
    return re.sub("[^actgnuwsmkrybdhvACTGNUWSMKRYBDHV\-]", "", sequence)

def get_gc(sequence):
    ''' Count GC content.
    '''
    length = len(sequence)
    if not length:
        return 0
    count_c = sequence.count('c')+sequence.count('C')
    count_g = sequence.count('g')+sequence.count('G')
    gc = float(count_c + count_g) / float(length)
    return float(gc)

def check_gapped(sequence):
    """ Check for N in sequence. Return n(with N) or w(whole).
    """
    w_regexp = re.compile('n|N')
    regexp_obj = w_regexp.search(sequence)
    if (regexp_obj):
        return True
    else:
        return False
    

def get_int_gc(sequence):
    ''' Get GC content from 0 to 100.
    '''
    gc = get_gc(sequence)
    return int(100*round(gc, 2))


def get_shifts_variants(sequence):
    '''
    '''
    shifts = set()
    for i in range(len(sequence)):
        shifts.add(sequence[i:]+sequence[:i])
    return list(shifts)


def sort_dictionary_by_value(d, reverse=False):
    ''' Sort dictionary by value. Retrun list of (v, k) pairs.'''
    result = [(v, k) for k, v in d.items()]
    result.sort(reverse=reverse, key=lambda x: x[0])
    return result


def remove_consensus_redundancy(trf_objs):
    ''' Take a minimal sequence from lexicographically sorted rotations of sequence and its reverse complement
        Example: ACT, underlined - reverse complement seqeunces
        ACT, AGT, CTA, GTA, TAC, TAG
        Find all possible multimers, e.g. replace GTAGTAGTA consensus sequence with ACT
        Return:
        1) list of sorted TRs
        2) list of (df, consensus) pairs
    '''
    # sort by length
    consensuses = [x.trf_consensus for x in trf_objs]
    consensuses = list(set(consensuses))
    consensuses.sort(key=lambda x: len(x))
    length2consensuses = {}
    # print "Group monomers by length and GC"
    for i, monomer in enumerate(consensuses):
        n = len(monomer)
        length2consensuses.setdefault(n, {})
        gc = get_int_gc(monomer)
        length2consensuses[n].setdefault(gc, [])
        length2consensuses[n][gc].append(i)
    # print "Iterate over consensuses"
    N = len(consensuses)
    result_rules = {}
    for i, monomer in enumerate(consensuses):
        if not monomer:
            continue
        # print i, N, "\r",
        if monomer in result_rules:
            continue
        gc = get_int_gc(monomer)
        base = len(monomer)
        n = base
        # maximal consensus length from TRF is 2000 bp
        variants = set(get_shifts_variants(monomer) + get_shifts_variants(get_revcomp(monomer)))
        if not variants:
            raise Exception("Wrong monomer sequence for '%s'" % monomer)
        lex_consensus = min(variants)        
        while n <= 2020:
            if n in length2consensuses and gc in  length2consensuses[n]: 
                for k in length2consensuses[n][gc]:
                    monomer_b = consensuses[k]
                    if monomer_b in result_rules:
                        continue
                    s = n // base
                    v = set()
                    for p in range(s):
                        v.add(monomer_b[p*base:(p+1)*base])            
                    if len(v) > 1:
                        continue
                    item = v.pop()
                    if item in variants:
                        # if consensuses[k] != lex_consensus:
                            # print
                            # print i, base, consensuses[k], "->", lex_consensus
                        result_rules[consensuses[k]] = lex_consensus
                        
            n += base
    # print
    # print "Fix momomers"
    variants2df = defaultdict(int)
    for i,trf_obj in enumerate(trf_objs):
        if not trf_obj.trf_consensus:
            trf_objs[i] = None
            continue
        if trf_obj.trf_consensus in result_rules:
            variants2df[result_rules[trf_obj.trf_consensus]] += 1
        else:
            message = "Error key with length", len(trf_obj.trf_consensus), trf_obj.trf_consensus
            # print str(trf_obj)
            # print message
            raise Exception(message)
        trf_obj.trf_consensus = result_rules[trf_obj.trf_consensus]
    # print "Sort families by df..."
    variants2df = sort_dictionary_by_value(variants2df, reverse=True)
    trf_objs = [x for x in trf_objs if x is not None]
    return trf_objs, variants2df