#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 23.02.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
'''
Functions related to ngram (k-mer).
'''
from collections import defaultdict
from sequence_tools import get_revcomp
import zlib
import logging


def generate_ngrams(text, n=12):
    ''' Yields all ngrams of length k from given text. 
    
    - n: ngram length
    '''
    for i in range(0, len(text) - n + 1):
        yield i, text[i:i + n]


def get_ngrams_freq(text, m=None, n=23, k=None):
    ''' Returns m most frequent (ngram of length n, fraction of possible ngrams) tuples for given text.
    
    - m: number of returned ngrams
    - n: ngram length
    '''
    if k:
        n = k
    ngrams = defaultdict(int)
    for pos, ngram in generate_ngrams(text, n=n):
        ngrams[ngram] += 1
    text_length = float(len(text) - n + 1)
    ngrams = [(key, value, value / text_length) for key, value in ngrams.items()]
    ngrams.sort(reverse=True, key=lambda x: x[1])
    if m:
        return ngrams[:m]
    return ngrams


def count_kmer_tfdf(sequence, tf_dict, df_dict, k):
    ''' Update tf and df data with k-mers from given sequence.
    '''
    seen = set()
    local_tf = defaultdict(int)
    local_df = defaultdict(int)
    sequence = sequence.lower()
    for (ngram, tf, nf) in get_ngrams_freq(sequence, n=k):
        if 'n' in ngram:
            continue
        seen.add(ngram)
        tf_dict[ngram] += tf
        local_tf[ngram] += tf
    for ngram in seen:
        df_dict[ngram] += 1
        local_df[ngram] += 1
    return tf_dict, df_dict, local_tf, local_df


def get_kmer_tf_df_for_data(data, k, docids=False, verbose=True):
    '''
    '''
    df = defaultdict(int)
    tf = defaultdict(int)
    kmer2ids = defaultdict(list)
    kmer2freq = defaultdict(list)
    N = len(data)
    if N>100:
        verbose = True
    for i, sequence in enumerate(data):
        if verbose:
            logging.info("Process tf/df: ", i, N, sep=" ")
        tf, df, local_tf, local_df = count_kmer_tfdf(sequence, tf, df, k)
        if docids:
            for key in local_tf:
                kmer2ids[key].append(i)
                kmer2freq[key].append(local_tf[key])
    if verbose:
        logging.info()
    if docids:
        return tf, df, kmer2ids, kmer2freq
    return tf, df


def process_list_to_kmer_index(data, k, docids=True, cutoff=None, verbose=True):
    ''' Get list of string.
    Return list of (kmer, revkmer, tf, df, None, None)
    OR
    Return list of (kmer, revkmer, tf, df, docids, freqs)
    '''
    if docids:
        (tf_dict, df_dict, doc_data, freq_data) = get_kmer_tf_df_for_data(data, k, docids=docids, verbose=verbose)
    else:
        (tf_dict, df_dict) = get_kmer_tf_df_for_data(data, k, docids=docids, verbose=verbose)
    result = []
    seen = set()
    skipped = 0
    added = 0
    if verbose:
        logging.info("Join data...")
    for key in df_dict:
        if verbose:
            logging.info(skipped, added, sep=" ")
        if key in seen:
            continue
        revkey = get_revcomp(key)
        if revkey in seen:
            continue
        if revkey in df_dict:
            df = df_dict[key] + df_dict[revkey]
            tf = tf_dict[key] + tf_dict[revkey]
            if docids:
                ids = doc_data[key] + doc_data[revkey]
                freqs = freq_data[key] + freq_data[revkey]
        else:
            df = df_dict[key]
            tf = tf_dict[key]
            if docids:
                ids = doc_data[key]
                freqs = freq_data[key]
        # skip by df
        if cutoff and df <= cutoff:
            skipped += 1
            continue
        added += 1
        if revkey < key:
            key, revkey = revkey, key
        if docids:
            result.append((key, revkey, tf, df, ids, freqs))
        else:
            result.append((key, revkey, tf, df, None, None))
        seen.add(key)
        seen.add(revkey)
    if verbose:        
        logging.info()
    result.sort(key=lambda x: x[-3], reverse=True)
    return result

def get_zlib_complexity(s):
    ''' Get simple complexity by zlib compression ratio.
    '''
    return float(len(zlib.compress(s)))/len(s)