#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 07.09.2011
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com
"""
Working with tab-delimited TRs datasets files.
TODO: check it.
"""
import os
from collections import defaultdict

from satellome.core_functions.io.tab_file import (sc_iter_simple_tab_file,
                                                 sc_iter_tab_file)
from satellome.core_functions.models.trf_model import (TRModel,
                                                      TRsClassificationModel)


def read_trid2ngrams(annotation_ngram_folder, trf_large_file):
    """Read trid -> [(ngram, tf), (rev_ngram, tf), ...] data."""

    trid2ngrams = {}
    for trf_obj in sc_iter_tab_file(trf_large_file, TRModel):
        file_name = os.path.join(annotation_ngram_folder, "%s.ngram" % trf_obj.trf_id)
        trid2ngrams.setdefault(trf_obj.trf_id, [])

        for data in sc_iter_simple_tab_file(file_name):
            ngram = data[0]
            rev_ngram = data[1]
            tf = float(data[2])
            trid2ngrams[trf_obj.trf_id].append((ngram, tf))
            trid2ngrams[trf_obj.trf_id].append((rev_ngram, tf))
    return trid2ngrams


def read_trid2meta(file_name):
    """Load trid to full index dictionary as string."""

    print("Load trid to full index dictionary")
    trid2meta = {}
    for trf_obj in sc_iter_tab_file(file_name, TRModel):
        trid2meta[trf_obj.trf_id] = str(trf_obj)
    return trid2meta


def get_all_trf_objs(trf_large_file):
    """Return list of trf_obj from given trf_large_file."""
    result = []
    for trf_obj in sc_iter_tab_file(trf_large_file, TRModel):
        result.append(trf_obj)
    return result


def get_all_class_objs(trf_class_file):
    """Return list of class_obj from given trf_class_file."""
    result = []
    for trf_obj in sc_iter_tab_file(trf_class_file, TRsClassificationModel):
        result.append(trf_obj)
    return result


def get_class_objs_dict(trf_class_file):
    """Return list of class_obj from given trf_class_file."""
    result = {}
    for trf_obj in sc_iter_tab_file(trf_class_file, TRsClassificationModel):
        result[trf_obj.trf_id] = trf_obj
    return result


def get_trf_objs_dict(trf_large_file):
    """Return dict of trf_obj from given trf_large_file."""
    result = {}
    for trf_obj in sc_iter_tab_file(trf_large_file, TRModel):
        result[trf_obj.trf_id] = trf_obj
    return result


def get_trfid_obj_dict(trf_large_file):
    """Return dicionary trf_id to trf_obj from given trf_large_file."""
    trs_dataset = get_all_trf_objs(trf_large_file)
    trid2obj = {}
    for trf_obj in trs_dataset:
        trid2obj[trf_obj.trf_id] = trf_obj
    return trid2obj


def save_trs_dataset(trs_dataset, output_file, dataset_id=None):
    """Save trs dataset to file."""
    if isinstance(trs_dataset, dict):
        trs_dataset = list(trs_dataset.items())
        trs_dataset.sort()
        trs_dataset = [x[1] for x in trs_dataset]
    if dataset_id is None:
        with open(output_file, "w") as fh:
            if trs_dataset:
                fh.write(trs_dataset[0].get_header_string())
            for trf_obj in trs_dataset:
                data = str(trf_obj)
                fh.write(data)
    else:
        with open(output_file, "a") as fh:
            if trs_dataset:
                fh.write(trs_dataset[0].get_header_string())
            for trf_obj in trs_dataset:
                trf_obj.giid = dataset_id
                data = str(trf_obj)
                fh.write(data)


def save_trs_class_dataset(tr_class_dataset, output_file):
    """Save TRs classification objects to file."""
    if isinstance(tr_class_dataset, dict):
        tr_class_dataset = list(tr_class_dataset.items())
        tr_class_dataset.sort()
        tr_class_dataset = [x[1] for x in tr_class_dataset]
    with open(output_file, "w") as fh:
        for class_obj in tr_class_dataset:
            data = str(class_obj)
            fh.write(data)


def save_trs_as_fasta(
    trf_file, fasta_file, project, add_project=False, skip_alpha=False
):
    """Save TRs dataset as one fasta file."""
    trf_objs = []
    for trf_obj in sc_iter_tab_file(trf_file, TRModel):
        trf_objs.append(trf_obj)
    with open(fasta_file, "w") as fh_fasta:
        for trf_obj in trf_objs:
            if skip_alpha:
                if trf_obj.trf_family == "ALPHA":
                    continue
            fh_fasta.write(trf_obj.get_fasta_repr(add_project=add_project))


def get_classification_dict(fam_kmer_index_file):
    """ """
    kmer2fam = defaultdict(list)
    with open(fam_kmer_index_file) as fh:
        for line in fh:
            fam, k, rk, tf, df = line.strip().split("\t")
            df = float(df)
            kmer2fam[k].append((fam, df))
            kmer2fam[rk].append((fam, df))
    return kmer2fam
