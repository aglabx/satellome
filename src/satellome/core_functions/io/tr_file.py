#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 07.09.2011
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com
"""
Utilities for loading, saving, and managing tandem repeat (TR) datasets.

Provides convenience functions for working with tab-delimited TR files,
including batch loading, dictionary conversions, FASTA export, and
classification data management. Simplifies common TR analysis workflows.

Functions:
    read_trid2ngrams: Load n-gram annotations for repeats from folder
    read_trid2meta: Load TR ID to metadata string mapping
    get_all_trf_objs: Load all TRF objects from file into list
    get_all_class_objs: Load all classification objects into list
    get_class_objs_dict: Load classification objects into ID-keyed dictionary
    get_trf_objs_dict: Load TRF objects into ID-keyed dictionary
    get_trfid_obj_dict: Alias for get_trf_objs_dict
    save_trs_dataset: Save TR dataset to tab-delimited file
    save_trs_class_dataset: Save classification dataset to file
    save_trs_as_fasta: Export TR dataset to FASTA format
    get_classification_dict: Load k-mer to family classification index

Key Features:
    - Batch loading of TR datasets into memory
    - Dictionary-based lookups by TR ID
    - FASTA export with optional project metadata
    - Classification object management
    - N-gram annotation integration

Example:
    >>> # Load all repeats into dictionary
    >>> trid2obj = get_trf_objs_dict("repeats.tab")
    >>> trf = trid2obj[12345]
    >>> print(f"{trf.trf_consensus}: {trf.trf_l_ind}-{trf.trf_r_ind}")
    >>>
    >>> # Export to FASTA
    >>> save_trs_as_fasta("repeats.tab", "repeats.fasta", "hg38")
    >>>
    >>> # Load classification data
    >>> kmer2fam = get_classification_dict("kmer_index.tab")

See Also:
    satellome.core_functions.io.tab_file: Low-level tab file I/O
    satellome.core_functions.models.trf_model: TRModel and classification models
"""
import os
from collections import defaultdict

from satellome.core_functions.io.tab_file import (sc_iter_simple_tab_file,
                                                 sc_iter_tab_file)
from satellome.core_functions.models.trf_model import (TRModel,
                                                      TRsClassificationModel)


def read_trid2ngrams(annotation_ngram_folder, trf_large_file):
    """
    Load n-gram annotations for all tandem repeats from folder.

    Reads per-repeat n-gram files (*.ngram) from a folder and builds a
    dictionary mapping TR IDs to n-gram term frequency data.

    Args:
        annotation_ngram_folder (str): Path to folder containing {trf_id}.ngram files
        trf_large_file (str): Path to tab-delimited TRF file containing TR metadata

    Returns:
        dict: Maps trf_id (int) to list of (n-gram, tf) tuples, where:
             - n-gram (str): Forward or reverse n-gram sequence
             - tf (float): Term frequency for this n-gram
             Each TR has both forward and reverse n-grams in the list

    Example:
        >>> trid2ngrams = read_trid2ngrams("ngrams/", "repeats.tab")
        >>> ngrams = trid2ngrams[12345]
        >>> for ngram, freq in ngrams[:5]:
        ...     print(f"{ngram}: {freq:.3f}")

    Note:
        - Expected n-gram file format: {ngram}<tab>{rev_ngram}<tab>{tf} per line
        - Both forward and reverse n-grams are stored in the result
        - File naming: {trf_id}.ngram (e.g., "12345.ngram")
    """

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
    """
    Load TR ID to tab-delimited metadata string mapping.

    Reads TRF file and creates dictionary mapping each TR ID to its
    complete tab-delimited string representation.

    Args:
        file_name (str): Path to tab-delimited TRF file

    Returns:
        dict: Maps trf_id (int) to tab-delimited metadata string

    Example:
        >>> trid2meta = read_trid2meta("repeats.tab")
        >>> metadata_line = trid2meta[12345]
        >>> print(metadata_line)  # Full tab-delimited row

    Note:
        - Result strings include newlines (from str(trf_obj))
        - Useful for quick lookups without reparsing entire objects
    """
    trid2meta = {}
    for trf_obj in sc_iter_tab_file(file_name, TRModel):
        trid2meta[trf_obj.trf_id] = str(trf_obj)
    return trid2meta


def get_all_trf_objs(trf_large_file):
    """
    Load all tandem repeat objects from file into list.

    Reads entire TRF file and returns all TRModel objects in a list.
    Not memory-efficient for large datasets.

    Args:
        trf_large_file (str): Path to tab-delimited TRF file

    Returns:
        list: All TRModel objects from file

    Example:
        >>> trf_objs = get_all_trf_objs("repeats.tab")
        >>> print(f"Total repeats: {len(trf_objs)}")
        >>> for trf in trf_objs[:10]:
        ...     print(trf.trf_consensus)

    Note:
        - Loads entire file into memory
        - For large files or streaming, use sc_iter_tab_file() directly
        - For dictionary lookup, use get_trf_objs_dict()
    """
    result = []
    for trf_obj in sc_iter_tab_file(trf_large_file, TRModel):
        result.append(trf_obj)
    return result


def get_all_class_objs(trf_class_file):
    """
    Load all classification objects from file into list.

    Reads entire classification file and returns all TRsClassificationModel
    objects in a list. Not memory-efficient for large datasets.

    Args:
        trf_class_file (str): Path to tab-delimited classification file

    Returns:
        list: All TRsClassificationModel objects from file

    Example:
        >>> class_objs = get_all_class_objs("classifications.tab")
        >>> for obj in class_objs:
        ...     print(f"{obj.trf_id}: {obj.class_name} / {obj.family_name}")

    Note:
        - Loads entire file into memory
        - For dictionary lookup, use get_class_objs_dict()
    """
    result = []
    for trf_obj in sc_iter_tab_file(trf_class_file, TRsClassificationModel):
        result.append(trf_obj)
    return result


def get_class_objs_dict(trf_class_file):
    """
    Load classification objects into TR ID-keyed dictionary.

    Reads classification file and returns dictionary mapping TR IDs to
    TRsClassificationModel objects for fast lookups.

    Args:
        trf_class_file (str): Path to tab-delimited classification file

    Returns:
        dict: Maps trf_id (int) to TRsClassificationModel object

    Example:
        >>> trid2class = get_class_objs_dict("classifications.tab")
        >>> class_obj = trid2class[12345]
        >>> print(f"{class_obj.class_name} / {class_obj.family_name}")

    Note:
        - Loads entire file into memory as dictionary
        - Enables O(1) lookup by TR ID
    """
    result = {}
    for trf_obj in sc_iter_tab_file(trf_class_file, TRsClassificationModel):
        result[trf_obj.trf_id] = trf_obj
    return result


def get_trf_objs_dict(trf_large_file):
    """
    Load TRF objects into TR ID-keyed dictionary.

    Reads TRF file and returns dictionary mapping TR IDs to TRModel objects
    for fast O(1) lookups.

    Args:
        trf_large_file (str): Path to tab-delimited TRF file

    Returns:
        dict: Maps trf_id (int) to TRModel object

    Example:
        >>> trid2obj = get_trf_objs_dict("repeats.tab")
        >>> trf = trid2obj[12345]
        >>> print(f"{trf.trf_consensus}: {trf.trf_l_ind}-{trf.trf_r_ind}")

    Note:
        - Loads entire file into memory as dictionary
        - Enables fast lookups by TR ID
        - Alias: get_trfid_obj_dict()
    """
    result = {}
    for trf_obj in sc_iter_tab_file(trf_large_file, TRModel):
        result[trf_obj.trf_id] = trf_obj
    return result


def get_trfid_obj_dict(trf_large_file):
    """
    Load TRF objects into TR ID-keyed dictionary (alias).

    Alias for get_trf_objs_dict(). Loads TRF file into dictionary
    mapping TR IDs to TRModel objects.

    Args:
        trf_large_file (str): Path to tab-delimited TRF file

    Returns:
        dict: Maps trf_id (int) to TRModel object

    See Also:
        get_trf_objs_dict: Preferred function name with identical functionality

    Note:
        - This function exists for backward compatibility
        - Uses get_all_trf_objs() internally (less efficient than get_trf_objs_dict)
    """
    trs_dataset = get_all_trf_objs(trf_large_file)
    trid2obj = {}
    for trf_obj in trs_dataset:
        trid2obj[trf_obj.trf_id] = trf_obj
    return trid2obj


def save_trs_dataset(trs_dataset, output_file, dataset_id=None):
    """
    Save tandem repeat dataset to tab-delimited file.

    Writes TRModel objects to file with optional header and dataset ID
    assignment. Can accept either list or dictionary of TR objects.

    Args:
        trs_dataset (list or dict): TRModel objects to save
                                    If dict, sorted by key before writing
        output_file (str): Path to output tab-delimited file
        dataset_id (int, optional): If provided:
                                   - Appends to existing file (mode 'a')
                                   - Sets giid field to this value for all objects
                                   If None:
                                   - Creates new file (mode 'w')
                                   - Writes header line (from first object)
                                   Defaults to None.

    Example:
        >>> # Save new dataset with header
        >>> trf_objs = get_all_trf_objs("input.tab")
        >>> save_trs_dataset(trf_objs, "output.tab")
        >>>
        >>> # Append to existing file with dataset ID
        >>> more_objs = get_all_trf_objs("more_input.tab")
        >>> save_trs_dataset(more_objs, "output.tab", dataset_id=2)
        >>>
        >>> # Save dictionary (sorted by key)
        >>> trid2obj = get_trf_objs_dict("input.tab")
        >>> save_trs_dataset(trid2obj, "output.tab")

    Note:
        - Header written only if dataset_id is None
        - Dictionary input is sorted by key before writing (deterministic order)
        - All objects get giid field set if dataset_id provided
        - Uses TRModel.__str__() for serialization
    """
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
    """
    Save TR classification dataset to tab-delimited file.

    Writes TRsClassificationModel objects to file. Can accept either
    list or dictionary of classification objects.

    Args:
        tr_class_dataset (list or dict): TRsClassificationModel objects to save
                                         If dict, sorted by key before writing
        output_file (str): Path to output tab-delimited file

    Example:
        >>> class_objs = get_all_class_objs("classifications.tab")
        >>> save_trs_class_dataset(class_objs, "output_classifications.tab")
        >>>
        >>> # Save dictionary (sorted by key)
        >>> trid2class = get_class_objs_dict("classifications.tab")
        >>> save_trs_class_dataset(trid2class, "output_classifications.tab")

    Note:
        - Creates new file (mode 'w')
        - No header line written
        - Dictionary input is sorted by key before writing
        - Uses TRsClassificationModel.__str__() for serialization
    """
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
    """
    Export tandem repeat dataset to FASTA format.

    Reads TRF file and writes tandem repeat sequences to FASTA format
    with optional project metadata and alpha satellite filtering.

    Args:
        trf_file (str): Path to input tab-delimited TRF file
        fasta_file (str): Path to output FASTA file
        project (str): Project name (currently unused in implementation)
        add_project (bool, optional): Add project metadata to FASTA headers.
                                     Defaults to False.
        skip_alpha (bool, optional): Skip repeats with trf_family == "ALPHA".
                                    Defaults to False.

    Example:
        >>> # Basic export
        >>> save_trs_as_fasta("repeats.tab", "repeats.fasta", "hg38")
        >>>
        >>> # With project metadata, skip alpha satellites
        >>> save_trs_as_fasta(
        ...     "repeats.tab",
        ...     "repeats.fasta",
        ...     "hg38",
        ...     add_project=True,
        ...     skip_alpha=True
        ... )

    Note:
        - Uses TRModel.get_fasta_repr() for sequence formatting
        - FASTA headers include TR ID and genomic coordinates
        - Alpha satellite filtering checks trf_family attribute
        - Loads all objects into memory before writing
    """
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
    """
    Load k-mer to family classification index from file.

    Reads tab-delimited k-mer index file and builds dictionary mapping
    k-mers (forward and reverse) to lists of (family, df_score) tuples
    for classification purposes.

    Args:
        fam_kmer_index_file (str): Path to k-mer index file
                                   Format: family<tab>kmer<tab>rev_kmer<tab>tf<tab>df

    Returns:
        defaultdict: Maps k-mer (str) to list of (family, df) tuples where:
                    - family (str): TR family name
                    - df (float): Document frequency score for this k-mer/family
                    Both forward and reverse k-mers are keys

    Example:
        >>> kmer2fam = get_classification_dict("kmer_family_index.tab")
        >>> families = kmer2fam["AATAT"]
        >>> for family, score in families:
        ...     print(f"{family}: {score:.3f}")
        SAT1: 0.850
        SAT2: 0.120

    Note:
        - Input file format: 5 tab-separated columns per line
        - Column order: family, k-mer, reverse k-mer, tf, df
        - Both forward and reverse k-mers added as keys
        - Returns defaultdict(list) - missing keys return empty list
        - TF (term frequency) column is read but not stored
    """
    kmer2fam = defaultdict(list)
    with open(fam_kmer_index_file) as fh:
        for line in fh:
            fam, k, rk, tf, df = line.strip().split("\t")
            df = float(df)
            kmer2fam[k].append((fam, df))
            kmer2fam[rk].append((fam, df))
    return kmer2fam
