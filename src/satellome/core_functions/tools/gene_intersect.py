#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.02.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

from collections import Counter
from intervaltree import IntervalTree
from satellome.core_functions.io.tr_file import save_trs_dataset
from satellome.core_functions.io.tab_file import sc_iter_tab_file
from satellome.core_functions.models.trf_model import TRModel
from satellome.core_functions.io.gff_file import sc_gff3_reader
from satellome.core_functions.tools.processing import count_lines_large_file
from tqdm import tqdm

def categorize_intervals(a, b, feature):
    if a[0] == b[0] and a[1] == b[1]:
        return f"TR equal {feature}"
    elif a[0] >= b[0] and a[1] <= b[1]:
        return f"TR in {feature}"
    elif b[0] >= a[0] and b[1] <= a[1]:
        return f"{feature} in TR"
    elif a[1] > b[0] and a[1] <= b[1] and a[0] < b[0]:
        return f"TR overlap left {feature}"
    elif a[0] < b[1] and a[0] >= b[0] and a[1] > b[1]:
        return f"TR overlap right {feature}"
    else:
        return "No overlap"
    
def interval_intersection(a, b):
    """
    Returns the intersection of two intervals.
    If no intersection, returns (0, 0).
    """
    return max(a[0], b[0]), min(a[1], b[1])

def interval_length(interval):
    """
    Returns the length of the interval.
    """
    return interval[1] - interval[0]


def filter_hits(hits):

    # print(hits)

    features = [x[2] for x in hits]

    hits = [x for x in hits if x[2] != "region"]

    if "CDS" in features:
        # assert features.count("CDS") == 1, features
        hits = [x for x in hits if x[2] not in ["gene", "mRNA", "exon"]]
    if "ncRNA" in features:
        # assert features.count("ncRNA") == 1, features
        # assert features.count("gene") == 1, features
        hits = [x for x in hits if x[2] not in ["gene", "exon"]]
    if "tRNA" in features:
        hits = [x for x in hits if x[2] not in ["gene", "exon"]]

    hits = list(set([tuple(x) for x in hits]))

    if len(hits) <= 1:
        return hits
    # else:
    #     print(hits)
    #     input("?")
    return hits


def _add_annotation(trf_file, gff_file, rm_file):
    ''' Add annotation to TRF file from GFF file.'''

    chrm2annotation = {}

    gff_file_lines = count_lines_large_file(gff_file)

    ### TODO: Make it less memory consuming
    for gff_record in tqdm(sc_gff3_reader(gff_file), total=gff_file_lines, desc="Load GFF"):
        chrm = gff_record.seqid
        if chrm not in chrm2annotation:
            chrm2annotation[chrm] = IntervalTree()
        chrm2annotation[chrm].addi(gff_record.start-1, gff_record.end, gff_record.as_dict())


    if rm_file:

        rm_file_lines = count_lines_large_file(rm_file)

        with open(rm_file) as fh:
            for line in tqdm(fh, total=rm_file_lines, desc="Load RM"):
                d = line.split()
                chrm = d[4]
                if chrm not in chrm2annotation:
                    chrm2annotation[chrm] = IntervalTree()
                start = int(d[5])
                end = int(d[6])
                fam = d[10]

                chrm2annotation[chrm].addi(start-1, end, {
                    "type": f"RM_{fam}",
                })

    trf_id2annotation = {}

    trf_file_lines = count_lines_large_file(trf_file)

    for j, trf_obj in tqdm(enumerate(sc_iter_tab_file(trf_file, TRModel)), total=trf_file_lines, desc="Load TRF"):

        chrm = trf_obj.trf_head.split()[0]
        start = trf_obj.trf_l_ind
        end = trf_obj.trf_r_ind

        trf_key = (trf_obj.trf_id, trf_obj.trf_l_ind, trf_obj.trf_r_ind)

        assert start < end, (start, end)

        if chrm not in chrm2annotation:
            print("No annotation for %s" % chrm)
            continue
        if chrm2annotation[chrm][start:end]:
            # print(trf_obj.as_dict())
            hits = [[x[0], x[1], x[-1]["type"]] for x in chrm2annotation[chrm][start:end]]
            for hit in hits:
                hit.append(categorize_intervals((start, end), (hit[0], hit[1]), hit[2]))
                
            hits = sorted(hits, key=lambda x: x[-1], reverse=True)

            hits = filter_hits(hits)
            
            trf_id2annotation[trf_key] = hits

        else:
            trf_id2annotation[trf_key] = None


    return trf_id2annotation


def add_annotation_from_gff(trf_file, gff_file, report_file, rm_file=None):
    ''' Add annotation to TRF file from GFF file.'''

    trf_id2annotation = _add_annotation(trf_file, gff_file, rm_file)

    c = Counter()
    for trf_id in trf_id2annotation:
        key = [x[-1] for x in trf_id2annotation[trf_id]]
        key.sort()
        key = tuple(key)
        c[key] += 1

    with open(report_file, "w") as fw:
        total = sum(c.values())
        for key, value in c.most_common(10000):
            if not key:
                key = ["No annotation"]
            fw.write(f"{'|'.join(key)}\t{value}\t{round(100.*value/total, 2)}%\n")

    all_trs = {}
    trf_ids = []
    for j, trf_obj in enumerate(sc_iter_tab_file(trf_file, TRModel)):
        all_trs[trf_obj.trf_id] = trf_obj
        trf_ids.append(trf_obj.trf_id)
    for trf_id, start, end in trf_id2annotation:
        key = [x[-1] for x in trf_id2annotation[(trf_id, start, end)]]
        if not key:
            key = ["No annotation"]
        key = '|'.join(key)
        all_trs[trf_id].trf_ref_annotation = key

    dataset = []
    for trf_id in trf_ids:
        dataset.append(all_trs[trf_id])
    save_trs_dataset(dataset, trf_file)


def get_gene_density(gff_file):

    chrm2annotation = {}

    gff_file_lines = count_lines_large_file(gff_file)

    for gff_record in tqdm(sc_gff3_reader(gff_file), total=gff_file_lines, desc="Load GFF"):
        chrm = gff_record.seqid
        if chrm not in chrm2annotation:
            chrm2annotation[chrm] = IntervalTree()
        if gff_record.type != "CDS":
            continue
        chrm2annotation[chrm].addi(gff_record.start-1, gff_record.end, "gene")
    
    
