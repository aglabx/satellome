#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# @created: 20.03.2023
# @author: Marina Popova
# @contact: marinaalexpopova@yandex.ru


import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from collections import defaultdict


def cluster(df, treshold=0.1, out_file="cluster_names.tsv"):
    cluster_df = pd.read_csv(df, sep="\t", names=['TR', 'TR length', 'period', range(0,1024)])
    X = cluster_df[:, 3:]
    clustering = AgglomerativeClustering(
        distance_threshold=treshold,
        metric="cosine",
        linkage="complete",
        n_clusters=None,
    ).fit(X)
    cluster_df["label"] = clustering.labels_
    return cluster_df


def naming(df, latin):
    df = df
    single = df.label.value_counts()[df.label.value_counts() == 1].index
    name = f"{latin.split()[0][:1]}{latin.split()[1][:3]}"
    family_median = defaultdict(list)
    for i in df.label:
        family_median[i] = list(df.loc[df.label == i].period)
    for i in family_median.keys():
        family_median[i].sort()
        family_median[i] = family_median[i][round(int(len(family_median[i])) / 2)]
    df["family"] = ["SING" if i in single else f"{name}_{i}_{}" for i in df.label]
    return df

