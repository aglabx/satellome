#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 11.01.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import plotly.graph_objects as go
import plotly.express as px
from trf_drawing import (
    scaffold_length_sort_length,
    read_trf_file,
    get_gaps_annotation,
)
import os
from collections import Counter
import pandas as pd
from intervaltree import IntervalTree
from trf_embedings import get_disances


class Graph:

    # init function to declare class variables
    def __init__(self, V):
        self.V = []
        self.id2node = {}
        self.node2id = {}
        for i, id1 in enumerate(V):
            self.id2node[id1] = i
            self.node2id[i] = id1
            self.V.append(i)
        self.adj = [[] for i in V]

    def DFSUtil(self, temp, v, visited):

        # Mark the current vertex as visited
        visited[v] = True

        # Store the vertex to list
        temp.append(v)

        # Repeat for all vertices adjacent
        # to this vertex v
        for i, dist in self.adj[v]:
            if visited[i] == False:
                # Update the list
                temp = self.DFSUtil(temp, i, visited)
        return temp

    # method to add an undirected edge
    def addEdge(self, v, w, dist):
        id1 = self.id2node[v]
        id2 = self.id2node[w]
        self.adj[id1].append((id2, dist))
        self.adj[id2].append((id1, dist))

    def remove_edges_by_distances(self, cutoff):
        for id1 in self.V:
            new_adj = []
            for id2, dist in self.adj[id1]:
                if dist < cutoff:
                    new_adj.append((id2, dist))
            self.adj[id1] = new_adj

    # Method to retrieve connected components
    # in an undirected graph
    def connectedComponents(self):
        visited = []
        cc = []
        for i in self.V:
            visited.append(False)
        for v in self.V:
            if visited[v] == False:
                temp = []
                cc.append(self.DFSUtil(temp, v, visited))
        return cc


def name_clusters(distances, tr2vector, df_trs, level=1):

    all_distances = list(distances.values())
    all_distances = list(set(map(int, all_distances)))
    all_distances.sort(reverse=True)
    start_cutoff = min(int(all_distances[0]), 50)

    G = Graph(list(tr2vector.keys()))
    for (id1, id2) in distances:
        G.addEdge(id1, id2, distances[(id1, id2)])

    for i in range(start_cutoff, level - 1, -1):
        G.remove_edges_by_distances(i)
        comps = G.connectedComponents()
        items = []
        singl = []
        for c in comps:
            ids = [(G.node2id[id1], df_trs.loc[G.node2id[id1]].period) for id1 in c]
            if len(ids) > 3:  # Why 3? It should be 1
                items.append(ids)
            else:
                singl += ids
        items.sort(key=lambda x: len(x))
        print(
            i, "->", len(comps), len(singl)
        )  # print(f"For distance {i} -> total number of components is {len(comps)} and {len(singl)} out of them are singletones")
        for class_name, d in enumerate(items):
            median_monomer = [x[1] for x in d]
            median_monomer.sort()
            median_monomer = median_monomer[int(len(median_monomer) / 2)]

            name = f"{class_name}_{median_monomer}"
            for id1, period in d:
                df_trs.at[id1, "family_name"] = name
                df_trs.at[id1, "locus_name"] = name

        for id1, period in singl:
            df_trs.at[id1, "family_name"] = "SING"

    return df_trs, tr2vector, distances, all_distances


def _draw_sankey(output_file_name, title_text, labels, source, target, value):
    fig = go.Figure(
        data=[
            go.Sankey(
                node=dict(
                    pad=15,
                    thickness=20,
                    line=dict(color="black", width=1),
                    label=labels,
                    color="blue",
                ),
                link=dict(
                    source=source,
                    target=target,
                    value=value,
                ),
            )
        ]
    )
    fig.update_layout(
        title_text=title_text,
        font_size=10,
        height=2000,
        width=2000,
    )
    fig.write_image(output_file_name)


def draw_sankey(
    output_file_name,
    title_text,
    df_trs,
    tr2vector,
    distances,
    all_distances,
    skip_singletons=True,
):

    G = Graph(list(tr2vector.keys()))
    for (id1, id2) in distances:
        G.addEdge(id1, id2, distances[(id1, id2)])

    steps = []

    start_cutoff = int(all_distances[0])

    name2trs = {}
    last_n_comp = 0
    id2names = {}

    for i in range(start_cutoff, 0, -1):
        G.remove_edges_by_distances(i)
        comps = G.connectedComponents()
        print(i, "->", len(comps))
        if len(comps) == last_n_comp:
            print("..skipped")
            continue
        last_n_comp = len(comps)

        items = []
        singl = []

        id2INstep = {}
        name2size = {}
        name2ids = {}
        name2id = {}

        for c in comps:
            ids = [(G.node2id[id1], df_trs.loc[G.node2id[id1]].period) for id1 in c]
            if len(ids) > 3:
                items.append(ids)
            else:
                singl += ids

        items.sort(key=lambda x: len(x))
        for class_name, d in enumerate(items):
            median_monomer = [x[1] for x in d]
            median_monomer.sort()
            median_monomer = median_monomer[int(len(median_monomer) / 2)]
            name = f"{i}_{class_name}_{median_monomer}"
            for id1, period in d:
                id2INstep[id1] = name
                name2id[name] = id1
            name2size[name] = len(d)
            name2ids[name] = d
            name2trs[name] = d

            for id_, _ in d:
                id2names.setdefault(id_, [])
                id2names[id_].append(name)

        if not skip_singletons and singl:
            name = f"{i}_SING"
            for id1, period in singl:
                id2INstep[id1] = name
                name2id[name] = id1

                id2names.setdefault(id1, [])
                id2names[id1].append(name)

            name2size[name] = len(singl)
            name2ids[name] = singl

        steps.append((id2INstep, name2size, name2ids, name2id))

    labels = []
    source = []
    target = []
    value = []
    name2monomers = {}
    name2lid = {}
    lid = 0
    prev_id2INstep, prev_name2size, name2ids, prev_name2id = steps[0]
    for name in prev_name2size:
        labels.append(name)
        name2lid[name] = lid
        lid += 1

        name2monomers[name] = name2ids[name]

    for id2INstep, name2size, name2ids, name2id in steps[1:]:

        print(name2size)

        for name in name2size:
            labels.append(name)
            name2lid[name] = lid
            lid += 1

            start = name2lid[prev_id2INstep[name2id[name]]]
            end = name2lid[name]

            source.append(start)
            target.append(end)
            value.append(name2size[name])

            name2monomers[name] = name2ids[name]

        prev_id2INstep = id2INstep

    _draw_sankey(output_file_name, title_text, labels, source, target, value)

    return name2monomers, name2lid, name2trs, id2names


def draw_spheres(output_file_name_prefix, title_text, df_trs):

    fig = px.scatter_3d(
        df_trs,
        x="gc",
        y="period",
        z="pmatch",
        color="family_name",
        size="log_length",
    )
    fig.update_layout(
        title={
            "text": title_text,
            "y": 0.99,
            "x": 0.5,
            "xanchor": "center",
            "yanchor": "top",
        }
    )
    fig.update_layout(width=800, height=800)
    output_file_name = output_file_name_prefix + ".3D.png"
    fig.write_image(output_file_name)

    fig = px.scatter_3d(
        df_trs[df_trs["family_name"] != "SING"],
        x="gc",
        y="period",
        z="pmatch",
        color="family_name",
        size="log_length",
    )
    fig.update_layout(
        title={
            "text": title_text + " No Singletons",
            "y": 0.99,
            "x": 0.5,
            "xanchor": "center",
            "yanchor": "top",
        }
    )
    fig.update_layout(width=800, height=800)
    output_file_name = output_file_name_prefix + ".3D.nosingl.png"
    fig.write_image(output_file_name)

    fig = px.scatter(df_trs, x="gc", y="period", color="family_name", size="log_length")
    output_file_name = output_file_name_prefix + ".2D.gc_period.png"
    fig.write_image(output_file_name)

    fig = px.scatter(df_trs, x="gc", y="period", color="family_name", size="log_length")
    output_file_name = output_file_name_prefix + ".2D.gc_period.png"
    fig.write_image(output_file_name)

    fig = px.scatter(df_trs, x="gc", y="pmatch", color="family_name", size="log_length")
    output_file_name = output_file_name_prefix + ".2D.gc_pmatch.png"
    fig.write_image(output_file_name)

    fig = px.scatter(
        df_trs, x="pmatch", y="period", color="family_name", size="log_length"
    )
    output_file_name = output_file_name_prefix + ".2D.period_period.png"
    fig.write_image(output_file_name)


def _draw_chromosomes(scaffold_for_plot, title_text, use_chrm=False):

    if use_chrm:
        scaffold_items = scaffold_for_plot["chrm"]
        yaxis_title = "Chromosome name"
    else:
        scaffold_items = scaffold_for_plot["scaffold"]
        yaxis_title = "Scaffold name"

    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=scaffold_for_plot["end"],
            y=scaffold_items,
            orientation="h",
            name="Scaffold",
            marker_color="#f3f4f7",
        )
    )
    fig.update_layout(barmode="overlay")
    fig.update_layout(
        title={
            "text": title_text,
            "y": 0.99,
            "x": 0.5,
            "xanchor": "center",
            "yanchor": "top",
        },
        xaxis_title="bp",
        yaxis_title=yaxis_title,
    )

    fig.update_layout(
        xaxis=dict(
            automargin=True,
            showline=True,
            showgrid=False,
            showticklabels=True,
            linecolor="rgb(204, 204, 204)",
            linewidth=1,
            ticks="outside",
            rangemode="nonnegative",
            tickfont=dict(
                family="Arial",
                size=15,
                color="rgb(82, 82, 82)",
            ),
        ),
        # Turn off everything on y axis
        yaxis=dict(
            showgrid=False,
            zeroline=False,
            showline=False,
            showticklabels=True,
            ticklabelstep=1,
            tickwidth=15,
            tickfont=dict(
                family="Arial",
                size=15,
                color="rgb(82, 82, 82)",
            ),
        ),
        width=1400,
        height=900,
        margin=dict(
            autoexpand=True,
            # l=150,
            # r=20,
            # t=110,
        ),
        showlegend=True,
        plot_bgcolor="white",
    )

    fig.update_layout(legend=dict(font=dict(family="Arial", size=15, color="black")))

    fig.update_xaxes(range=[0, max(scaffold_for_plot["end"]) + 1000])

    return fig


def draw_karyotypes(
    output_file_name_prefix,
    title_text,
    df_trs,
    scaffold_for_plot,
    gaps_df,
    repeats_with_gap,
    repeats_without_gaps,
    use_chrm=False,
    enhance=2000000,
    gap_cutoff=1000,
):

    if use_chrm:
        _df_trs = df_trs[df_trs["chrm"].isin(scaffold_for_plot["chrm"])]
    else:
        _df_trs = df_trs[df_trs["chrm"].isin(scaffold_for_plot["scaffold"])]

    ### 1. Raw gaps
    title_text_ = title_text + "(raw gaps)"
    fig = _draw_chromosomes(scaffold_for_plot, title_text_, use_chrm=use_chrm)
    fig.add_trace(
        go.Bar(
            base=gaps_df["start"],
            x=gaps_df["length"],
            y=gaps_df["scaffold"],
            orientation="h",
            name="gaps",
            marker_color="rgba(0, 0, 0)",
        )
    )
    output_file_name = output_file_name_prefix + ".gaps.png"
    fig.write_image(output_file_name)

    ### 2. Enhanced gaps
    title_text_ = title_text + "(enlarged gaps)"
    fig = _draw_chromosomes(scaffold_for_plot, title_text_, use_chrm=use_chrm)
    _gaps_df = gaps_df[gaps_df["length"] > gap_cutoff]
    _gaps_df["length"] = [max(x["length"], enhance) for i, x in _gaps_df.iterrows()]
    fig.add_trace(
        go.Bar(
            base=_gaps_df["start"],
            x=_gaps_df["length"],
            y=_gaps_df["scaffold"],
            orientation="h",
            name="gaps",
            marker_color="rgba(0, 0, 0)",
        )
    )
    output_file_name = output_file_name_prefix + f".gaps.{gap_cutoff}bp.enhanced.png"
    fig.write_image(output_file_name)

    ### #3. Enhanced repeats_with_gap
    title_text_ = title_text + "(enlarged TRs with gaps)"
    size = enhance
    fig = _draw_chromosomes(scaffold_for_plot, title_text_, use_chrm=use_chrm)
    repeats_with_gap_df = pd.DataFrame(
        repeats_with_gap,
        columns=["chrm", "start", "end", "family_name", "gap_type", "length"],
    )
    aN = repeats_with_gap_df.loc[
        (repeats_with_gap_df.gap_type == "aN") & (repeats_with_gap_df.length < size)
    ]
    Na = repeats_with_gap_df.loc[
        (repeats_with_gap_df.gap_type == "Na") & (repeats_with_gap_df.length < size)
    ]
    aNa = repeats_with_gap_df.loc[
        (repeats_with_gap_df.gap_type == "aNa") & (repeats_with_gap_df.length < size)
    ]

    fig.add_trace(
        go.Bar(
            base=aN["start"],
            x=[size] * len(aN),
            y=aN["chrm"],
            orientation="h",
            name="Tandem Repeat_with gap aN",
            marker_color="#FF00FF",
        )
    )

    fig.add_trace(
        go.Bar(
            base=aN["start"] + size,
            x=[size] * len(aN),
            y=aN["chrm"],
            orientation="h",
            name="gaps aN",
            marker_color="#663399",
        )
    )

    fig.add_trace(
        go.Bar(
            base=Na["start"],
            x=[size] * len(Na),
            y=Na["chrm"],
            orientation="h",
            name="Tandem Repeat_with gap Na",
            marker_color="#00CED1",
        )
    )

    fig.add_trace(
        go.Bar(
            base=Na["start"] + size,
            x=[size] * len(Na),
            y=Na["chrm"],
            orientation="h",
            name="gaps Na",
            marker_color="#00BFFF",
        )
    )

    fig.add_trace(
        go.Bar(
            base=aNa["start"],
            x=[size * 2 / 3] * len(aNa),
            y=aNa["chrm"],
            orientation="h",
            name="Tandem Repeat_with gap aNa",
            marker_color="#00FF7F",
        )
    )

    fig.add_trace(
        go.Bar(
            base=aNa["start"] + (size * 2 / 3),
            x=[size * 2 / 3] * len(aNa),
            y=aNa["chrm"],
            orientation="h",
            name="gaps aNa",
            marker_color="#228B22",
        )
    )

    fig.add_trace(
        go.Bar(
            base=aNa["start"] + +size * 2 / 3 + +size * 2 / 3,
            x=[size * 2 / 3] * len(aNa),
            y=aNa["chrm"],
            orientation="h",
            marker_color="#00FF7F",
        )
    )
    output_file_name = output_file_name_prefix + ".repeats.with.gaps.enhanced.png"
    fig.write_image(output_file_name)

    ### #4. Enhanced TRs without gaps
    _title_text = title_text + " (TRs without gaps)"
    fig = _draw_chromosomes(scaffold_for_plot, _title_text, use_chrm=use_chrm)
    names = set([x["family_name"] for x in repeats_without_gaps])
    for name in names:
        items = [x for x in repeats_without_gaps if x["family_name"] == name]
        fig.add_trace(
            go.Bar(
                base=[x.start for x in repeats_without_gaps],
                x=[max(x.length, size) for x in repeats_without_gaps],
                y=[x.chrm for x in repeats_without_gaps],
                orientation="h",
                name=name,
            )
        )
    output_file_name = output_file_name_prefix + ".repeats.nogaps.enhanced.png"
    fig.write_image(output_file_name)

    ### #5. Raw all TRs

    _title_text = title_text + " (all)"
    fig = _draw_chromosomes(scaffold_for_plot, _title_text, use_chrm=use_chrm)
    names = set([x["family_name"] for i, x in _df_trs.iterrows()])
    for name in names:
        items = _df_trs[_df_trs["family_name"] == name]
        fig.add_trace(
            go.Bar(
                base=items["start"],
                x=items["length"],
                y=items["chrm"],
                orientation="h",
                name=name,
            )
        )
    output_file_name = output_file_name_prefix + ".raw.png"
    fig.write_image(output_file_name)

    ### #6. Raw all TRs no singletons

    _title_text = title_text + " (no singletons)"
    fig = _draw_chromosomes(scaffold_for_plot, _title_text, use_chrm=use_chrm)
    names = set([x["family_name"] for i, x in _df_trs.iterrows()])
    for name in names:
        if name == "SING":
            continue
        items = _df_trs[_df_trs["family_name"] == name]
        fig.add_trace(
            go.Bar(
                base=items["start"],
                x=items["length"],
                y=items["chrm"],
                orientation="h",
                name=name,
            )
        )
    output_file_name = output_file_name_prefix + ".nosing.png"
    fig.write_image(output_file_name)

    ### #7. Raw TRs singletons

    _title_text = title_text + " (no singletons)"
    fig = _draw_chromosomes(scaffold_for_plot, _title_text, use_chrm=use_chrm)
    names = set([x["family_name"] for i, x in _df_trs.iterrows()])
    for name in names:
        if name != "SING":
            continue
        items = _df_trs[_df_trs["family_name"] == name]
        fig.add_trace(
            go.Bar(
                base=items["start"],
                x=items["length"],
                y=items["chrm"],
                orientation="h",
                name=name,
            )
        )
    output_file_name = output_file_name_prefix + ".sing.png"
    fig.write_image(output_file_name)

    ### #8. Raw all TRs (enhanced)

    _title_text = title_text + " (enlarged, all)"
    fig = _draw_chromosomes(scaffold_for_plot, _title_text, use_chrm=use_chrm)
    names = set([x["family_name"] for i, x in _df_trs.iterrows()])
    for name in names:
        items = _df_trs[_df_trs["family_name"] == name]
        items["length"] = [max(x["length"], enhance) for i, x in items.iterrows()]
        fig.add_trace(
            go.Bar(
                base=items["start"],
                x=items["length"],
                y=items["chrm"],
                orientation="h",
                name=name,
            )
        )
    output_file_name = output_file_name_prefix + ".raw.enhanced.png"
    fig.write_image(output_file_name)

    ### #9. Raw all TRs no singletons (enhanced)

    _title_text = title_text + " (enlarged, no singletons)"
    fig = _draw_chromosomes(scaffold_for_plot, _title_text, use_chrm=use_chrm)
    names = set([x["family_name"] for i, x in _df_trs.iterrows()])
    for name in names:
        if name == "SING":
            continue
        items = _df_trs[_df_trs["family_name"] == name]
        items["length"] = [max(x["length"], enhance) for i, x in items.iterrows()]
        fig.add_trace(
            go.Bar(
                base=items["start"],
                x=items["length"],
                y=items["chrm"],
                orientation="h",
                name=name,
            )
        )
    output_file_name = output_file_name_prefix + ".nosing.enchanced.png"
    fig.write_image(output_file_name)

    ### #10. Raw TRs singletons (enhanced)

    _title_text = title_text + " (enlarged, only singletons)"
    fig = _draw_chromosomes(scaffold_for_plot, _title_text, use_chrm=use_chrm)
    names = set([x["family_name"] for i, x in _df_trs.iterrows()])
    for name in names:
        if name != "SING":
            continue
        items = _df_trs[_df_trs["family_name"] == name]
        items["length"] = [max(x["length"], enhance) for i, x in items.iterrows()]
        fig.add_trace(
            go.Bar(
                base=items["start"],
                x=items["length"],
                y=items["chrm"],
                orientation="h",
                name=name,
            )
        )
    output_file_name = output_file_name_prefix + ".sing.enchanced.png"
    fig.write_image(output_file_name)


def draw_all(
    trf_file,
    fasta_file,
    chm2name,
    output_folder,
    taxon,
    lenght_cutoff=10000000,
    level=1,
    enhance=1000000,
    gap_cutoff=1000,
):

    print("Loading chromosomes...")
    scaffold_df = scaffold_length_sort_length(fasta_file, lenght_cutoff=lenght_cutoff)

    print("Loading trs...")
    df_trs = read_trf_file(trf_file)
    df_trs = df_trs.loc[df_trs.period > 5]
    print(f"Quantity of TRs: {len(df_trs)}")

    if len(df_trs) > 1999:
        print("To many TRs")
        print("Filtering them...")
        df_trs = df_trs.sort_values(
            by="length", axis=0, ascending=False, ignore_index=True
        )[:2000]
        print(f"Updated quantity of TRs: 2000")

    print("Computing distances...")
    distances, tr2vector = get_disances(df_trs)
    print("Naming repeats...")
    df_trs, tr2vector, distances, all_distances = name_clusters(
        distances, tr2vector, df_trs, level=level
    )

    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)

    ### TODO: save distances and tr2vector

    output_file_name = os.path.join(output_folder, f"{taxon}.trs_flow.png")
    title_text = f"Tandem repeats flow in {taxon}"
    print("Draw sankey...")
    name2monomers, name2lid, name2ids, id2names = draw_sankey(
        output_file_name,
        title_text,
        df_trs,
        tr2vector,
        distances,
        all_distances,
        skip_singletons=True,
    )

    print("Draw spheres...")
    output_file_name_prefix = os.path.join(output_folder, f"{taxon}.spheres")
    title_text = f"Tandem repeats distribution in {taxon}"
    draw_spheres(output_file_name_prefix, title_text, df_trs)

    print("Draw gaps...")
    gaps_data = get_gaps_annotation(fasta_file, lenght_cutoff=lenght_cutoff)
    gaps_lengths = Counter([x[-1] for x in gaps_data])

    print("Gaps distribution:", gaps_lengths)

    gaps_df = pd.DataFrame(gaps_data, columns=["scaffold", "start", "end", "length"])

    df_trs["chrm"] = [x["scaffold"] for i, x in df_trs.iterrows()]

    chrms = set([x[0] for x in gaps_data])
    print(chrms)

    chrm2gapIT = {}
    for chrm in chrms:
        chrm2gapIT[chrm] = IntervalTree()

    for chrm, start, end, length in gaps_data:
        chrm2gapIT[chrm].addi(start, end, "gap")

    repeats_with_gap = []
    repeats_without_gaps = []
    for i, d in df_trs.iterrows():
        if d.chrm not in chrm2gapIT:
            continue
        found_gaps = chrm2gapIT[d.chrm][d.start - 10000 : d.end + 10000]
        if not found_gaps:
            repeats_without_gaps.append(d)
            continue
        found_gaps = list(found_gaps)
        start_gap = min([x.begin for x in found_gaps])
        end_gap = max([x.end for x in found_gaps])

        if d.start < start_gap and d.end < end_gap:
            gap_type = "aN"
        elif d.start < start_gap and d.end > end_gap:
            gap_type = "aNa"
        elif start_gap < d.start:
            gap_type = "Na"
        else:
            print("Unknown")
            print(d.family_name, d.start, d.end, start_gap, end_gap)
        repeats_with_gap.append(
            [
                d.chrm,
                min(d.start, start_gap),
                max(d.end, end_gap),
                d.family_name,
                gap_type,
                abs(min(d.start, start_gap) - max(d.end, end_gap)),
            ]
        )

    ### TODO: save gaps

    print("Draw karyotypes...")
    output_file_name_prefix = os.path.join(output_folder, f"{taxon}.karyo")
    title_text = f"Tandem repeats in {taxon}"
    draw_karyotypes(
        output_file_name_prefix,
        taxon,
        df_trs,
        scaffold_df,
        gaps_df,
        repeats_with_gap,
        repeats_without_gaps,
        use_chrm=chm2name,
        enhance=enhance,
        gap_cutoff=gap_cutoff,
    )
