#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 11.01.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

from trseeker.tools.sequence_tools import get_revcomp
import math
import plotly.graph_objects as go
import plotly.express as px
from satelome.trf_drawing import scaffold_length_sort_length, read_trf_file, get_gaps_annotation
import os
from collections import Counter
import pandas as pd

class Graph:
 
    # init function to declare class variables
    def __init__(self, V):
        self.V = []
        self.id2node = {}
        self.node2id = {}
        for i,id1 in enumerate(V):
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


def get_pentatokens():
    token2id = {}
    token2revtoken = {}
    i = 0
    for n1 in "ACGT":
        for n2 in "ACGT":
            for n3 in "ACGT":
                for n4 in "ACGT":
                    for n5 in "ACGT":
                        token = "".join([n1,n2,n3,n4,n5])
                        token2id[token] = i
                        i += 1
                        token2revtoken[token] = get_revcomp(token)
    return token2id, token2revtoken

def fill_vectors(df_trs, token2id, token2revtoken, k=5):
    tr2vector = {}
    for id1, x in df_trs.iterrows():
        vector = [0] * len(token2id)
        seq = x['array'].upper()
        N = len(seq)-k+1
        for i in range(N):
            token = seq[i:i+k]
            if "N" in token:
                continue
            vector[token2id[token]] += 1
            vector[token2id[token2revtoken[token]]] += 1
        vector = [100.*x/(2*N) for x in vector]
        tr2vector[id1] = vector
    return tr2vector

def compute_distances(tr2vector):
    distances = {}
    for id1 in tr2vector:
        for id2 in tr2vector:
            distances[(id1,id2)] = math.dist(tr2vector[id1], tr2vector[id2])
    return distances

def name_clusters(df_trs, level=1):
    token2id, token2revtoken = get_pentatokens()
    tr2vector = fill_vectors(df_trs, token2id, token2revtoken, k=5)
    distances = compute_distances(tr2vector)

    all_distances = list(distances.values())
    all_distances = list(set(map(int, all_distances)))
    all_distances.sort(reverse=True)

    G = Graph(list(tr2vector.keys()))
    for (id1,id2) in distances:
        G.addEdge(id1, id2, distances[(id1,id2)])

    for i in range(level, level-1, -1):
        G.remove_edges_by_distances(i)
        comps = G.connectedComponents()
        print(i, "->", len(comps))
        items = []
        singl = []
        for c in comps:
            ids = [(G.node2id[id1], df_trs.loc[G.node2id[id1]].period) for id1 in c]
            if len(ids) > 3:
                items.append(ids)
            else:
                singl += ids
        items.sort(key=lambda x: len(x))
        print("Singletons", len(singl))
        for class_name, d in enumerate(items):
            print(len(d), [x[1] for x in d])
            
            median_monomer = [x[1] for x in d]
            median_monomer.sort()
            median_monomer = median_monomer[int(len(median_monomer)/2)]
            
            name = f"{class_name}_{median_monomer}" 
            for id1, period in d:
                df_trs.at[id1, 'family_name'] = name
        for id1, period in singl:
            df_trs.at[id1, 'family_name'] = "SING"

    return df_trs, tr2vector, distances, all_distances
        

def _draw_sankey(output_file_name, title_text, labels, source, target, value):
    fig = go.Figure(data=[go.Sankey(
        node = dict(
        pad = 15,
        thickness = 20,
        line = dict(color = "black", width = 1),
        label = labels,
        color = "blue",
        
        ),
        link = dict(
        source = source,
        target = target,
        value = value,
    ))])
    fig.update_layout(
            title_text=title_text, 
            font_size=10,
            height=2000,
            width=2000,
    )
    fig.write_image(output_file_name)


def draw_sankey(output_file_name, title_text, df_trs, tr2vector, distances, all_distances, skip_singletons=True):

    G = Graph(list(tr2vector.keys()))
    for (id1,id2) in distances:
        G.addEdge(id1, id2, distances[(id1,id2)])

    steps = []

    start_cutoff = int(all_distances[0])
    
    name2trs = {}
    last_n_comp = 0
    
    
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
            median_monomer = median_monomer[int(len(median_monomer)/2)]
            name = f"{i}_{class_name}_{median_monomer}" 
            for id1, period in d:
                id2INstep[id1] = name
                name2id[name] = id1
            name2size[name] = len(d)
            name2ids[name] = d
            name2trs[name] = d
            
            
        if not skip_singletons and singl:
            name = f"{i}_SING" 
            for id1, period in singl:
                id2INstep[id1] = name
                name2id[name] = id1
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
    
    return name2monomers, name2lid, name2trs


def draw_spheres(output_file_name_prefix, title_text, df_trs):
    
    fig = px.scatter_3d(df_trs, 
                        x='gc', 
                        y='period', 
                        z='pmatch', 
                        color='family_name', 
                        size='log_length',
                    )
    fig.update_layout(title={
            'text': title_text,
            'y':0.99,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top'})
    fig.update_layout(width=800, height=800)
    output_file_name = output_file_name_prefix + ".3D.png"
    fig.write_image(output_file_name)


    fig = px.scatter_3d(df_trs[df_trs["family_name"] != 'SING'], 
                        x='gc', 
                        y='period', 
                        z='pmatch', 
                        color='family_name', 
                        size='log_length',
                    )
    fig.update_layout(title={
            'text': title_text + " No Singletons",
            'y':0.99,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top'})
    fig.update_layout(width=800, height=800)
    output_file_name = output_file_name_prefix + ".3D.nosingl.png"
    fig.write_image(output_file_name) 

    fig = px.scatter(df_trs, 
                 x="gc", 
                 y="period", 
                 color="family_name",
                 size='log_length'
                )
    output_file_name = output_file_name_prefix + ".2D.gc_period.png"
    fig.write_image(output_file_name) 

    fig = px.scatter(df_trs, 
                 x="gc", 
                 y="period", 
                 color="family_name",
                 size='log_length'
                )
    output_file_name = output_file_name_prefix + ".2D.gc_period.png"
    fig.write_image(output_file_name)

    fig = px.scatter(df_trs, 
                 x="gc", 
                 y="pmatch", 
                 color="family_name",
                 size='log_length'
                )
    output_file_name = output_file_name_prefix + ".2D.gc_pmatch.png"
    fig.write_image(output_file_name)

    fig = px.scatter(df_trs, 
                 x="pmatch", 
                 y="period", 
                 color="family_name",
                 size='log_length'
                )
    output_file_name = output_file_name_prefix + ".2D.period_period.png"
    fig.write_image(output_file_name)


def _draw_chromosomes(scaffold_for_plot, title_text, use_chrm=False):

    if use_chrm:
        scaffold_items = scaffold_for_plot['chrm']
        yaxis_title='Chromosome name'
    else:
        scaffold_items = scaffold_for_plot['scaffold']
        yaxis_title='Scaffold name'
    
    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=scaffold_for_plot['end'],
        y=scaffold_items,
        orientation='h',
        name = 'Scaffold',
        marker_color='#f3f4f7',
    ))
    fig.update_layout(barmode='overlay')
    fig.update_layout(title={
            'text': title_text,
            'y':0.99,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top'}, 
            xaxis_title='bp', yaxis_title=yaxis_title)

    fig.update_layout(
        xaxis=dict(
            automargin=True,
            showline=True,
            showgrid=False,
            showticklabels=True,
            linecolor='rgb(204, 204, 204)',
            linewidth=1,
            ticks='outside',
            rangemode="nonnegative",
            tickfont=dict(
                family='Arial',
                size=15,
                color='rgb(82, 82, 82)',
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
                family='Arial',
                size=15,
                color='rgb(82, 82, 82)',
            ),
        ),
        width=1400,
        height=900,
        margin=dict(
            autoexpand=True,
            #l=150,
            #r=20,
            #t=110,
        ),
        showlegend=True,
        plot_bgcolor='white'
    )

    fig.update_layout(legend = dict(font = dict(family = "Arial", size = 15, color = "black")))

    fig.update_xaxes(range=[0, max(scaffold_for_plot['end'])+1000])

    return fig

def draw_karyotypes(output_file_name_prefix, title_text, df_trs, scaffold_for_plot, gaps_df, use_chrm=False, enhance=2500000, gap_cutoff=1000):

        
    if use_chrm:
        _df_trs = df_trs[df_trs['chrm'].isin(scaffold_for_plot['chrm'])]
    else:
        _df_trs = df_trs[df_trs['chrm'].isin(scaffold_for_plot['scaffold'])]

        
    fig = _draw_chromosomes(scaffold_for_plot, title_text, use_chrm=use_chrm)
    fig.add_trace(go.Bar(
            base=gaps_df['start'],
            x=gaps_df['length'],
            y=gaps_df['scaffold'],
            orientation="h",
            name="gaps",
            marker_color='rgba(0, 0, 0)',
        ))
    output_file_name = output_file_name_prefix + ".gaps.png"
    fig.write_image(output_file_name)
    
    fig = _draw_chromosomes(scaffold_for_plot, title_text, use_chrm=use_chrm)
    _gaps_df = gaps_df[gaps_df["length"] > gap_cutoff]
    _gaps_df["length"] = [max(x["length"], 2000000) for i,x in _gaps_df.iterrows()]
    fig.add_trace(go.Bar(
            base=_gaps_df['start'],
            x=_gaps_df['length'],
            y=_gaps_df['scaffold'],
            orientation="h",
            name="gaps",
            marker_color='rgba(0, 0, 0)',
        ))
    output_file_name = output_file_name_prefix + ".gaps.1kb.enhanced.png"
    fig.write_image(output_file_name)
    
    
        
    fig = _draw_chromosomes(scaffold_for_plot, title_text, use_chrm=use_chrm)
    _df_trs['chrm'] = [x['scaffold'] for i, x in _df_trs.iterrows()]
    names = set([x["family_name"] for i,x in _df_trs.iterrows()])
    for name in names:
        items = _df_trs[_df_trs["family_name"]==name]
        fig.add_trace(go.Bar(
            base=items['start'],
            x=items['length'],
            y=items['chrm'],
            orientation="h",
            name=name
        ))
    output_file_name = output_file_name_prefix + ".raw.png"
    fig.write_image(output_file_name)
    
    


    _title_text = title_text + " (enhanced)"
    fig = _draw_chromosomes(scaffold_for_plot, title_text, use_chrm=use_chrm)
    names = set([x["family_name"] for i,x in _df_trs.iterrows()])
    for name in names:
        items = _df_trs[_df_trs["family_name"]==name]
        items["length"] = [max(x["length"], enhance) for i,x in items.iterrows()]
        fig.add_trace(go.Bar(
            base=items['start'],
            x=items['length'],
            y=items['chrm'],
            orientation="h",
            name=name
        ))
    output_file_name = output_file_name_prefix + ".enchanced.png"
    fig.write_image(output_file_name)

    _title_text = title_text + " (enhanced, no singletons)"
    fig = _draw_chromosomes(scaffold_for_plot, title_text, use_chrm=use_chrm)
    names = set([x["family_name"] for i,x in _df_trs.iterrows()])
    for name in names:
        if name == "SING":
            continue
        items = _df_trs[_df_trs["family_name"]==name]
        items["length"] = [max(x["length"], enhance) for i,x in items.iterrows()]
        fig.add_trace(go.Bar(
            base=items['start'],
            x=items['length'],
            y=items['chrm'],
            orientation="h",
            name=name
        ))
    output_file_name = output_file_name_prefix + ".nosing.enchanced.png"
    fig.write_image(output_file_name)
    
    
    _title_text = title_text + " (no singletons)"
    fig = _draw_chromosomes(scaffold_for_plot, title_text, use_chrm=use_chrm)
    names = set([x["family_name"] for i,x in _df_trs.iterrows()])
    for name in names:
        if name == "SING":
            continue
        items = _df_trs[_df_trs["family_name"]==name]
        fig.add_trace(go.Bar(
            base=items['start'],
            x=items['length'],
            y=items['chrm'],
            orientation="h",
            name=name
        ))
    output_file_name = output_file_name_prefix + ".nosing.raw.png"
    fig.write_image(output_file_name)


def draw_all(trf_file, fasta_file, chm2name, output_folder, taxon, lenght_cutoff=10000000, level=1, enhance=1000000):

    scaffold_df = scaffold_length_sort_length(fasta_file, lenght_cutoff=lenght_cutoff)
    df_trs = read_trf_file(trf_file)

    df_trs, tr2vector, distances, all_distances  = name_clusters(df_trs, level=level)

    output_file_name = os.path.join(output_folder, f"{taxon}.trs_flow.png")
    title_text = f"Tandem repeats flow in {taxon}"
    name2monomers, name2lid, name2ids = draw_sankey(output_file_name, title_text, df_trs, tr2vector, distances, all_distances, skip_singletons=True)

    output_file_name_prefix = os.path.join(output_folder, f"{taxon}.spheres")
    title_text = f"Tandem repeats distribution in {taxon}"
    draw_spheres(output_file_name_prefix, title_text, df_trs)

    gaps_data = get_gaps_annotation(fasta_file, lenght_cutoff=lenght_cutoff)
    gaps_lengths = Counter([x[-1] for x in gaps_data])

    print("Gaps distribution:", gaps_lengths)

    gaps_df = pd.DataFrame(gaps_data, columns =['scaffold', 'start', 'end', 'length'])

    df_trs['chrm'] = [x['scaffold'] for i, x in df_trs.iterrows()]

    output_file_name_prefix = os.path.join(output_folder, f"{taxon}.karyo")
    title_text = f"Tandem repeats in {taxon}"
    draw_karyotypes(output_file_name_prefix, title_text, df_trs, scaffold_df, gaps_df, enhance=enhance, use_chrm=False, gap_cutoff=50)
