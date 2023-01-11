#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 11.01.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

from trseeker.tools.sequence_tools import get_revcomp
import math

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

def name_clusters(df_trs):
    token2id, token2revtoken = get_pentatokens()
    tr2vector = fill_vectors(df_trs, token2id, token2revtoken, k=5)
    distances = compute_distances(tr2vector)

    all_distances = list(distances.values())
    all_distances = set(map(all_distances, int))
    all_distances.sort(reverse=True)

    G = Graph(list(tr2vector.keys()))
    for (id1,id2) in distances:
        G.addEdge(id1, id2, distances[(id1,id2)])

    for i in range(1, 0, -1):
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
            name = f"{class_name}_{d[0][1]}" 
            for id1, period in d:
                df_trs.at[id1, 'family_name'] = name
        for id1, period in singl:
            df_trs.at[id1, 'family_name'] = "SING"

    return df_trs
        
        