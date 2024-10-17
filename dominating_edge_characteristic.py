"""
    Sagemath mamba installation guide found on https://doc.sagemath.org/html/en/installation/conda.html#sec-installation-conda
    Python 3.12 was too new - fell back to Python 3.10 so the install command was: 'mamba create -n sage sage python=3.10'
"""
import numpy as np
import math 
import networkx as nx
import itertools
import time
from sage.all import *
import h5py
import os


#Graphs and complements are put into list for easy ordering/sorting
def find_graphs(n):
    Graphs = []
    Comps = []
    for G in graphs(n):
        if G.diameter()==3:
            H=G.complement()
            if H.diameter()==3:
                Graphs.append(G)
                Comps.append(H)
    return Graphs, Comps

"""Class for Observing the Fiedler entries on a graph"""
class GraphObserver:

    def __init__(self, G):
        self.G = G.copy()
        self.fiedler_vec = sorted(self.G.kirchhoff_matrix().eigenvectors_right(),reverse=True)[-2][1][0]

    def show_graph(self):
        labels = []
        for i, label in enumerate(self.G.vertices()):
            labels.append(str(label) + ': ' + '{0:.2f}'.format(float(self.fiedler_vec[i])))
        self.G.plot(figsize=7, vertex_size=1800.0, vertex_labels=dict(zip(self.G, labels))).save('dominating_edge_characteristic_results/G_plot.png')

    def dom_edges(self):
        return [i for i in self.G.dominating_sets() if self.G.has_edge(i)]
    
    def comp_dom_edges(self):
        # Check to see if each pair has different sign
        pairs = self.dom_edges()
        found_counterexample = False
        for pair in pairs:
            a = pair[0]
            b = pair[1]
            vec_a = self.fiedler_vec[a]
            vec_b = self.fiedler_vec[b]
            if (vec_a > 0 and vec_b > 0) or (vec_a < 0 and vec_b < 0):
                print(pair)
                found_counterexample = True
                break
        return found_counterexample

if __name__ == "__main__":
    n = 8
    for i in range(n,9):
        Graphs, Comps = find_graphs(i)
        print("n:", i)
        for j, graph in enumerate(Graphs):
            print(j)
            observer = GraphObserver(graph)
            found_counterexample = observer.comp_dom_edges()
            if found_counterexample:
                observer.show_graph()