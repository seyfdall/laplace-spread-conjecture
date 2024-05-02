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


def make_comp_graph(a,b,c):
    """ Creates a complement family graph based on three sets of vertices in each of the three groups (with a bottom left, b bottom right, c top) in the complement family

        Returns
        ------
        G: Graph
            Complement family graph with groups a,b,c
    """
    if type(a) == int:
        tmp_a = [i for i in range(a)]
        tmp_b = [i + a for i in range(b)]
        tmp_c = [i + a + b for i in range(c)]
        a = tmp_a
        b = tmp_b
        c = tmp_c

    n = len(a)+len(b)+len(c)
    G = graphs.CompleteGraph(n)

    for v in a:
        G.add_edge(n,v)

    for v in b:
        G.add_edge(n+1,v)

    return G

def equit_part_eigs(a,b,c):
    """ Creates equitable parititon matrix
    """
    if type(a) != int:
        a = len(a)
        b = len(b)
        c = len(c)
    equitable = np.array([[a,-a,0,0,0],[-1,b+c+1,-c,-b,0],[0,-a,a+b,-b,0],[0,-a,-c,a+c+1,-1],[0,0,0,-b,b]])
    eigvals, eigvecs = np.linalg.eig(equitable)
    return eigvals,eigvecs

def dom_sets(G):
    return [i for i in G.dominating_sets() if G.has_edge(i)]

def test_eigenvector(a,b,c):
    """Tests the largest eigenvector of a graph in the complement family to see if the entries in c's partition are 0

       Returns
       -------
       sum(entries) : float
           The sum of c's entries in the largest eigenvector of the graph
       entries : list(float)
           Each of c's entries in the eigenvector
    """
    G = make_comp_graph(a,b,c)
    G_eig_vec = sorted(G.kirchhoff_matrix().eigenvectors_right(),reverse=True)[0][1][0]
    entries = [0] * len(c)
    for i in range(len(c)):
        entries[i] = G_eig_vec[c[i]]
    return sum(entries), entries



"""One of our goals is, given any graph with diameter 3 whose complement also has diameter 3 to determine which complement family graph it would fall under after transformations."""

def compare_comp_family(G,plots = False):
    """Function to determine relative comp graph Family"""

    # First determine a pair of eccentricity vertices
    x, y = -1, -1
    for v1 in G.vertices():
        for v2 in G.vertices():
            if v1 != v2 and len(G.shortest_path(v1, v2)) == 3+1: # Counts nodes instead of edges in path
                x = v1
                y = v2
                break

    #Determine the dominating edges
    dom_edges = dom_sets(G)

    # TODO: Begin forming groups a,b,c,d,e from remaining vertices
    # Get remaining vertices
    remaining_vertices = set(G.vertices()) - set([x,y])
    for edge in dom_edges:
        remaining_vertices = remaining_vertices - set(edge)

    # Cycle through remaining vertices and assign their groups
    a = set()
    b = set()
    c = set()

    # Determine a and b groups (dominating nodes)
    for edge in dom_edges:
        v1, v2 = edge
        if len(G.shortest_path(v1,x)) == 2:
            a.add(v1)
            b.add(v2)
        else:
            a.add(v2)
            b.add(v1)

    vertices_to_remove = set()
    # Determine c group (shared between dominating edges)
    for node in remaining_vertices:
        if len(G.shortest_path(node,x)) > 1 and len(G.shortest_path(node,y)) > 1:
            c.add(node)
            vertices_to_remove.add(node)
    remaining_vertices = remaining_vertices.difference(vertices_to_remove)
    vertices_to_remove.clear()

    # Create comparison complement family and truss graph
    comp_graph = make_comp_graph(len(a), len(b), len(c))

    # Relabel the comp graph to match passed in graph
    relabel_map = {len(a)+len(b)+len(c): x, len(a)+len(b)+len(c)+1: y}
    for i, vertex in enumerate(a):
        relabel_map[i] = vertex
    for i, vertex in enumerate(b):
        j = i + len(a)
        relabel_map[j] = vertex
    for i, vertex in enumerate(c):
        j = i + len(a) + len(b)
        relabel_map[j] = vertex
    comp_graph.relabel(relabel_map)

    if plots:
        G.plot(figsize=5).show()
        comp_graph.plot().show()
    G_eign = G.spectrum(laplacian=True)[0]
    comp_eign = comp_graph.spectrum(laplacian=True)[0]

    print(comp_eign)
    if G_eign > comp_eign:
        print("G_eign > comp_eign", G_eign, truss_eign)
        G.plot().show()

    G_eign, G_eigvec, mult = sorted(G.eigenvectors(laplacian=True))[-1]
    comp_eigvec = comp_graph.eigenvectors(laplacian=True)[0]


"""
Main Method
"""

if __name__ == "__main__":
    print("Start")
    G = make_comp_graph(int(3),int(2),int(1))
    print(dom_sets(G))
    G_plot = G.plot(figsize=5)
    G_plot.save('comp_test.png')