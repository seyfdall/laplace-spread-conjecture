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


def generate_diameter_3_graphs(n, return_comps=False):
    """ Cycles through graphs on n vertices finding the ones with diameter 3 whose complements have diameter 3 and saves both.
        This generation is deterministic so we can trust the ordering for reruns.
        Parameters
        ----------
        n: int 
            Number of vertices to build the graphs from

        Returns
        -------
        Graphs: list(graphs)
            List of Sagemath graphs whose diameter is 3 and whose complement is also diameter 3
        Comps: list(graphs) (optional)
            The complement of the graphs above
    """
    if return_comps:
        Graphs, Comps = [], []
        for G in graphs(n):
            if G.diameter() == 3:
                H = G.complement()
                if H.diameter() == 3:
                    Graphs.append(G)
                    Comps.append(H)
        return Graphs, Comps
    
    # Split this way to save on RAM
    Graphs = []
    for G in graphs(n):
        if G.diameter() == 3:
            H = G.complement()
            if H.diameter() == 3:
                Graphs.append(G)
    return Graphs

def make_dc_graph(a,b,c):
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
    G = make_dc_graph(a,b,c)
    G_eig_vec = sorted(G.kirchhoff_matrix().eigenvectors_right(),reverse=True)[0][1][0]
    entries = [0] * len(c)
    for i in range(len(c)):
        entries[i] = G_eig_vec[c[i]]
    return sum(entries), entries


# TODO: Create generic relabelling function to relabel all graphs according to a desired order
# i.e. pass in a graph on n vertices and a map of each vertex mapping to what it should be
def relabel_graph(G, v2dc_group_map=None):
    """ Relabel graph G according to the desired mapping.  The purpose of this is to help with the data analysis
        as in, we'll be able to evaluate averages, maxes, mins and sign changes across eigenvector entries since
        each will be aligned with a different eccentricity instead of being random.

        Parameters
        ----------
        G: Sagemath graph on n vertices
            Graph to be relabelled.
        v2dc_group_map: list(n)
            Map of eccentricities of each node
            The map should have two vertices mapping to eccentricity 3 
                The graph may have more vertices of eccentricity 3, but these will just be placed in the 2 slot
            The rest of the vertices should be mapped to 1 if on the dominating edge and 2 otherwise
    """

    # First determine a pair of eccentricity 3 vertices
    x, y = -1, -1
    for v1 in G.vertices():
        for v2 in G.vertices():
            if v1 != v2 and len(G.shortest_path(v1, v2)) == 3+1: # Counts nodes instead of edges in path
                x = v1
                y = v2
                break

    #Determine the dominating edges
    dom_edges = dom_sets(G)

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
    dc_graph = make_dc_graph(len(a), len(b), len(c))

    relabel_map = {}

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
    dc_graph.relabel(relabel_map)

    pass # This is unfinished for now


"""One of our goals is, given any graph with diameter 3 whose complement also has diameter 3 to determine which complement family graph it would fall under after transformations."""

def compare_dc_family(G, plots=False):
    """Function to determine relative comp graph Family"""

    # First determine a pair of eccentricity vertices
    x, y = -1, -1
    for v1 in G.vertices():
        for v2 in G.vertices():
            if v1 != v2 and len(G.shortest_path(v1, v2)) == 3+1: # Counts nodes instead of edges in path
                x = v1
                y = v2
                break

    # Determine the dominating edges
    dom_edges = dom_sets(G)

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
    dc_graph = make_dc_graph(len(a), len(b), len(c))


    # Relabel the comp graph to match passed in graph
    v2dc_group_map = {len(a)+len(b)+len(c): 0, len(a)+len(b)+len(c)+1: 4}
    relabel_map = {len(a)+len(b)+len(c): x, len(a)+len(b)+len(c)+1: y}
    for i, vertex in enumerate(a):
        relabel_map[i] = vertex
        v2dc_group_map[i] = 1
    for i, vertex in enumerate(b):
        j = i + len(a)
        relabel_map[j] = vertex
        v2dc_group_map[j] = 3
    for i, vertex in enumerate(c):
        j = i + len(a) + len(b)
        relabel_map[j] = vertex
        v2dc_group_map[j] = 2
    dc_graph.relabel(relabel_map)

    if plots:
        G_plot = G.plot(figsize=5)
        Comp_plot = dc_graph.plot(figsize=5)
        G_plot.save('plots/G_plot.png')
        Comp_plot.save('plots/dc_plot.png')

    # TODO: What happens if multiplicity increases?
    G_spec = np.array(G.spectrum(laplacian=True), dtype=float)
    dc_spec = np.array(dc_graph.spectrum(laplacian=True), dtype=float)
    G_eign = np.abs(float(G_spec[0]))
    dc_eign = np.abs(float(dc_spec[0]))

    if G_eign > dc_eign:
        print("Counterexample found: Bad News")
        G_plot = G.plot(figsize=5)
        dc_plot = dc_graph.plot(figsize=5)
        G_plot.save('plots/G_plot.png')
        dc_plot.save('plots/dc_plot.png')

    # TODO: This is only getting the first eigenvector associated with spec radius,
    # What should we do if the spec radius multiplicity is greater than 1?
    G_eigvec = np.array(sorted(G.eigenvectors(laplacian=True))[-1][1][0], dtype=float)
    dc_eigvec = np.array(sorted(dc_graph.eigenvectors(laplacian=True))[-1][1][0], dtype=float)

    # Normalize the vectors for comparison
    G_eigvec = G_eigvec / np.linalg.norm(G_eigvec)
    dc_eigvec = dc_eigvec / np.linalg.norm(dc_eigvec)

    return G_eign, dc_eign, G_eigvec, dc_eigvec, G_spec, dc_spec, v2dc_group_map


def compare_sym_dc_family(G, plots=False):
    """Function to determine relative symmetric comp graph Family"""
    # TODO: Absorb this function into the one above if possible for code simplification

    # First determine a pair of eccentricity vertices
    x, y = -1, -1
    for v1 in G.vertices():
        for v2 in G.vertices():
            if v1 != v2 and len(G.shortest_path(v1, v2)) == 3+1: # Counts nodes instead of edges in path
                x = v1
                y = v2
                break

    # Determine the dominating edges
    dom_edges = dom_sets(G)

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

    # End if not an even amount of vertices
    # TODO: How should we handle an odd total of dominating vertices?  Skipping for now
    dom_vert_count = len(a) + len(b)
    if dom_vert_count % 2 != 0:
        return
    
    vertices_to_remove = set()
    # Determine c group (shared between dominating edges)
    for node in remaining_vertices:
        if len(G.shortest_path(node,x)) > 1 and len(G.shortest_path(node,y)) > 1:
            c.add(node)
            vertices_to_remove.add(node)
    remaining_vertices = remaining_vertices.difference(vertices_to_remove)
    vertices_to_remove.clear()

    # Create comparison complement family and truss graph
    dc_graph = make_dc_graph(len(a), len(b), len(c))
    dc_sym_graph = make_dc_graph(dom_vert_count // 2, dom_vert_count // 2, len(c))

    # Relabel the comp graph to match passed in graph
    v2dc_group_map = {len(a)+len(b)+len(c): 0, len(a)+len(b)+len(c)+1: 4}
    relabel_map = {len(a)+len(b)+len(c): x, len(a)+len(b)+len(c)+1: y}
    for i, vertex in enumerate(a):
        relabel_map[i] = vertex
        v2dc_group_map[i] = 1
    for i, vertex in enumerate(b):
        j = i + len(a)
        relabel_map[j] = vertex
        v2dc_group_map[j] = 3
    for i, vertex in enumerate(c):
        j = i + len(a) + len(b)
        relabel_map[j] = vertex
        v2dc_group_map[j] = 2
    dc_graph.relabel(relabel_map)

    # TODO: May need to update this to account for the switching around to fit symmetry
    v2dc_sym_group_map = {len(a)+len(b)+len(c): 0, len(a)+len(b)+len(c)+1: 4}
    relabel_sym_map = {len(a)+len(b)+len(c): x, len(a)+len(b)+len(c)+1: y}
    for i, vertex in enumerate(a):
        relabel_sym_map[i] = vertex
        v2dc_sym_group_map[i] = 1
    for i, vertex in enumerate(b):
        j = i + len(a)
        relabel_sym_map[j] = vertex
        v2dc_sym_group_map[j] = 3
    for i, vertex in enumerate(c):
        j = i + len(a) + len(b)
        relabel_sym_map[j] = vertex
        v2dc_sym_group_map[j] = 2
    dc_sym_graph.relabel(relabel_sym_map)

    if plots:
        G_plot = G.plot(figsize=5)
        Comp_plot = dc_graph.plot(figsize=5)
        Comp_sym_plot = dc_sym_graph.plot(figsize=5)
        G_plot.save('plots/G_plot.png')
        Comp_plot.save('plots/dc_plot.png')
        Comp_sym_plot.save('plots/dc_sym_plot.png')

    # TODO: What happens if multiplicity increases?
    G_spec = np.array(G.spectrum(laplacian=True), dtype=float)
    dc_spec = np.array(dc_graph.spectrum(laplacian=True), dtype=float)
    dc_sym_spec = np.array(dc_graph.spectrum(laplacian=True), dtype=float)
    G_eign = np.abs(float(G_spec[0]))
    dc_eign = np.abs(float(dc_spec[0]))
    dc_sym_eign = np.abs(float(dc_sym_spec[0]))

    if G_eign > dc_eign:
        print("Asymmetric Counterexample found: Bad News")
        G_plot = G.plot(figsize=5)
        dc_plot = dc_graph.plot(figsize=5)
        G_plot.save('plots/G_plot.png')
        dc_plot.save('plots/dc_plot.png')

    if dc_eign > dc_sym_eign:
        print("Symmetric Counterexample found: Bad News")
        dc_plot = dc_graph.plot(figsize=5)
        dc_sym_plot = dc_sym_graph.plot(figsize=5)
        dc_plot.save('plots/dc_plot.png')
        dc_sym_plot.save('plots/dc_sym_plot.png')

    # TODO: This is only getting the first eigenvector associated with spec radius,
    # What should we do if the spec radius multiplicity is greater than 1?
    G_eigvec = np.array(sorted(G.eigenvectors(laplacian=True))[-1][1][0], dtype=float)
    dc_eigvec = np.array(sorted(dc_graph.eigenvectors(laplacian=True))[-1][1][0], dtype=float)
    dc_sym_eigvec = np.array(sorted(dc_sym_graph.eigenvectors(laplacian=True))[-1][1][0], dtype=float)

    # Normalize the vectors for comparison
    G_eigvec = G_eigvec / np.linalg.norm(G_eigvec)
    dc_eigvec = dc_eigvec / np.linalg.norm(dc_eigvec)
    dc_sym_eigvec = dc_sym_eigvec / np.linalg.norm(dc_sym_eigvec)

    return G_eign, dc_eign, dc_sym_eign, G_eigvec, dc_eigvec, dc_sym_eigvec, G_spec, dc_spec, dc_sym_spec, v2dc_group_map, v2dc_sym_group_map


"""
Graph search methods
"""

def dc_search(ns, tf, dir_path):
    """Search and compare between graphs and their complement family versions"""
    t0 = time.time()

    for n in ns:
        file_path = dir_path + f'/results/results_on_{n}.h5'
        with h5py.File(file_path, 'w') as file:
            Graphs = generate_diameter_3_graphs(n)
            for i, G in enumerate(Graphs):
                t1 = time.time()
                if t1 - t0 > tf:
                    print("Break in Combo")
                    break
                G_eign, dc_eign, G_eigvec, dc_eigvec, G_spec, dc_spec, v2dc_group_map = compare_dc_family(G)
                v2dc_group_arr = np.array([v2dc_group_map[i] for i in range(n)])

                group = file.create_group(f"Graph_{i}_on_{n}_vertices")
                group.attrs['G_spec_radius'] = G_eign
                group.attrs['DC_spec_radius'] = dc_eign
                group.create_dataset('G_eigvec', data=G_eigvec)
                group.create_dataset('DC_eigvec', data=dc_eigvec)
                group.create_dataset('G_spec', data=G_spec)
                group.create_dataset('DC_spec', data=dc_spec)
                group.create_dataset('v2dc_group_map', data=v2dc_group_arr)
            
        if t1 - t0 > tf:
            print("Break in Combo")
            break
            

def dc_sym_search(ns, tf, dir_path):
    """Search and compare between graphs, their complement family versions, and their symmetric versions"""
    t0 = time.time()

    for n in ns:
        file_path = dir_path + f'/sym_results/sym_results_on_{n}.h5'
        print(file_path)
        with h5py.File(file_path, 'w') as file:
            Graphs = generate_diameter_3_graphs(n)
            for i, G in enumerate(Graphs):
                t1 = time.time()
                if t1 - t0 > tf:
                    print("Break in Combo")
                    break

                results = compare_sym_dc_family(G)
                if not results:
                    continue

                G_eign, dc_eign, dc_sym_eign, G_eigvec, dc_eigvec, dc_sym_eigvec, G_spec, dc_spec, dc_sym_spec, v2dc_group_map, v2dc_sym_group_map = results
                v2dc_group_arr = np.array([v2dc_group_map[i] for i in range(n)])
                v2dc_sym_group_arr = np.array([v2dc_sym_group_map[i] for i in range(n)])

                group = file.create_group(f"Graph_{i}_on_{n}_vertices")
                group.attrs['G_spec_radius'] = G_eign
                group.attrs['DC_spec_radius'] = dc_eign
                group.attrs['DC_sym_spec_radius'] = dc_sym_eign
                group.create_dataset('G_eigvec', data=G_eigvec)
                group.create_dataset('DC_eigvec', data=dc_eigvec)
                group.create_dataset('DC_sym_eigvec', data=dc_sym_eigvec)
                group.create_dataset('G_spec', data=G_spec)
                group.create_dataset('DC_spec', data=dc_spec)
                group.create_dataset('DC_sym_spec', data=dc_sym_spec)
                group.create_dataset('v2dc_group_map', data=v2dc_group_arr)
                group.create_dataset('v2dc_sym_group_map', data=v2dc_sym_group_arr)
            
        if t1 - t0 > tf:
            print("Break in Combo")
            break


"""
Main Method
"""

if __name__ == "__main__":
    # DC := Dandelion Complement
    ns = range(6,9)
    tf = 85000
    dir_path = os.path.dirname(os.path.realpath(__file__))
    dc_search(ns, tf, dir_path)
    # dc_sym_search(ns, tf, dir_path)
    


