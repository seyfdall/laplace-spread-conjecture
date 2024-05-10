import numpy as np
import math 
import networkx as nx
import itertools
import time
from sage.all import *
import h5py


"""
Main Method
"""

if __name__ == "__main__":
    # Example of how to read in data from an h5df file in Python
    with h5py.File('results_on_6.h5', 'r') as file:
        for i in range(len(Graphs)):
            group = file[f'Graph_{i}_on_{n}_vertices']
            print([attribute for attribute in group.attrs])
            print(group.attrs['DC_spec_radius'])
            print(group['G_eigvec'][:])
            print('\n')
