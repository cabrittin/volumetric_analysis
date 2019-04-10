"""
test_split_graph.py

Splits graph based on left right nodes

@author Christoper Brittin
@date 09 April 2019
"""

import sys
sys.path.append('./volumetric_analysis')
sys.path.append('.')
import networkx as nx

from mat_loader import MatLoader
from connectome.load import from_db


if __name__=="__main__":
    _db = 'N2U'
    
    M = MatLoader()
    M.load_left()
    M.load_right()
    M.load_lrmap()
    C = from_db(_db,adjacency=True,chemical=True,
                electrical=True,dataType='networkx')
    C.split_left_right(M.left,M.right) 
    
    print(C.A.number_of_edges())
    print(C.Al.number_of_edges())
    print(C.Ar.number_of_edges())
    print(C.C.number_of_edges())
    print(C.Cl.number_of_edges())
    print(C.Cr.number_of_edges())
    print(C.E.number_of_edges())
    print(C.El.number_of_edges())
    print(C.Er.number_of_edges())
    
    C.map_right_graphs(M.lrmap)
    print(C.Cr.edges())

    
    

