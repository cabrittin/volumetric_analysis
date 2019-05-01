"""
sbe_cluster.py

Look for genes that are over represented relative to random synapse partners


@author Christopher Brittin
@date 17 April 2019
"""

import sys
sys.path.append('.')
sys.path.append('./volumetric_analysis')
import argparse
import numpy as np
from tqdm import tqdm
from random import shuffle
import matplotlib.pyplot as plt

from mat_loader import MatLoader
from cam.expression import Matrix
import cam.cam_predict as predict


cam = 'mat/cam_isoforms.txt'
FOUT = 'cam_analysis/results/cam_over_represent_deg%d.csv'

if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('matrix',
                        action = 'store',
                        help = 'Path to matrix file')

    parser.add_argument('deg',
                        action = 'store',
                        type= int,
                        help = 'Conserved degree')

   
    parser.add_argument('cell',action='store',help='Cell name')

    params = parser.parse_args()


    M = MatLoader()
    M.load_left()
    C = M.load_consensus_graphs(params.deg)
    S = M.load_consensus_chemical_synapse(params.deg)
    S1 = M.load_consensus_chemical_synapse(1)
    C1 = M.load_consensus_graphs(1)
    print(C.C.number_of_edges(),C1.C.number_of_edges())
    
    e = Matrix(cam,params.matrix)
    e.load_genes()
    e.load_cells(sorted(C.A.nodes()))
    e.assign_expression()
    e.binarize()
    print(len(e.gene_idx),len(e.cells_idx),e.E.shape)

    syn,neigh,cneigh = predict.get_synapse_data(S[params.cell],e,cpartners=set(C.C.neighbors(params.cell)))

    gene_sig = predict.gene_differential(e.E,syn,neigh)        
        

    
