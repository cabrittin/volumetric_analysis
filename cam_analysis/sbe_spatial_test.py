"""
sbe_spatial_test.py

Show that spatial constraints are required for SBE


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
FOUT = 'cam_analysis/results/sbe_spatial_pre_train_adult_test_l4_deg%d.csv'

DEG = [1,2,3,4]
rand_iter = 1000


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('matrix',
                        action = 'store',
                        help = 'Path to matrix file')

    params = parser.parse_args()
    
    for deg in DEG:
        print('Degree: %d' %deg)

        M = MatLoader()
        M.load_left()
        C = M.load_consensus_graphs(deg)
        S = M.load_consensus_chemical_synapse(deg)
        S1 = M.load_consensus_chemical_synapse(1)
        C1 = M.load_consensus_graphs(1)
    
        e = Matrix(cam,params.matrix)
        e.load_genes()
        e.load_cells(sorted(C.A.nodes()))
        e.assign_expression()
        e.binarize()
        
        data = {}
        IDSCORE1,IDSCORE2 = [],[]
        RSCORE1,RSCORE2 = [],[]
        for cell in tqdm(M.left,desc='Cells: '):
            if not C.C.has_node(cell): continue
            syn,neigh = predict.get_synapse_data(S[cell],e,
                cpartners=set(C.C.neighbors(cell)),screen='N2U')
            gene_sig = predict.gene_differential(e.E,syn,neigh)  
            if not gene_sig:
                data[cell] = [-1,-1]
                continue
            jsyn,jneigh = predict.get_synapse_data(S1[cell],e,screen='JSH') 
            if not jsyn:
                data[cell] = [-2,-2]
                continue
            ssig,nsig,idscore1 = predict.get_overlap(gene_sig,e.E,jsyn,jneigh)
            
            rscore1 = []
            for i in range(rand_iter):
                rscore1.append(predict.get_random_overlap(jsyn,jneigh))

            jsyn,jneigh = predict.get_synapse_data(S1[cell],e,
                    cpartners=set(C1.C.neighbors(cell)),screen='JSH')
            ssig,nsig,idscore2 = predict.get_overlap(gene_sig,e.E,jsyn,jneigh)
            
            rscore2 = []
            for i in range(rand_iter):
                rscore2.append(predict.get_random_overlap(jsyn,jneigh))

            
            data[cell] = [idscore1,idscore2,
                    np.mean(rscore1),np.std(rscore1),
                    np.mean(rscore2),np.std(rscore2)]
            IDSCORE1.append(idscore1)
            RSCORE1.append(np.mean(rscore1))
            IDSCORE2.append(idscore2)
            RSCORE2.append(np.mean(rscore2))
            

        print(np.mean(IDSCORE1),np.std(IDSCORE1),np.mean(RSCORE1),np.std(RSCORE1))
        print(np.mean(IDSCORE2),np.std(IDSCORE2),np.mean(RSCORE2),np.std(RSCORE2))
