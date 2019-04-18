"""
over_represented_genes.py

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

from mat_loader import MatLoader
from cam.expression import Matrix



cam = 'mat/cam_isoforms.txt'
FOUT = 'cam_analysis/results/cam_over_represent_deg%d.csv'

def get_gene_count(E,syn,neigh):
    (n,m) = E.shape
    k = len(syn)
    syn_count = np.zeros(m)
    neigh_count = np.zeros(m)
    for i in range(k):
        sdx = syn[i]
        ndx = neigh[i]
        ssum = np.sum(E[sdx,:],axis=0)
        ssum[ssum > 0] = 1
        syn_count += ssum
        nsum = np.sum(E[ndx,:],axis=0)
        nsum[nsum > 0] = 1
        neigh_count += nsum
        diff = ssum - nsum
        tmp = E[sdx,:] - nsum 
        #print('gene diff',np.where(diff > 0))
        #print(sdx,ndx,np.where(diff > 1))

    diff = syn_count - neigh_count
    print('Sum gene diff',np.where(diff > 0))
    print(k)
    
    #reldff = 0.5*(syn_count - neigh_count) / (syn_count + neigh_count)
    #print('Rel diff', np.where(reldff))

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

    
    e = Matrix(cam,params.matrix)
    e.load_genes()
    e.load_cells(sorted(C.A.nodes()))
    e.assign_expression()
    e.binarize()
    print(len(e.gene_idx),len(e.cells_idx),e.E.shape)
    syn,neigh = [],[]
    for cell in [params.cell]:#S:
        cneigh = set(C.C.neighbors(cell))
        for cont in S[cell]:
            partners = set(S[cell][cont]['partners'])
            neighbors = set(S[cell][cont]['neighbors'])
            nonsyn = neighbors - partners - cneigh
            syn.append([e.cells[n] for n in partners])
            neigh.append([e.cells[n] for n in nonsyn])


    get_gene_count(e.E,syn,neigh)        

    

