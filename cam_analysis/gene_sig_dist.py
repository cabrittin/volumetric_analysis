"""
gene_sig_dist.py

Gets the gene signatures for all cell classes

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
import networkx as nx


from mat_loader import MatLoader
from cam.expression import Matrix
import cam.cam_predict as predict
import aux

cam = 'mat/cam_isoforms.txt'
FOUT = 'cam_analysis/gene_sigs/gene_sig_%s_deg%d.csv'

DEG = [1,2,3,4]
if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('matrix',
                        action = 'store',
                        help = 'Path to matrix file')
    
    parser.add_argument('-m','--mode',
                        action='store',
                        dest = 'mode',
                        choices = ['pre','post','gap'],
                        default = 'pre',
                        required = False,
                        help = 'Specify pre, post or gap junction synaptic')
   
    params = parser.parse_args()

    M = MatLoader()
    M.load_left()
    
    for deg in DEG:
        print('Degre: %d' %deg)
        C = M.load_consensus_graphs(deg)
        if params.mode == 'pre':
            S = M.load_consensus_chemical_synapse(deg)
        elif params.mode == 'post':
            S = M.load_consensus_chemical_post_synapse(deg)
        else:
            S = M.load_consensus_gap_junctions(deg)
    
        e = Matrix(cam,params.matrix)
        e.load_genes()
        e.load_cells(sorted(C.A.nodes()))
        e.assign_expression()
        e.binarize()
        
        D = nx.DiGraph()
        SSIG,NSIG,IDSCORE = [],[],[]
        GENE_SIG = {}
        for cell in tqdm(M.left,desc='Cells'):
            if not C.C.has_node(cell): continue
            syn,neigh = predict.get_synapse_data(S[cell],e,cpartners=set(C.C.neighbors(cell)))
            if not syn: continue
            gene_sig = predict.gene_differential(e.E,syn,neigh)        
            if not gene_sig: 
                GENE_SIG[cell] = ['-999']
                continue
            
            ssig,nsig,idscore = predict.get_overlap(gene_sig,e.E,syn,neigh)
            SSIG += ssig
            NSIG += nsig
            IDSCORE.append(idscore)
            GENE_SIG[cell] = [e.gene_idx[i] for i in gene_sig]

        fout = FOUT%(params.mode,deg)
        aux.write.from_dict(fout,GENE_SIG)

        print(np.mean(IDSCORE),np.std(IDSCORE),np.min(IDSCORE),np.max(IDSCORE))    

    
