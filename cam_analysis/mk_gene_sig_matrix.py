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
import networkx as nx

from mat_loader import MatLoader
from cam.expression import Matrix
import cam.cam_predict as predict
import aux

cam = 'mat/cam_isoforms.txt'
FOUT = 'cam_analysis/gene_sigs/sig_%s_deg%d.graphml'

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
    M.load_lrmap()
    
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
        for cell in tqdm(M.left,desc='Cells'):
            if not C.C.has_node(cell): continue
            syn,neigh = predict.get_synapse_data(S[cell],e,
                                cpartners=set(C.C.neighbors(cell)))
            if not syn: continue
            gene_sig = set(predict.gene_differential(e.E,syn,neigh))        
            if not gene_sig: 
                print('Skippin cell: ' + cell)
                continue
            
            for n in C.A.neighbors(cell):
                idx = e.cells[n]
                ssig = set(np.where(e.E[idx,:] > 0)[0].tolist())
                score = predict.score_overlap(gene_sig,ssig)
                D.add_edge(cell,n,weight=score) 
                if params.mode == 'post':
                    D.add_edge(M.lrmap[cell],n,weight=score)

        fout = FOUT%(params.mode,deg)
        nx.write_graphml(D,fout)

    
