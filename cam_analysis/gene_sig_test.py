"""
gene_sig_test.py

Test gene signature model

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


DEG = [1,2,3,4]
if __name__=='__main__':

    M = MatLoader()
    M.load_left()

    for deg in DEG:
        C = M.load_consensus_graphs(deg)
        S = M.load_consensus_chemical_synapse(deg)
        pre,post = M.load_gene_sig_graph(deg)
        SCORE = [[],[]]
        for cell in tqdm(M.left,desc='Cells'):
            if not C.C.has_node(cell): continue
            if cell not in S: continue
            if pre.out_degree(cell) == 0: continue
            cneighbors = set(C.C.neighbors(cell))
            k,s1,s2 = 0.,0.,0.
            for cont in S[cell]:
                partners = set(S[cell][cont]['partners'])
                neighbors = set(S[cell][cont]['neighbors'])
                nonsyn = neighbors - partners - cneighbors
                #if not nonsyn: continue
                pscore = [0] + [pre[cell][p]['weight'] for p in partners]
                nscore = [0] + [pre[cell][p]['weight'] for p in nonsyn]
                k += 1
                if np.max(pscore) > np.max(nscore): s1 += 1
                
                pscore,nscore = [0],[0]
                for p in partners:
                    if not post.has_edge(p,cell):
                        pscore.append(pre[cell][p]['weight'])
                    else:
                        pscore.append(pre[cell][p]['weight'] + post[p][cell]['weight'])

                for p in nonsyn:
                    if not post.has_edge(p,cell):
                        nscore.append(pre[cell][p]['weight'])
                    else:
                        nscore.append(pre[cell][p]['weight'] + post[p][cell]['weight'])
                if np.max(pscore) > np.max(nscore): s2 += 1

            if k > 0:
                SCORE[0].append(s1 / k)
                SCORE[1].append(s2 / k)
            else: 
                print(cell)
        
        print(np.mean(SCORE[0]),np.std(SCORE[0]),np.min(SCORE[0]),np.max(SCORE[0]))
        print(np.mean(SCORE[1]),np.std(SCORE[1]),np.min(SCORE[1]),np.max(SCORE[1]))

    
