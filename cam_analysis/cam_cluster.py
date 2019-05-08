"""
cam_cluster.py

Looks for combinations of frequently occuring CAM gene combinations

"""

import sys
sys.path.append('./volumetric_analysis')
sys.path.append('.')
import argparse
import numpy as np
from itertools import combinations
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt
import networkx as nx

from cam.expression import Matrix
from mat_loader import MatLoader
from connectome.load import from_db
import aux

idx_gene = {'cad':range(85,98),'lrr':range(54,85),
            'igsf':range(54),'nrx':range(98,106),
            'all':range(106)}

def cad_class(i):
    if i < 15: return 0
    elif i < 54: return 1
    else: return 2

def igsf_class(i):
    if i < 11: return 0
    elif i < 21: return 1
    elif i < 29: return 2
    elif i < 38: return 3
    elif i < 48: return 4
    elif i < 57: return 5
    else: return 6

def lrr_class(i):
    if i < 10: return 0
    elif i < 17: return 1
    elif i < 35: return 2
    elif i < 40: return 3
    elif i < 51: return 4
    else: return 5

def nrx_class(i):
    if i < 6: return 0
    elif i < 22: return 1
    elif i < 37: return 2
    elif i < 46: return 3
    else: return 4

def all_class(i):
    if i < 8: return 0
    elif i < 16: return 1
    elif i < 21: return 2
    elif i < 41: return 3
    elif i < 50: return 4
    elif i < 56: return 5
    elif i < 62: return 6
    else: return 7

cam_class = {'cad':cad_class,'igsf':igsf_class,'lrr':lrr_class,
            'nrx':nrx_class,'all':all_class}

def overlap(u,v):
    w = u + v
    m = float(len(np.where(w > 0)[0]))
    n = len(np.where(w > 1)[0])
    if m > 0: return 1 - n / m
    return 1

REMOVE = ['VB01', 'VD01']
FOUT = 'mat/cam_class/consensus_cam_class_%s_%s.csv'

if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('matrix',
                        action = 'store',
                        help = 'Path to matrix file')
    
    parser.add_argument('metric',
                        action = 'store',
                        default='pearsonr',
                        choices=['pearsonr','jaccard','hamming'],
                        help = 'Specify metric')

    parser.add_argument('camtype',
                        action = 'store',
                        default='all',
                        choices=['all','cad','lrr','igsf','nrx'],
                        help = 'Specify CAM choice')

    params = parser.parse_args()
    
    ML = MatLoader()
    ML.load_lrmap()
    #nodes = sorted(ML.load_reduced_nodes())
    nodes = sorted(ML.load_all_tissue())

    e = Matrix(ML.cam,params.matrix)
    e.load_genes()
    e.load_cells(nodes)
    e.assign_expression()
    e.clean_expression()
    e.binarize()
    #M = e.E[:,idx_gene]
    #for i in range(M.shape[0]): np.random.shuffle(M[i,:])

    D = e.distance_matrix(gdx=idx_gene[params.camtype],metric=params.metric)
        
    Y = sch.linkage(D, method='ward')
    Z = sch.dendrogram(Y, orientation='right')

    fig,axmatrix = plt.subplots(1,1,figsize=(10,10))

    index = Z['leaves']
    D = D[index,:]
    D = D[:,index]
    im = axmatrix.matshow(D, aspect='auto', origin='lower')
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])
    plt.colorbar(im)

    data = []
    idx = 0
    for i in index: 
        print(idx,e.cells_idx[i],cam_class[params.camtype](idx),i)
        data.append([e.cells_idx[i],cam_class[params.camtype](idx)])
        idx += 1
    fout = FOUT%(params.metric,params.camtype)
    #aux.write.from_list(fout,data) 
    plt.show()
    
