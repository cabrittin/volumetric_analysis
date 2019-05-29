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
import matplotlib as mpl
import networkx as nx
import seaborn as sns
from collections import defaultdict

from cam.expression import Matrix
from mat_loader import MatLoader
from connectome.load import from_db
import aux

mpl.rcParams['xtick.labelsize'] = 4

idx_gene = {'cad':range(85,98),'lrr':range(54,85),
            'igsf':range(54),'nrx':range(98,106),
            'all':range(106),
            'cad_lron':list(range(69,84)) + list(range(85,98)),
            'non_cad_lron':list(range(69)) + [85] + list(range(98,106))}
CAD_CLASS = [7,21,34,46,74,82,103]
CAD_COLOR = ['#EC360F','#fac127','#f0fc3c','#70fc3c','#3cfcf0','#3c70fc','#ad3cfc','#fc3cdf']

IGSF_CLASS = [10,23,45,60,65,83,89]
IGSF_COLOR = ['#EC360F','#fac127','#f0fc3c','#70fc3c','#3cfcf0','#3c70fc','#ad3cfc','#fc3cdf']

LRR_CLASS = [39,64,84,93,103,109,112]
LRR_COLOR = ['#EC360F','#fac127','#f0fc3c','#70fc3c','#3cfcf0','#3c70fc','#ad3cfc','#fc3cdf']

NRX_CLASS = [25,41,66,77]
NRX_COLOR = ['#EC360F','#f8d94d','#4df864','#55c0fa','#8255fa','#000000']

ALL_CLASS = [11,20,26,39,63,84,98]
ALL_COLOR = ['#EC360F','#fac127','#f0fc3c','#70fc3c','#3cfcf0','#3c70fc','#ad3cfc','#fc3cdf']

CAD_LRON_CLASS = [28,46,52,70,91,97]
CAD_LRON_COLOR = ['#EC360F','#fac127','#f0fc3c','#70fc3c','#3cfcf0','#3c70fc','#ad3cfc']  

NON_CAD_LRON_CLASS = [19,30,33,43,58,83,98]   
NON_CAD_LRON_COLOR = ['#EC360F','#fac127','#f0fc3c','#70fc3c','#3cfcf0','#3c70fc','#ad3cfc','#fc3cdf']   

def assign_cam_class(i,cam_bounds):
    for j in range(len(cam_bounds)):
        if i < cam_bounds[j]: return j
    return j+1


cam_class = {'cad':CAD_CLASS,'igsf':IGSF_CLASS,'lrr':LRR_CLASS,
    'nrx':NRX_CLASS,'all':ALL_CLASS,'cad_lron':CAD_LRON_CLASS,
    'non_cad_lron':NON_CAD_LRON_CLASS}

cam_color = {'cad':CAD_COLOR,'igsf':IGSF_COLOR,'lrr':LRR_COLOR,
    'nrx':NRX_COLOR,'all':ALL_COLOR,'cad_lron':CAD_LRON_COLOR,
    'non_cad_lron':NON_CAD_LRON_COLOR}

def overlap(u,v):
    w = u + v
    m = float(len(np.where(w > 0)[0]))
    n = len(np.where(w > 1)[0])
    if m > 0: return 1 - n / m
    return 1

REMOVE = ['VB01', 'VD01']
FOUT = 'mat/cam_class/consensus_cam_class_all_tissue_%s_%s.csv'
NODE_SCREEN = ['NSM','MC']

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
                        choices=['all','cad','lrr','igsf','nrx','cad_lron','non_cad_lron'],
                        help = 'Specify CAM choice')

    parser.add_argument('-o','--output',
                        action='store',
                        dest='fout',
                        required=False,
                        default=None,
                        help='Output file')

    params = parser.parse_args()
    
    ML = MatLoader()
    ML.load_lrmap()
    #nodes = sorted(ML.load_reduced_nodes())
    nodes = sorted(ML.load_all_tissue())
    neurons = sorted(ML.load_reduced_nodes()) + NODE_SCREEN
   
    e = Matrix(ML.cam,params.matrix)
    e.load_genes()
    e.load_cells(nodes)
    e.assign_expression()
    e.clean_expression()
    e.binarize()
    #M = e.E[:,idx_gene]
    #for i in range(M.shape[0]): np.random.shuffle(M[i,:])
    
    #for (i,g) in  e.gene_idx.items(): print(i,g)
    
    D = e.distance_matrix(gdx=idx_gene[params.camtype],metric=params.metric)

    Y = sch.linkage(D, method='ward')
    Z = sch.dendrogram(Y, orientation='right')

    #fig,axmatrix = plt.subplots(1,1,figsize=(10,10))

    index = Z['leaves']
    k = len(e.cells)
    ordered_cells = [None]*k
    tclass = ['k']*k
    cclass = [0]*k
    ccolor = ['k']*k
    idx = 0
    for i in index:
        ordered_cells[i] = e.cells_idx[i]
        if e.cells_idx[i] in neurons: tclass[i] = '#E6E6E6'
        cclass[i] = assign_cam_class(idx,cam_class[params.camtype])
        ccolor[i] = cam_color[params.camtype][cclass[i]]
        idx += 1
    D2 = D[index,:]
    D2 = D2[:,index]
    im = sns.clustermap(D,row_linkage=Y,col_linkage=Y,xticklabels=ordered_cells,yticklabels=[],
            row_colors=tclass,col_colors=ccolor,figsize=(12,12))
    #im.ax_row_dendrogram.set_visible(False)
    #im = axmatrix.matshow(D2, aspect='auto', origin='lower')
    #axmatrix.set_xticks([])
    #axmatrix.set_yticks([])
    #plt.colorbar(im)

    if params.fout: 
        print(params.fout)
        im.savefig(params.fout,dpi=400)

    data = []
    idx = 0
    for i in index: 
        print(idx,e.cells_idx[i],cclass[i],i)
        data.append([e.cells_idx[i],cclass[i]])
        idx += 1
    fout = FOUT%(params.metric,params.camtype)
    aux.write.from_list(fout,data) 
    plt.show()
    
