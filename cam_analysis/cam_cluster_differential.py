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

mpl.rcParams['xtick.labelsize'] = 16 
mpl.rcParams['xtick.labelsize'] = 14 

idx_gene = {'cad':range(85,98),'lrr':range(54,85),
            'igsf':range(54),'nrx':range(98,106),
            'all':range(106)}

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

def assign_cam_class(i,cam_bounds):
    for j in range(len(cam_bounds)):
        if i < cam_bounds[j]: return j
    return j+1


cam_class = {'cad':CAD_CLASS,'igsf':IGSF_CLASS,'lrr':LRR_CLASS,
            'nrx':NRX_CLASS,'all':ALL_CLASS}

cam_color = {'cad':CAD_COLOR,'igsf':IGSF_COLOR,'lrr':LRR_COLOR,
            'nrx':NRX_COLOR,'all':ALL_COLOR}


REMOVE = ['VB01', 'VD01']
FOUT = 'mat/cam_class/consensus_cam_class_all_tissue_%s_%s.csv'

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
    neurons = sorted(ML.load_reduced_nodes())
    camclass = ML.load_cam_class(params.metric,params.camtype)
    
    e = Matrix(ML.cam,params.matrix)
    e.load_genes()
    e.load_cells(nodes)
    e.assign_expression()
    e.clean_expression()
    e.binarize()
   
    gdx = idx_gene[params.camtype]
    cdx = [[] for i in range(max(camclass.values()) + 1)]
    for (k,v) in camclass.items(): 
        cdx[v].append(e.cells[k])
    
    M = e.M[:,gdx]
    data = []
    maxdata = []
    mindata = []
    std = []
    for c in cdx: 
        data.append(np.mean(M[c,:],axis=0))
        maxdata.append(np.max(M[c,:],axis=0))
        mindata.append(np.min(M[c,:],axis=0))
        std.append(np.std(M[c,:],axis=0))

    genes = [e.gene_idx[i] for i in gdx]
    k = len(cdx)
    width = 1. / (k+1)
    idx = np.arange(M.shape[1])
    fig,ax = plt.subplots(1,1,figsize=(17,10))
    ax.axvspan(3-width,5-width,facecolor='#969696')
    ax.axvspan(12-width,13-width,facecolor='#969696')
    ax.axvspan(6-width,7-width,facecolor='#e4e4e4')
    for i in range(k):
        ax.bar(idx + i*width,data[i],width,yerr=[mindata[i],maxdata[i]],
                color=cam_color[params.camtype][i],label='Cluster %d'%i )
    ax.set_xticks(idx + k/2*width)
    ax.set_xticklabels(genes)
    ax.set_title('Average cadherin expression by cluster',fontsize=18)
    ax.set_ylabel('Average counts',fontsize=18)
    ax.legend(fontsize=18)
    plt.tight_layout()
    plt.show()
    
