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
from matplotlib.patches import Rectangle
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
            'all':range(106),
            'cad_lron':list(range(69,84)) + list(range(85,98)),
            'non_cad_lron':list(range(69)) + [85] + list(range(98,106))}
CAD_CLASS = [7,21,34,46,74,82,103]
CAD_COLOR = ['#EC360F','#fac127','#f0fc3c','#70fc3c','#3cfcf0','#3c70fc','#ad3cfc','#fc3cdf']

IGSF_CLASS = [10,23,45,60,65,83,89]
IGSF_COLOR = ['#EC360F','#fac127','#f0fc3c','#70fc3c','#3cfcf0','#3c70fc','#ad3cfc','#fc3cdf']

LRR_CLASS = [39,64,84,93,103,109,112]
LRR_COLOR = ['#EC360F','#fac127','#f0fc3c','#70fc3c','#3cfcf0','#3c70fc','#ad3cfc','#fc3cdf']

#NRX_CLASS = [25,41,66,77]
#NRX_COLOR = ['#EC360F','#f8d94d','#4df864','#55c0fa','#8255fa','#000000']
NRX_CLASS = [15,30,42,54,65]
NRX_COLOR = ['#EC360F','#f8d94d','#4df864','#107902','#55c0fa','#8255fa']

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
    neurons = sorted(ML.load_reduced_nodes())
    camclass = ML.load_cam_class(params.metric,params.camtype)
    
    e = Matrix(ML.cam,params.matrix)
    e.load_genes()
    e.load_cells(nodes)
    e.assign_expression()
    e.clean_expression()
    e.binarize()
    
    _tmp = aux.read.into_list('mat/non_cadlron_differential.txt')
    gdx = [e.genes[g].idx for g in _tmp]
    #gdx = idx_gene[params.camtype]
    cdx = [[] for i in range(max(camclass.values()) + 1)]
    for (k,v) in camclass.items(): 
        cdx[v].append(e.cells[k])
    M = e.M[:,gdx]
    matrix = np.zeros((len(gdx),len(cdx)))
    data = []
    maxdata = []
    mindata = []
    std = []
    i = 0
    for c in cdx:
        matrix[:,i] = np.mean(M[c,:],axis=0)
        data.append(np.mean(M[c,:],axis=0))
        maxdata.append(np.max(M[c,:],axis=0))
        mindata.append(np.min(M[c,:],axis=0))
        std.append(np.std(M[c,:],axis=0))
        i += 1

    #print(matrix)
    #martrix = np.log10(matrix)
    #print(matrix)
    matrix[matrix<1] = 1
    matrix = np.log10(matrix)
    genes = [e.gene_idx[i] for i in gdx]
    clusters = range(1,len(cdx)+1)

    cmap = sns.cubehelix_palette(np.max(matrix)*100,rot=-.3, reverse=True,
                                    start=0)
    g= sns.clustermap(matrix,row_cluster=False,col_cluster=False,col_colors=cam_color[params.camtype],
                    yticklabels=genes,xticklabels=clusters,cmap=cmap,linewidth=1,#figsize=(13,15),
                    cbar_kws={'label': 'log(average adjusted counts)'})
    
    g.cax.yaxis.label.set_size(20)
    g.cax.set_position([.15, .2, .03, .45])
    #g.ax_heatmap.add_patch(Rectangle((1, 3), 2, 2, fill=False, edgecolor='#ffe400', lw=3))
    if params.fout: g.savefig(params.fout)
    plt.show()
