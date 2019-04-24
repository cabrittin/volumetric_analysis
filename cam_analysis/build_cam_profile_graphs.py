"""
build_cam_profile_graphs.py

Builds expression profile graphs

@author Christopher Brittin
@date 24 April 2019
"""

import sys
sys.path.append('.')
sys.path.append('./volumetric_analysis')
import argparse
from tqdm import tqdm
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

from mat_loader import MatLoader
from cam.expression import Matrix
import cam.cam_predict as predict

cam = 'mat/cam_isoforms.txt'
#FOUT = 'cam_analysis/data/profile_graphs/cam_profile_deg%d.graphml'
FOUT = '/home/cabrittin/Dropbox/PhD/sr_vol/figs2/cam_profiles/profile_cell%s_deg%d.png'


def plot_profile(ax,pre_exp,post_exp,title=None):
    idx = range(len(pre_exp))
    ax.bar(idx,pre_exp,color='r')
    ax.bar(idx,post_exp,color='b')
    ax.set_ylim([0,1])
    #ax.set_xlabel('CAM gene index')
    #ax.set_ylabel('Expression profile')
    if title: ax.set_title(title)
 


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
    (n,m) = e.E.shape

    cell = params.cell
    syn,neigh,cneigh = predict.get_synapse_data(S[cell],e,cpartners=set(C.C.neighbors(cell)))
    gene_sig = predict.gene_differential(e.E,syn,neigh)
    gene_mask = np.zeros(m)
    gene_mask[gene_sig] = 1
    gene_sig = predict.gene_profile(e.E,syn,neigh)
    
    profile = {}
    mean_profile = np.zeros(m)
    for c in gene_sig:
        profile[e.cells_idx[c]] = np.multiply(gene_mask,gene_sig[c])
        tmp = np.multiply(gene_mask,gene_sig[c])
        mean_profile += tmp 

    mean_profile /= len(gene_sig)
  
    jdx = e.cells[cell]
    cexp = e.E[jdx,:]
    M = int(np.ceil(C.C.out_degree(cell) / 2))
    fig,ax = plt.subplots(M,2,figsize=(20,10),sharex=True,sharey=True)
    ax = ax.flatten()
    idx = 0
    for (k,v) in profile.items():
        title = '%s -> %s' %(cell,k)
        plot_profile(ax[idx],cexp,v,title=title)
        if idx % 2 == 0: ax[idx].set_ylabel('Expression profile') 
        idx += 1
    ax[-1].set_xlabel('CAM gene index')
    ax[-2].set_xlabel('CAM gene index')
    plt.tight_layout()
    plt.savefig(FOUT%(cell,params.deg))
    plt.show()


