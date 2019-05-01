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
from mpl_toolkits.mplot3d import Axes3D
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
 
def pca_basis(profile,syn,nonsyn):
    from sklearn.decomposition import PCA
    n = len(profile)
    X = np.array([v for (k,v) in profile.items()])
    pca = PCA(n_components=5,svd_solver='full')
    pca.fit(X)
    print(pca.explained_variance_ratio_)
    X = pca.transform(X)
    Y = pca.transform(nonsyn)
    Z = pca.transform(syn)
    fig = plt.figure(1,figsize=(15,10))
    print(Z.shape,X.shape)
    ax = Axes3D(fig)#rect=[0,0,0.95,1],elev=48,azim=134)
    ax.scatter(X[:,0],X[:,1],X[:,2],facecolor='b')
    ax.scatter(Y[:,0],Y[:,1],Y[:,2],facecolor='r')
    ax.scatter(Z[:,0],Z[:,1],Z[:,2],facecolor='g')
    #plt.scatter(X[:,0],X[:,1],facecolor='b')
    #plt.scatter(Y[:,0],Y[:,1],facecolor='r') 
    plt.show()

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
    #e.shuffle_rows()
    (n,m) = e.E.shape

    cell = params.cell
    syn,neigh,cneigh = predict.get_synapse_data(S[cell],e,
                                cpartners=set(C.C.neighbors(cell)),remove_partners=True)
    gene_sig = predict.gene_differential(e.E,syn,neigh)
    gene_mask = np.zeros(m)
    gene_mask[gene_sig] = 1
    gene_sig = predict.gene_mean_profile(e.E,syn,neigh)
    
    profile = np.multiply(gene_mask,gene_sig)
    
    x,cost_rec = predict.run_optimation(e.E,syn,neigh,alpha=[0.4,0.4,0.2])
    x *= 0.5
    #for c in gene_sig:
    #    profile[e.cells_idx[c]] = np.multiply(gene_mask,gene_sig[c])

    jdx = e.cells[cell]
    cexp = e.E[jdx,:]
    fig,ax = plt.subplots(3,1,figsize=(10,15))
    plot_profile(ax[1],cexp,profile)
    ax[1].set_xlabel('CAM gene index')
    ax[1].set_ylabel('Expression profile')
    plot_profile(ax[2],cexp,x)
    ax[2].set_xlabel('CAM gene index')
    ax[2].set_ylabel('Optimized profile')
    ax[0].plot(cost_rec,'b-')
    ax[0].set_ylabel('Cost')
    ax[0].set_xlabel('Iteration')
    
    plt.show()

    """
    syn = set(C.C.neighbors(cell))
    nonsyn = set(C.A.neighbors(cell)) - syn
    ndata = np.zeros((len(nonsyn),m))
    i = 0
    for n in nonsyn:
        idx = e.cells[n]
        ndata[i,:] = np.multiply(e.E[idx,:],gene_mask)
        #ndata[i,:] = e.E[idx,:]
        i += 1

    print(syn)
    sdata = np.zeros((len(syn),m))
    i = 0 
    for n in syn:
        idx = e.cells[n]
        sdata[i,:] = np.multiply(e.E[idx,:],gene_mask)
        #sdata[i,:] = e.E[idx,:]
        i =+ 1
    
    pca_basis(profile,sdata,ndata)
    
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
    """

