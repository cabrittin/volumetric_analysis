"""
cam_contact.py

Look for CAM expression that is over represented in larger contact.

@author Christopher Brittin 
@date 28 April 2019
"""

import sys
sys.path.append('.')
sys.path.append('./volumetric_analysis')
import argparse
import numpy as np
import networkx as nx
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import umap

from mat_loader import MatLoader
from cam.expression import Matrix
import cam.cam_predict as predict
from connectome.load import from_db


def contact_profile(e,A,cell,alpha=0.33):
    neigh = sorted(A[cell].items(), key=lambda edge: edge[1]['weight'])
    k = len(neigh)
    _low =  [n[0] for n in neigh[:int(k*alpha)]]
    _high =  [n[0] for n in neigh[-int(k*alpha):]]
    low = [e.cells[n[0]] for n in neigh[:int(k*alpha)]]
    high = [e.cells[n[0]] for n in neigh[-int(k*alpha):]]
    elow = np.sum(e.E[low,:],axis=0)
    ehigh = np.sum(e.E[high,:],axis=0)

    num = ehigh - elow
    den = 0.5*(ehigh + elow)
    idx = np.where(den > 0)
    ratio = -9 * np.ones(e.E.shape[1])
    ratio[idx] = num[idx] / den[idx]
    jdx = np.where(ratio > 1.)[0].tolist()
    genes = [e.gene_idx[j] for j in jdx]
    print(_low)
    print(_high)
    print(genes) 

def generate_umap_plot(ax,M,clr,nodes=None,colorbar=False):
    print(M.shape)
    model = umap.UMAP(n_neighbors = 15, n_components=2,min_dist=0.10,metric='jaccard',
            random_state=12)
    res = model.fit_transform(M)

    if nodes:
        for i in range(len(nodes)): print(nodes[i],res[i,0],res[i,1])
    scat = ax.scatter(res[:,0],res[:,1],c=clr,
                cmap=plt.get_cmap('jet'),
                vmin=0,vmax=1)
    if colorbar:
        cbar = plt.colorbar(scat)

DB = ['JSH','N2U']
REMOVE = ['VB01', 'VD01']
cmap = {'SMN':0.0,'HMNa':0.1, 'HMNp':0.2,'I1':0.4, 'I2':0.5, 
        'Sa':0.7, 'Sp1':0.9, 'Sp2':0.9, 'XMN':1.0}

#cam_test = range(85,98)
#cam_test = range(54,85)
#cam_test = range(54)
#cam_test = range(98,106)
cam_test = range(106)

if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('matrix',
                        action = 'store',
                        help = 'Path to matrix file')
    
    #parser.add_argument('db',
    #                    action = 'store',
    #                    help = 'Database')

    #parser.add_argument('cell',action='store',help='Cell name')

    params = parser.parse_args()
    #cell = params.cell

    M = MatLoader()
    nclass = M.load_nerve_ring_classes()

    #N = from_db('N2U',adjacency=True,dataType='networkx',remove=REMOVE)
    #J = from_db('JSH',adjacency=True,dataType='networkx',remove=REMOVE)

    #nodes = sorted(set(N.A.nodes()) & set(J.A.nodes()))
    nodes = M.load_reduced_nodes()
    nrclass = [nclass[n] for n in nodes]
    print(sorted(set(nrclass)))

    e = Matrix(M.cam,params.matrix)
    e.load_genes()
    e.load_cells(nodes)
    e.assign_expression()
    e.binarize()
    #np.random.shuffle(e.M)
    
    
    #Gn = nx.to_numpy_matrix(N.A,nodelist=nodes,weight='weight')
    #Gj = nx.to_numpy_matrix(J.A,nodelist=nodes,weight='weight')

    #G = 0.5*(Gn + Gj)
    #G[G<500] = 0

    #model = TSNE(n_components=2, random_state=0)
    clr = [cmap[nclass[n]] for n in nodes]

    M = e.E[:,cam_test]
    #np.random.shuffle(M)
    fig,ax = plt.subplots(1,2,figsize=(15,10))
    generate_umap_plot(ax[1],M,clr,nodes=nodes,colorbar=True)
    np.random.shuffle(M)
    generate_umap_plot(ax[0],M,clr)
    ax[0].set_title('Shuffled genes')
    plt.show()
    
    #contact_profile(e,C.A,cell)
    

