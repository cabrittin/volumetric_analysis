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


from cam.expression import Matrix
from mat_loader import MatLoader


idx_gene = range(85,98)
#idx_gene = range(54,85)
#idx_gene = range(54)
#idx_gene = range(98,106)
#idx_gene = range(106)

def overlap(u,v):
    w = u + v
    m = float(len(np.where(w > 0)[0]))
    n = len(np.where(w > 1)[0])
    if m > 0: return 1 - n / m
    return 1


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('matrix',
                        action = 'store',
                        help = 'Path to matrix file')
    
    parser.add_argument('db',
                        action = 'store',
                        help = 'Database')
    
    params = parser.parse_args()
    
    ML = MatLoader()
    ML.load_lrmap()
    nodes = sorted(ML.load_reduced_nodes())
    sa = ML.load_sa(params.db)
        
    #C = from_db(params.db,adjacency=True,dataType='networkx',remove=REMOVE)
    
    e = Matrix(ML.cam,params.matrix)
    e.load_genes()
    e.load_cells(nodes)
    e.assign_expression()
    e.binarize()
    #np.random.shuffle(e.M)
   
    i1 = e.cells['IL2VL']
    i2 = e.cells['OLQVL']

    M = e.E[:,idx_gene]
    print(M[i1,:] + M[i2,:])

    k = len(nodes)
    D = np.ones((k,k))
    np.fill_diagonal(D,0)
    ndict = {nodes[i]:i for i in range(k)}
    comb = combinations(nodes,2)
    for (u,v) in comb:
        i = e.cells[u]
        j = e.cells[v]
        r = overlap(M[i,:],M[j,:])
        if np.isnan(r): continue
        D[ndict[u],ndict[v]] = r
        D[ndict[v],ndict[u]] = r

   
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

    idx = 0
    for i in index: 
        print(idx,nodes[i],i)
        idx += 1
    plt.show()
    
