"""
cam_predict.py

Predict CAM expression patterns that result in synapses

@author Christopher Brittin
@date 06 April 2019
"""
import sys
sys.path.append('./volumetric_analysis')
import argparse
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE

import db
from connectome.load import from_db
from networks.stats import get_adj_deg,get_deg 
from cam.expression import Matrix
from connectome.skeleton import Skeleton
import aux

TCELL = 'AIAR'
SYNTYPE = 'chemical'

def degree_coefficients(C,vertices=None):
    #Create degree array.
    #Used to scale CAM connectivity graph
    #Rows ordered by vertice names  
    #Column order: [adj_deg,pre_deg,post_deg,gap_deg]

    if not vertices: vertices = sorted(C.neurons)
    adjdeg = get_adj_deg(C,vertices = vertices)
    predeg = get_deg(C.C,vertices = vertices,mode='OUT')
    pstdeg = get_deg(C.C,vertices = vertices,mode='IN')
    gapdeg = get_deg(C.E,vertices = vertices,mode='ALL')

    n = len(vertices)
    D = np.zeros((n,4))
    for i in range(n):
        v = vertices[i]
        D[i,0] = adjdeg[v]
        if v in predeg: D[i,1] = predeg[v]
        if v in pstdeg: D[i,2] = pstdeg[v]
        if v in gapdeg: D[i,3] = gapdeg[v]
    
    return D



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('db',
                        action="store",
                        help="Database name")

    parser.add_argument('matrix',
                        action = 'store',
                        help = 'Path to matrix file')

    parser.add_argument('fsize',
                        action = 'store',
                        help = 'Path to size file')

    params = parser.parse_args()

    END = None
    if params.db == 'N2U': END = 325 

    dsize = aux.read.into_list2(params.fsize)
    dtable = dict([(d[1],d) for d in dsize])

    _remove = ['VC01','VD01','VB01','VB02']
    C = from_db(params.db,adjacency=True,chemical=True,
                electrical=True,remove=_remove)
    C.remove_self_loops()
    C.reduce_to_adjacency()

    vertices = sorted(C.neurons)
    vdict = dict([(vertices[i],i) for i in range(len(vertices))])
    D = degree_coefficients(C,vertices=vertices)
    N = len(vertices)

    e = Matrix(params.matrix,cells=vertices)
    e.binarize()
    e.E = e.E.T
    E = np.copy(e.E)
    
    idx = vertices.index(TCELL)
    alpha = float(-(D[idx,0] - D[idx,1]))
    beta = float(D[idx,1])

    con = db.connect.default(params.db)
    cur = con.cursor()
    contins = db.mine.get_presynapse_contins(cur,TCELL,end=END)
    G = np.zeros((E.shape[0],1))
    Gp = np.zeros((E.shape[0],E.shape[0]))
    Gn = np.zeros((E.shape[0],E.shape[0]))
    Mp,Mn = 0.,0.
    objs = {}
    for _c in contins:
        syn = db.mine.get_synapse_data_by_contin(cur,_c)
        post,neigh = [],[]
        objtmp = []
        for s in syn:
            _post = [db.mine.get_object_contin_name(cur,p) for p in s[2]]
            _neigh = db.mine.get_object_adjacency(cur,s[1])
            post += _post
            neigh += _neigh
            objtmp.append(int(s[1]))

        post = set(post)
        mp = len(post)
        Mp += mp
        neigh = set(neigh)
        nosyn = neigh - post
        mn = len(nosyn)
        Mn += mn
        s = np.zeros((N,1))
        for v in post: 
            if v not in vdict: continue
            s[vdict[v]] = 1.  
        for v in nosyn: 
            if v not in vdict: continue
            s[vdict[v]] = -1.
   
        g = np.dot(E,s)
        gp = np.where(g > 0)
        gn = np.where(g < 0)
        tmp = np.zeros((E.shape[0],1))
        tmp[gp] = mp
        Gp += np.outer(tmp,tmp)
        tmp = np.zeros((E.shape[0],1))
        tmp[gn] = mn
        Gn += np.outer(tmp,tmp)
        for o in objtmp:
            objs[o] = gp[0]

    Gp /= Mp
    Gn /= Mn
   
    model = TSNE(n_components=2, random_state=0)
    res = model.fit_transform(Gp)
  
    ot = list(objs.keys())
    contin = str(db.mine.get_object_contin(cur,ot[0]))
    skel = Skeleton(contin)
    edges = db.mine.get_contin_edges(cur,contin,end=END)
    skel.add_edges_from(edges)
    for o in skel.nodes():
        xyz = np.array(list(map(int,db.mine.get_object_xyz(cur,o))))
        skel.node[o]['loc'] = xyz

    skel.set_distances()
    skel.length = float(dtable[contin][2])
    skel.endpts = [int(dtable[contin][3]),int(dtable[contin][4])]
    gl = np.zeros((E.shape[0],2))
    for o in ot:
        try:
            l = skel.get_node_distance(o)
            for j in  objs[o]:
                gl[j,0] += l
                gl[j,1] += 1
        except:
            print('Skipping: %d ' %o)
    
    clr = gl[:,0] / gl[:,1]

    plt.figure()
    scat = plt.scatter(res[:,0],res[:,1],c=clr,
                        cmap=plt.get_cmap('jet'),
                        vmin=0,vmax=1)
    cbar = plt.colorbar(scat)
    plt.show()

    #Ap = nx.from_numpy_matrix(np.asmatrix(Gp))
    
    
    

