import sys
sys.path.append('./volumetric_analysis')
sys.path.append('.')
import networkx as nx
from tqdm import tqdm
import numpy as np

from connectome.load import from_db
from mat_loader import MatLoader
from cam.expression import Matrix


def consensus_graph(G,H,deg,nodes):
    for n in nodes:
        neigh = set([])
        for h in H:
            if not h.has_node(n): continue
            neigh = neigh.union(set(list(h.neighbors(n))))
        for m in neigh:
            _deg = 0
            _w = 0.
            for h in H:
                if h.has_edge(n,m):
                    _deg += 1
                    _w += h[n][m]['weight']
            if _deg >= deg:
                G.add_edge(n,m,weight=_w/_deg)
    return G

def filter_graph(A,C,E,pct=50):
    H = nx.Graph()
    nodes = sorted(A.nodes())
    weights = [A[u][v]['weight'] for (u,v) in A.edges()]
    thresh = np.percentile(weights,pct)

    for (u,v) in A.edges():
        c1 = not C.has_edge(u,v)
        c2 = not C.has_edge(v,u)
        c3 = not E.has_edge(u,v)
        c4 = thresh > A[u][v]['weight']
        if c1 and c2 and c3 and c4: continue
        H.add_edge(u,v,weight=A[u][v]['weight'])
    return H

def filter_graph_simple(A,pct=50,thresh_high=True):
    H = nx.Graph()
    nodes = sorted(A.nodes())
    weights = [A[u][v]['weight'] for (u,v) in A.edges()]
    thresh = np.percentile(weights,pct)

    for (u,v) in A.edges():
        if thresh_high: 
            c = thresh > A[u][v]['weight']
        else:
            c = thresh <= A[u][v]['weight']
        if c: continue
        H.add_edge(u,v,weight=A[u][v]['weight'])
    return H


 
def apply_node_class(G,nclass,attr,default=-1):
    for n in G.nodes():
        G.node[n][attr] = default 
        if n in nclass: G.node[n][attr] = nclass[n]
 

def split_lateral(G,left,right):
    """
    Split ipsilater from contralateral connections
    """
    H = nx.Graph()
    if G.is_directed(): H = nx.DiGraph()

    for (u,v) in G.edges():
        c1 = (u in left) and (v in right)
        c2 = (u in right) and (v in left)
        if c1 or c2: H.add_edge(u,v,weight=G[u][v]['weight'])
    
    G.remove_nodes_from(right)
    return G,H

DOUT = './cam_analysis/mat/consensus_graphs/consensus_graph_%s_p%d_deg%d.graphml'
REMOVE = ['VB01', 'VD01','VC01']
CAM_TYPES = {'cam_all':'all','cam_cad':'cad','cam_igsf':'igsf',
            'cam_lrr':'lrr','cam_nrx':'nrx'}
METRIC = 'pearsonr'
deg = 3
P = 66
Q = 100 - P

if __name__=="__main__":
    M = MatLoader()
    M.load_left()
    M.load_right()
    M.load_lrmap()

    nclass = M.load_nerve_ring_classes()
    nodes = M.load_reduced_nodes()
    single = set(nodes) - set(M.left) - set(['ASER'])

    N = from_db('N2U',adjacency=True,chemical=True,electrical=True,
            remove=REMOVE,dataType='networkx')
    Nsa = N.split_graph(N.A,single)
    Nsc = N.split_graph(N.C,single)
    Nse = N.split_graph(N.E,single)
    N.split_left_right(M.left,M.right)
    N.map_right_graphs(M.lrmap)

    J = from_db('JSH',adjacency=True,chemical=True,
                electrical=True,dataType='networkx')
    Jsa = J.split_graph(J.A,single)
    Jsc = J.split_graph(J.C,single)
    Jse = J.split_graph(J.E,single)
    J.split_left_right(M.left,M.right)
    J.map_right_graphs(M.lrmap)
 
    A = nx.Graph()
    C = nx.DiGraph()
    E = nx.Graph()
    A = consensus_graph(A,[N.Al,N.Ar,J.Al,J.Ar],deg,M.left)
    A = consensus_graph(A,[Nsa,Jsa],min(deg,2),single)
    C = consensus_graph(C,[N.Cl,N.Cr,J.Cl,J.Cr],deg,M.left)
    C = consensus_graph(C,[Nsc,Jsc],min(deg,2),single)
    E = consensus_graph(E,[N.El,N.Er,J.El,J.Er],deg,M.left)
    E = consensus_graph(E,[Nse,Jse],min(deg,2),single)
        
    Ai,Ac = split_lateral(A,M.left,M.right)
    Ci,Cc = split_lateral(C,M.left,M.right)
    Ei,Ec = split_lateral(E,M.left,M.right)
        
    Aip = filter_graph_simple(Ai,pct=P)
    Acp = filter_graph_simple(Ac,pct=P)
    Aiq = filter_graph_simple(Ai,pct=Q,thresh_high=False)
    Acq = filter_graph_simple(Ac,pct=Q,thresh_high=False)
       
    #for G in [A,C,E]: G.remove_nodes_from(M.right)
    for G in [Aip,Aiq,Acp,Acq]: apply_node_class(G,nclass,'cell_class',default='NA')
        
    for cam in CAM_TYPES: 
        cdict = M.load_cam_class(METRIC,CAM_TYPES[cam])
        for (l,r) in zip(M.left,M.right):
            if l in cdict: cdict[r] = cdict[l]
        for G in [Aip,Aiq,Acp,Acq]: apply_node_class(G,cdict,cam,default=-1)
        
    nx.write_graphml(Aip,DOUT%('ipsi_adj',P,deg))
    nx.write_graphml(Acp,DOUT%('cont_adj',P,deg))
    nx.write_graphml(Aiq,DOUT%('ipsi_adj',Q,deg))
    nx.write_graphml(Acq,DOUT%('cont_adj',Q,deg))


