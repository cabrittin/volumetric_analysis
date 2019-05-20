"""
cam_tissue_mix.py

calculate the average diversity of cell clusters

"""

import sys
sys.path.append('./volumetric_analysis')
sys.path.append('.')
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import networkx as nx
import igraph
import seaborn as sns
from sklearn.metrics.cluster import normalized_mutual_info_score
from itertools import combinations

from networks.random import edge_switch

mpl.rcParams['xtick.labelsize'] = 16 
mpl.rcParams['xtick.labelsize'] = 14 

cam_class = ['cam_cad','cam_igsf','cam_lrr','cam_nrx','cam_all']
METRIC = 'pearsonr'
NODE_SCREEN = ['NSM','MC']
DIN = './cam_analysis/mat/consensus_graphs/consensus_graph_%s_p%d_deg%d.graphml'
EIN = './cam_analysis/mat/cam_distance_graphs/cam_distance_%s_%s.graphml'  

DEG = 4 
P = 66
Q = 100 - P

AIP = DIN%('ipsi_adj',P,DEG)
AIQ = DIN%('ipsi_adj',Q,DEG)
ACP = DIN%('cont_adj',P,DEG)
ACQ = DIN%('cont_adj',Q,DEG)

AFILES = [AIP,AIQ,ACP,ACQ]

def binarize_data(G,E,thresh=0.5):
    """
    Returns vectors of edges (1) and non-edges (0)
    for graphs G and E. 

    Paramters:
    ----------
    G : adjacency graph (networkx)
    E : cam expression graph (networkx)
    thresh: threshold for cam correaltion. Values above thresh set to 1
          Values below thresh set to 0.
    """

    nodes = set(G.nodes()) & set(E.nodes())
    k = len(nodes)
    m = k*(k-1)
    g,e = np.zeros(m),np.zeros(m)
    idx = 0
    for (u,v) in combinations(nodes,2):
        if G.has_edge(u,v): g[idx] = 1
        if E[u][v]['weight'] > thresh: e[idx] = 1
        idx += 1
    return g,e


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
 
    params = parser.parse_args()
    

    for cam in cam_class:
        E = nx.read_graphml(EIN%(METRIC,cam))
        mi = []
        for a in AFILES:
            G = nx.read_graphml(a)
            g,e = binarize_data(G,E,thresh=.4)
            _mi = normalized_mutual_info_score(g,e)
            mi.append(_mi)
        
        print(cam,mi) 

    
