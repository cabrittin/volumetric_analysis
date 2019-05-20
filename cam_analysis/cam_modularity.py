"""
cam_tissue_mix.py

calculate the average diversity of cell clusters

"""

import sys
sys.path.append('./volumetric_analysis')
sys.path.append('.')
import argparse
import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt
import matplotlib as mpl
import networkx as nx
import igraph
import seaborn as sns
from collections import defaultdict

from sklearn.metrics.cluster import normalized_mutual_info_score
from networks.random import edge_switch

mpl.rcParams['xtick.labelsize'] = 16 
mpl.rcParams['xtick.labelsize'] = 14 

cam_class = ['cam_cad','cam_igsf','cam_lrr','cam_nrx','cam_all']
METRIC = 'pearsonr'
NODE_SCREEN = ['NSM','MC']
DIN = './cam_analysis/mat/consensus_graphs/consensus_graph_%s_p%d_deg%d.graphml'

DEG = 4 
P = 66
Q = 100 - P

AIP = DIN%('ipsi_adj',P,DEG)
AIQ = DIN%('ipsi_adj',Q,DEG)
ACP = DIN%('cont_adj',P,DEG)
ACQ = DIN%('cont_adj',Q,DEG)

AFILES = [AIP,AIQ,ACP,ACQ]

def get_communities(G,cam):
    camclass = defaultdict(lambda : [])
    for n in G.nodes():
        _camclass = G.node[n][cam]
        camclass[_camclass].append(n)

    return camclass.values()

def remove_non_classifed_nodes(G,cam):
    nodes = [n for n in G.nodes() if G.node[n][cam] == -1]
    G.remove_nodes_from(nodes)
    return G

def get_membership(G,cam):
    membership = [v[cam] for v in G.vs]
    k = max(membership) + 1
    for i in range(len(membership)):
        if membership[i] == -1: membership[i] = k + 1
    return membership

if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
 
    params = parser.parse_args()
    

    for cam in cam_class:
        mod = []
        for a in AFILES:
            A = igraph.Graph.Read_GraphML(a)
            mem = get_membership(A,cam)
            k = int(max(mem))
            #_mod = A.modularity(mem)
            #mod.append(_mod)
            vc = A.community_leading_eigenvector(clusters=k,weights='weight')
            #vc = A.community_infomap(edge_weights='weight')
            #vc = A.community_multilevel(weights='weight')
            print(vc.membership[2:8])
            mi = normalized_mutual_info_score(mem,vc.membership)
            mod.append(mi)
            #_mod = A.modularity(vc)
            #mod.append(_mod)
            edge_switch(A)
            vc = A.community_leading_eigenvector(clusters=7,weights='weight')
            #_mod = A.modularity(mem)
            mi = normalized_mutual_info_score(mem,vc.membership)
            #mod.append(mi)
            #mod.append(_mod)
        print(cam,mod)

