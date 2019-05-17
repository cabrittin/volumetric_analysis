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
        Aip = igraph.Graph.Read_GraphML(AIP)
        Aiq = igraph.Graph.Read_GraphML(AIQ)
        Acp = igraph.Graph.Read_GraphML(ACP)
        Acq = igraph.Graph.Read_GraphML(ACQ)
        
        mem_ip = get_membership(Aip,cam)
        mem_iq = get_membership(Aiq,cam)
        mem_cp = get_membership(Acp,cam)
        mem_cq = get_membership(Acq,cam)

        mip = Aip.modularity(mem_ip)
        miq = Aiq.modularity(mem_iq)
        mcp = Acp.modularity(mem_cp)
        mcq = Acq.modularity(mem_cq)

        edge_switch(Aip)
        edge_switch(Aiq)
        edge_switch(Acp)
        edge_switch(Acq)

        mipr = Aip.modularity(mem_ip)
        miqr = Aiq.modularity(mem_iq)
        mcpr = Acp.modularity(mem_cp)
        mcqr = Acq.modularity(mem_cq)
        
        print(Aip.ecount(),Aiq.ecount(),Acp.ecount(),Acq.ecount())
        print(cam,mip,mipr,miq,miqr,mcp,mcpr,mcq,mcqr)




