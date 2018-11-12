"""
tpair_syn_adj_ratio.py

Distribution of the synapse-to-adjacency ratios for the adult and
L4. Most synaptic contacts are less than 0.5 (red dashed line).

Author: Christopher Brittin
Created: 07 February 2018

"""

import sys
sys.path.append(r'./volumetric_analysis')
import matplotlib.pyplot as plt
import matplotlib as mpl

#Brittin modules
from connectome.load import from_db
import connectome.synspecificity as synspec
from networks.stats import *
from figures.stats import *
import aux

mpl.rcParams['xtick.labelsize'] = 32
mpl.rcParams['ytick.labelsize'] = 32

PRESYN_COL = '#FA5252'
POSTSYN_COL = '#F9E530'
GAPJUNC_COL = '#47A705'

db = 'N2U'
_remove = ['VC01','VD01','VB01','VB02']
P_THRESH = 0.05

lr_dict = './mat/lr_dict.txt'
homologs = './mat/homologs.txt'
left_nodes = './mat/left_nodes.txt'
right_nodes = './mat/right_nodes.txt'



def get_ratio(C,S):
    data = {'pre':[],'post':[],'gap':[]}
    
    for (u,v) in C.C.edges():
        adj = float(C.A[u][v]['count'])
        w = C.C[u][v]['weight']
        if u in S and S[u][1] < P_THRESH: data['pre'].append(w/adj)
        if v in S and S[v][2] < P_THRESH: data['post'].append(w/adj)

    for (u,v) in C.E.edges():
        adj = float(C.A[u][v]['count'])
        w = C.E[u][v]['weight']
        if (u in S) and (v in S) and (S[u][0] < P_THRESH or S[v][0] < P_THRESH):
            data['gap'].append(w/adj)

    return data


def run(fout=None):
    N2U = 'N2U'
    JSH = 'JSH'
    _remove = ['VC01','VD01','VB01','VB02']
    lrd = aux.read.into_lr_dict(lr_dict)
    left = aux.read.into_list(left_nodes)
    left.remove('CEHDL')
    left.remove('CEHVL')
    left.remove('HSNL')
    left.remove('PVNL')
    left.remove('PLNL')
    
    N2U = from_db(N2U,adjacency=True,chemical=True,
                  electrical=True,remove=_remove,dataType='networkx')
    N2U.reduce_to_adjacency()
    JSH = from_db(JSH,adjacency=True,chemical=True,
                  electrical=True,remove=_remove,dataType='networkx')
    JSH.reduce_to_adjacency()
    
    n2uspec = synspec.get_bilateral_specificity(N2U,lrd,left)
    jshspec = synspec.get_bilateral_specificity(JSH,lrd,left)   
    
    ndata = get_ratio(N2U,n2uspec)
    jdata = get_ratio(JSH,jshspec)

    data = [ndata['gap'],jdata['gap'],
            ndata['pre'],jdata['pre'],
            ndata['post'],jdata['post']]

    fig,ax = plt.subplots(1,1,figsize=(12,10))
    tpair_syn_adj_ratio(ax,data,fout=fout)
    plt.show()

if __name__ == '__main__':
    run()


