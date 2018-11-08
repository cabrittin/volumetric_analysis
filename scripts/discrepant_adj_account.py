"""
discrepant_adj_account.py

Fraction of discrepant synaptic connections that occur at 
discrepant adjacency connections.

Fraction of discrepant adjacency connections with a 
synaptic connections.

created: Christopher Brittin
date: 01 November 2018

"""
import sys
sys.path.append(r'./volumetric_analysis')
import matplotlib.pyplot as plt
import matplotlib as mpl
import networkx as nx

#Brittin modules
from connectome.load import from_db
import connectome.synspecificity as synspec
from networks.stats import *
from figures.stats import *
import aux

mpl.rcParams['xtick.labelsize'] = 32
mpl.rcParams['ytick.labelsize'] = 32

ADULT_COL = '#FFC300'
L4_COL = '#3380FF'
AL_COL = '#14FA29'

left_nodes = './mat/left_nodes.txt'
right_nodes = './mat/right_nodes.txt'
lr_dict = './mat/lr_dict.txt'

def get_discrepant_contra_directed(G1,G2,lr,left,right):
    no_adj,total = 0,0.
    for (u,v) in G1.edges():
        if (u not in left) and (u not in right): continue
        ur,vr = lr[u],lr[v]
        if not G1.has_edge(ur,vr):
            total += 1
            if not G2.has_edge(ur,vr):
                no_adj += 1
                
    return [no_adj / total, no_adj, total]

def get_discrepant_contra_undirected(G1,G2,lr,left,right):
    no_adj,total = 0,0
    for (u,v) in G1.edges():
        if (((u not in left) and (u not in right)) or
            ((v not in left) and (v not in right))):
            continue
        ur,vr = lr[u],lr[v]
        if not G1.has_edge(ur,vr):
            total += 1
            if not G2.has_edge(ur,vr):
                no_adj += 1
    return [no_adj / total, no_adj, total]

def get_discrepant_dev(C1,C2,A2):
    no_adj,total = 0,0.
    for (u,v) in C1.edges():
        if not C2.has_edge(u,v):
            total += 1
            if not A2.has_edge(u,v):
                no_adj += 1
    return [no_adj / total, no_adj, total]

def get_discrepant_adj_with_chem(A,C,lr,left,right):
    has_chem,total = 0,0.
    for (u,v) in A.edges():
        if (((u not in left) and (u not in right)) or
            ((v not in left) and (v not in right))):
            continue
        ur,vr = lr[u],lr[v]
        if not A.has_edge(ur,vr):
            total += 1
            if C.has_edge(u,v) or C.has_edge(v,u):
                has_chem += 1
    return [has_chem / total, has_chem, total]

def get_discrepant_adj_with_gap(A,E,lr,left,right):
    has_gap,total = 0,0.
    for (u,v) in A.edges():
        if (((u not in left) and (u not in right)) or
            ((v not in left) and (v not in right))):
            continue
        ur,vr = lr[u],lr[v]
        if not A.has_edge(ur,vr):
            total += 1
            if E.has_edge(u,v):
                has_gap += 1
    return [has_gap / total, has_gap, total]

def get_discrepant_adj_with_chem_dev(A1,A2,C1):
    has_chem,total = 0,0.
    for (u,v) in A1.edges():
        if not A2.has_edge(u,v):
            total += 1
            if C1.has_edge(u,v) or C1.has_edge(v,u):
                has_chem += 1
    return [has_chem / total, has_chem, total]

def get_discrepant_adj_with_gap_dev(A1,A2,E1):
    has_gap,total = 0,0.
    for (u,v) in A1.edges():
        if not A2.has_edge(u,v):
            total += 1
            if E1.has_edge(u,v):
                has_gap += 1
    return [has_gap / total, has_gap, total]

def make_plot(ax,data,width=0.15,xlim=[-0.5,1.5],
              ylim=[0,0.3],title=None):
    ind = np.arange(2)
    ax.axvspan(-0.5,0.5,facecolor='#C3C3C3')
    ax.axvspan(0.5,1.5,facecolor='#D8D7D7')
    ax.grid(axis='y',color='0.9',linestyle='-',linewidth=1)
    rects0 = ax.bar(ind-width, data[0], width, color=ADULT_COL)
    rects1 = ax.bar(ind, data[1], width, color=L4_COL)
    rects2 = ax.bar(ind+width, data[2], width, color=AL_COL)
    ax.set_xticks(ind)
    ax.set_xticklabels(('gap j.','chemical'))
    ax.set_xlim(xlim)
    ax.set_yticks([0.1,0.2,0.3])
    ax.set_ylim(ylim)
    if title: ax.set_title(title,fontsize=28)
    ax.legend((rects0[0],rects1[0],rects2[0]),
              ('Adult L/R', 'L4 L/R', 'Adult/L4'),
              fontsize=18)
    
def run(fout=None):
    left = aux.read.into_list(left_nodes)
    right = aux.read.into_list(right_nodes)
    _lrd = aux.read.into_dict(lr_dict)
    lrd = {}
    for key,val in _lrd.items():
        lrd[key] = val
        lrd[val] = key

    _remove = ['VC01','VD01','VB01','VB02']

    n2u = from_db('N2U',adjacency=True,chemical=True,electrical=True,
                  dataType='networkx',remove=_remove)
    jsh = from_db('JSH',adjacency=True,chemical=True,electrical=True,
                  dataType='networkx',remove=_remove)

    nsyn = get_discrepant_contra_directed(n2u.C,n2u.A,lrd,left,right)
    ngap = get_discrepant_contra_undirected(n2u.E,n2u.A,lrd,left,right)
    
    jsyn = get_discrepant_contra_directed(jsh.C,jsh.A,lrd,left,right)
    jgap = get_discrepant_contra_undirected(jsh.E,jsh.A,lrd,left,right)

    bsyn = get_discrepant_dev(n2u.C,jsh.C,jsh.A)
    bgap = get_discrepant_dev(n2u.E,jsh.E,jsh.A)

    nadjc = get_discrepant_adj_with_chem(n2u.A,n2u.C,lrd,left,right)
    nadjg = get_discrepant_adj_with_gap(n2u.A,n2u.E,lrd,left,right)

    jadjc = get_discrepant_adj_with_chem(jsh.A,jsh.C,lrd,left,right)
    jadjg = get_discrepant_adj_with_gap(jsh.A,jsh.E,lrd,left,right)

    badjc = get_discrepant_adj_with_chem_dev(n2u.A,jsh.A,n2u.C)
    badjg = get_discrepant_adj_with_gap_dev(n2u.A,jsh.A,n2u.E)


    data1 = [[ngap[0],nsyn[0]],
             [jgap[0],jsyn[0]],
             [bgap[0],bsyn[0]]]

    data2 = [[nadjg[0],nadjc[0]],
             [jadjg[0],jadjc[0]],
             [badjg[0],badjc[0]]]


    fig,ax = plt.subplots(1,1,figsize=(10,10))
    make_plot(ax,data1,title=("Fraction of discrepant synaptic connections\n"
                              "with a discrepant adjacency connection"))
    plt.tight_layout()
    if fout: plt.savefig(fout.replace('.','_1.'))

    fig,ax = plt.subplots(1,1,figsize=(10,10))
    make_plot(ax,data2,title=("Fraction of discrepant adjacency connections\n"
                              "with a synaptic connection"))
    plt.tight_layout()
    if fout: plt.savefig(fout.replace('.','_2.'))
    plt.show()

if __name__ == '__main__':
    run()
