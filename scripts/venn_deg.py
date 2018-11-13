"""
venn_deg.py

Venn diagrams showing the reproducibility of connections between the
adult and L4 data sets for chemical synapses (left) gap junctions (center) 
and adjacency (right) connections.(b) Distribution of synaptic degree 
differences for homologous neurons.

created: Christopher Brittin
date: 01 November 2018 

"""
import sys
sys.path.append(r'./volumetric_analysis')
import matplotlib.pyplot as plt
import matplotlib as mpl
import networkx as nx
from matplotlib_venn import venn2
from matplotlib_venn import venn2

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

def get_venn_data(G1,G2):
    g1only = 0
    both = 0
    for (u,v) in G1.edges():
        if G2.has_edge(u,v):
            both += 1
        else:
            g1only += 1
    g2only = G2.number_of_edges() - both
    return [g1only,g2only,both]
    

_remove = ['VC01','VD01','VB01','VB02']

def run(fout=None):
    n2u = from_db('N2U',adjacency=True,chemical=True,electrical=True,
                  dataType='networkx',remove=_remove)

    jsh = from_db('JSH',adjacency=True,chemical=True,electrical=True,
                  dataType='networkx',remove=_remove)

    print(n2u.A.number_of_nodes(),jsh.A.number_of_nodes())
    print(n2u.A.number_of_edges() - jsh.A.number_of_edges())
    print(n2u.C.number_of_edges() - jsh.C.number_of_edges())
    print(n2u.E.number_of_edges() - jsh.E.number_of_edges())
    print(set(n2u.A.nodes()) - set(jsh.A.nodes()))

    syn = get_venn_data(n2u.C,jsh.C)
    gap = get_venn_data(n2u.E,jsh.E)
    adj = get_venn_data(n2u.A,jsh.A)
    print(syn,adj)

    fig,ax = plt.subplots(1,3,figsize=(12,5))
    
    v = venn2(subsets=syn,ax=ax[0],set_labels=('Adult','L4'))
    v.get_patch_by_id('100').set_color(ADULT_COL)
    v.get_patch_by_id('010').set_color(L4_COL)
    v.get_patch_by_id('11').set_color(AL_COL)
    for text in v.set_labels: text.set_fontsize(24)
    for text in v.subset_labels: text.set_fontsize(14)
    ax[0].set_title('# of chemical\nsynaptic connections',fontsize=24)
    v = venn2(subsets=gap,ax=ax[1],set_labels=('Adult','L4'))
    v.get_patch_by_id('100').set_color(ADULT_COL)
    v.get_patch_by_id('010').set_color(L4_COL)
    v.get_patch_by_id('11').set_color(AL_COL)
    for text in v.set_labels: text.set_fontsize(24)
    for text in v.subset_labels: text.set_fontsize(14)
    ax[1].set_title('# of gap junction\nconnections',fontsize=24)
    v = venn2(subsets=adj,ax=ax[2],set_labels=('Adult','L4'))
    v.get_patch_by_id('100').set_color(ADULT_COL)
    v.get_patch_by_id('010').set_color(L4_COL)
    v.get_patch_by_id('11').set_color(AL_COL)
    for text in v.set_labels: text.set_fontsize(24)
    for text in v.subset_labels: text.set_fontsize(14)
    ax[2].set_title('# of adjacency\nconnections',fontsize=24)
    if fout: plt.savefig(fout)
    plt.show()

if __name__=="__main__":
    run()
