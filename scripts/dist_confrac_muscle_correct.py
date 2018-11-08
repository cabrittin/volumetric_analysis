"""
dist_confrac_muscle_correct.py

Connectivity fractions corrected for muscles. A significant fraction 
of motor neuron output and some sensory/interneuron output is onto 
muscle which is not included in our segmentation. To correct for this, 
we artifically added muscles to the adjacency list. If neuron A synapses 
onto muscle M, then the adjacency (A, M) was added to the adjacency list.
Hence, every adjacency added has a corresponding synaptic connection. 
Therefore, the corrected connectivity fractions are likely overestimates 
because there are adjacencies with muscles that do not result in synapses.


created: Christopher Brittin
date: 01 November 2018

"""

import sys
sys.path.append(r'./volumetric_analysis')
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib as mpl

from connectome.load import from_db
import aux
import db
from networks.stats import *
from figures.stats import *


mpl.rcParams['xtick.labelsize'] = 32
mpl.rcParams['ytick.labelsize'] = 32

_nclass = "./mat/simmu_cell_class.txt"

def get_cf(C,vertices=None):
    data =  {'pre':get_cpre(C,vertices=vertices),
             'post':get_cpost(C,vertices=vertices),
             'gap':get_cgap(C,vertices=vertices)}
    return data

def get_cpre(C,vertices=None):
    if vertices:
        v = vertices
    else:
        v = C.A.nodes()
    con = dict(C.C.out_degree())
    adj = dict(C.A.degree())
    return dict([(n,con[n] / float(adj[n])) for n in adj if n in con])
    
def get_cpost(C,vertices=None):
    if vertices:
        v = vertices
    else:
        v = C.A.nodes()
    con = dict(C.C.in_degree())
    adj = dict(C.A.degree())
    return dict([(n,con[n] / float(adj[n])) for n in adj if n in con])

def get_cgap(C,vertices=None):
    if vertices:
        v = vertices
    else:
        v = C.A.nodes()
    con = dict(C.E.degree())
    adj = dict(C.A.degree())
    return dict([(n,con[n] / float(adj[n])) for n in adj if n in con])

def run(fout=None):
    _db = 'N2U'
    con = db.connect.default(_db)
    cur = con.cursor()

    muscle = ['dBWML1', 'dBWML2', 'dBWML3', 'dBWML5',
              'dBWML6', 'dBWML7', 'dBWMR1', 'dBWMR2',
              'dBWMR3', 'dBWMR4', 'dBWMR5', 'dBWMR6',
              'dBWMR7', 'exc_gl', 'hyp', 'pharyngealepithelium',
              'seamL', 'vBWML1', 'vBWML2', 'vBWML4',
              'vBWML5', 'vBWML6', 'vBWML7', 'vBWML8',
              'vBWMR1', 'vBWMR2', 'vBWMR3', 'vBWMR4',
              'vBWMR5', 'vBWMR6', 'vBWMR7', 'vBWMR8']

    remove = ['VA02','VC01','VD01','VB01','VB02']
    nclass = aux.read.into_dict(_nclass)
    C = from_db(_db,
                chemical=True,electrical=True,
                remove=remove,dataType='networkx')

    neurons = sorted(db.mine.get_adjacency_cells(cur))
    adjacency = db.mine.get_adjacency_data(cur)
    C.load_adjacency(adjacency,directed=False)


    for (u,v) in C.C.edges():
        if v in muscle:
            C.A.add_edge(u,v)

    for (u,v) in C.E.edges():
        if u in muscle or v in muscle:
            C.A.add_edge(u,v)

    nx.set_node_attributes(C.C,nclass,'name')
    nx.set_node_attributes(C.E,nclass,'name')
    nx.set_node_attributes(C.A,nclass,'name')

    for n in C.A.nodes():
        if 'name' not in C.A.node[n]:
            C.A.node[n]['name'] = 'Mu'

    data = []
    stype = ['gap','pre','post']
    s_idx = [n for n in C.A.nodes if C.A.node[n]['name'] == 'S']
    i_idx = [n for n in C.A.nodes if C.A.node[n]['name'] == 'I']
    m_idx = [n for n in C.A.nodes if C.A.node[n]['name'] == 'M']

    cf = get_cf(C)

    for s in stype:
        #print(sorted(cf[s].keys()))
        data.append(list(cf[s][n] for n in s_idx if n in cf[s]))
        data.append(list(cf[s][n] for n in i_idx if n in cf[s]))
        data.append(list(cf[s][n] for n in m_idx if n in cf[s]))

    #print(data)
    fig,ax = plt.subplots(1,1,figsize=(12,10))
    plot_confrac_subgroups(ax,data,fout=fout)

    plt.tight_layout()
    if fout: plt.savefig(fout)
    plt.show()

if __name__ == "__main__":
    run()
