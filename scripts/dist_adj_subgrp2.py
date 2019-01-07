"""
dist_adj_subgrp2.py

Plots adjacency degree distributions broken down by anatomical groups

created: Christopher Brittin
date: 01 November 2018 

"""

import sys
sys.path.append(r'./volumetric_analysis')
import matplotlib.pyplot as plt
import matplotlib as mpl
#Brittin modules
from connectome.load import from_db
from networks.stats import *
from figures.stats import *
import aux


mpl.rcParams['xtick.labelsize'] = 32
mpl.rcParams['ytick.labelsize'] = 32

neuron_class = './mat/nerve_ring_classes.txt'

def get_group_index(A,grp):
    idx = []
    if isinstance(grp,list):
        for g in grp:
            idx += [v.index for v in A.vs.select(group=g)]
    else:
        idx = [v.index for v in A.vs.select(group=grp)]
    return idx

def write_degrees(db,_neuron_class,fout):
    _remove = ['VC01','VD01','VB01','VB02']
    nclass = aux.read.into_dict(_neuron_class)

    C = from_db(db,adjacency=True,remove=_remove)
    C.A.assign_membership_dict(nclass,key='group')

    data = []
    for v in C.A.vs():
        data.append([v['name'],v['group'],C.A.degree(v)])

    aux.write.from_list(fout,data)

def group_degrees(db,_neuron_class):
    _remove = ['VC01','VD01','VB01','VB02']
    nclass = aux.read.into_dict(_neuron_class)

    C = from_db(db,adjacency=True,remove=_remove)
    C.A.assign_membership_dict(nclass,key='group')
    
    sp_idx = get_group_index(C.A,['Sp1','Sp2'])
    i1_idx = get_group_index(C.A,'I1')
    i2_idx = get_group_index(C.A,'I2')
    sa_idx = get_group_index(C.A,'Sa')
    smn_idx = get_group_index(C.A,'SMN')
    hmnp_idx = get_group_index(C.A,'HMNp')
    hmna_idx = get_group_index(C.A,'HMNa')

    degrees = [C.A.degree(sp_idx),
               C.A.degree(sa_idx),
               C.A.degree(i1_idx),
               C.A.degree(i2_idx),
               C.A.degree(smn_idx),
               C.A.degree(hmnp_idx),
               C.A.degree(hmna_idx)]

    return degrees

def run(fout=None,source_data=None):
    if source_data:
        fsplit = source_data.split('.')
        l4out = fsplit[0] + '_l4.' + fsplit[1]
        adultout = fsplit[0] + '_adult.' + fsplit[1]
        write_degrees('N2U',neuron_class,l4out)
        write_degrees('JSH',neuron_class,adultout)
        
    n2u = group_degrees('N2U',neuron_class)
    jsh = group_degrees('JSH',neuron_class)
    
    data = []
    for i in range(len(n2u)):
        data.append(n2u[i])
        data.append(jsh[i])

    fig,ax = plt.subplots(1,1,figsize=(15,10))
    dist_adj_subgroups2(ax,data,fout=fout)
    ax.xaxis.set_tick_params(labelsize=20) 
    plt.show()
        

if __name__ == '__main__':
    run()
