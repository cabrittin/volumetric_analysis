"""
dist_confrac_subgrp.py

Plot distribution of connectivity fractions for functional neuron classes.

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

_nclass = "./mat/simmu_cell_class.txt"

def run(fout=None):
    db = 'N2U'
    _remove = ['VC01','VD01','VB01','VB02']
    nclass = aux.read.into_dict(_nclass)

    C = from_db(db,adjacency=True,chemical=True,electrical=True,remove=_remove)
    C.C.assign_membership_dict(nclass,key='group')

    for e in C.C.es:
        if e['weight'] > 50:
            print(e['weight'],C.C.vs[e.source]['name'],C.C.vs[e.target]['name'])


    print(C.C.ecount(),C.E.ecount())
    cf = get_cf(C)
    
    data = []
    stype = ['gap','pre','post']
    grp = ['S','I','M']

    s_idx = [v.index for v in C.C.vs.select(group='S')]
    i_idx = [v.index for v in C.C.vs.select(group='I')]
    m_idx = [v.index for v in C.C.vs.select(group='M')]

    for s in stype:
        data.append(cf[s][s_idx])
        data.append(cf[s][i_idx])
        data.append(cf[s][m_idx])

    fig,ax = plt.subplots(1,1,figsize=(12,10))
    plot_confrac_subgroups(ax,data,fout=fout)

    plt.tight_layout()
    if fout: plt.savefig(fout)
    plt.show()

    """
    db = 'JSH'
    
    _remove = ['VC01','VD01','VB01','VB02']
    #loader = LUT.loader(db)
    nclass = aux.read.into_dict(loader['nclass'])

    C = from_db(db,adjacency=True,chemical=True,electrical=True,remove=_remove)
    C.C.assign_membership_dict(nclass,key='group')
    print(C.C.ecount(),C.E.ecount())
    """

if __name__ == '__main__':
    run()
