import matplotlib.pyplot as plt
import matplotlib as mpl
#Brittin modules
from Connectome.load import from_db
from Networks.stats import *
from Figures.stats import *
import LUT
import aux


mpl.rcParams['xtick.labelsize'] = 32
mpl.rcParams['ytick.labelsize'] = 32

fout = '/home/cabrittin/Dropbox/PhD/sr_vol/figs2/dist_adj_subgrp2.png'
neuron_class = '../Data/brittin_classv3.txt'

def get_group_index(A,grp):
    idx = []
    if isinstance(grp,list):
        for g in grp:
            idx += [v.index for v in A.vs.select(group=g)]
    else:
        idx = [v.index for v in A.vs.select(group=grp)]
    return idx

def group_degrees(db):
    _remove = ['VC01','VD01','VB01','VB02']
    loader = LUT.loader(db)
    nclass = aux.read.into_dict(neuron_class)

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

n2u = group_degrees('N2U')
jsh = group_degrees('JSH')

data = []
for i in range(len(n2u)):
    data.append(n2u[i])
    data.append(jsh[i])

fig,ax = plt.subplots(1,1,figsize=(12,10))
dist_adj_subgroups2(ax,data,fout=fout)
plt.show()
"""
db = 'JSH'

C = from_db(db,adjacency=True,remove=_remove)
C.A.assign_membership_dict(nclass,key='group')

s_idx = [v.index for v in C.A.vs.select(group='S')]
i_idx = [v.index for v in C.A.vs.select(group='I')]
m_idx = [v.index for v in C.A.vs.select(group='M')]

data.append(C.A.degree(s_idx))
data.append(C.A.degree(i_idx))
data.append(C.A.degree(m_idx))



fig,ax = plt.subplots(1,1,figsize=(12,10))
dist_adj_subgroups(ax,data,fout=fout,)

plt.tight_layout()
plt.savefig(fout)
plt.show()
"""

