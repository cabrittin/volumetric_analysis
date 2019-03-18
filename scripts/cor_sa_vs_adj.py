"""
cor_sa_vs_adj.py

Plots the correlation of surface area vs adjacency degree

@author Christopher Brittin
@date 2019-03-17
"""

import sys
sys.path.append('./volumetric_analysis')
import argparse
import numpy as np
import matplotlib.pyplot as plt

from connectome.load import from_db
import aux
from figures.stats import *


def plot_cor(ax,A):
    xmu = np.mean(A[:,0])
    ymu = np.mean(A[:,1])
    A[:,0] = (A[:,0] - xmu) / xmu
    A[:,1] = (A[:,1] - ymu) / ymu
    ax.plot(A[:,0],A[:,1],'bo')
    ax.set_xlim([-2,2])
    ax.set_ylim([-2,2])

SCALE = 450e-6

if __name__== '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('db',
                        action="store",
                        help="Database")

    parser.add_argument('stats',
                        action="store",
                        help="Path to stats file")
    
    parser.add_argument('nclass',
                        action="store",
                        help="Path the class file")

    parser.add_argument('-o','--output',
                        dest = 'fout',
                        action="store",
                        required= False,
                        default = None,
                        help="Output file, if you wish to save the output."
                        )

    params = parser.parse_args()
    
    stats = aux.read.into_list2(params.stats)
    nclass = aux.read.into_dict(params.nclass)
    _remove = ['VC01','VD01','VB01','VB02']
    C = from_db(params.db,adjacency=True,chemical=True,electrical=True,remove=_remove,dataType='networkx')
  
    sa = {}
    A = np.zeros((len(stats),2))
    i  = 0
    for d in stats:
        if not C.A.has_node(d[0]): continue
        nc = nclass[d[0]]
        if nc not in sa: sa[nc] = []
        A[i][1] = C.A.degree(d[0])
        A[i][0] = int(d[1])
        sa[nc].append(int(d[1])*SCALE)
        i += 1
   
    SA = [sa['Sp1'] + sa['Sp2'],
            sa['Sa'],sa['I1'],sa['I2'],sa['SMN'],
            sa['HMNp'],sa['HMNa']]
    
    data = []
    for i in range(len(SA)):
        data.append(SA[i])
        data.append(SA[i])
    fig,ax = plt.subplots(1,1,figsize=(20,10))
    #plot_cor(ax[0],A)
    dist_adj_subgroups2(ax,data)
    ax.set_ylim([0,500])    
    plt.show()


