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
from numpy.polynomial.polynomial import polyfit
from scipy.stats import pearsonr


from connectome.load import from_db
import aux
from figures.stats import *


def plot_cor(ax,A):
    xmu = np.mean(A[:,0])
    ymu = np.mean(A[:,1])
    A[:,0] = (A[:,0] - xmu) / xmu
    A[:,1] = (A[:,1] - ymu) / ymu
    ax.plot(A[:,0],A[:,1],'bo')
    ax.plot(np.unique(A[:,0]), np.poly1d(np.polyfit(A[:,0], A[:,1], 1))(np.unique(A[:,0])),'r-')
    r =pearsonr(A[:,0],A[:,1])
    ax.text(0.05,0.9,r'$n$ = %d'%A.shape[0],
            transform=ax.transAxes,fontsize=24)
    ax.text(0.05,0.83,r'$r^2$ = %1.2f' %r[0]**2,
            transform=ax.transAxes,fontsize=24)



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
    
    parser.add_argument('-o','--output',
                        dest = 'fout',
                        action="store",
                        required= False,
                        default = None,
                        help="Output file, if you wish to save the output."
                        )
    
    parser.add_argument('--title',
                        dest = 'title',
                        required=False,
                        default=None,
                        help="Plot title")

    params = parser.parse_args()
    
    stats = aux.read.into_list2(params.stats)
    _remove = ['VC01','VD01','VB01','VB02']
    C = from_db(params.db,adjacency=True,chemical=True,electrical=True,remove=_remove,dataType='networkx')
  
    A = np.zeros((len(stats),2))
    i  = 0
    for d in stats:
        if not C.A.has_node(d[0]): continue
        A[i][1] = C.A.degree(d[0])
        A[i][0] = int(d[1])
        i += 1
   
    
    fig,ax = plt.subplots(1,1,figsize=(15,10))
    plot_cor(ax,A)
    ax.set_xlim([-2,2])
    ax.set_ylim([-2,2])
    ax.set_xlabel('Normalized cell surface area',fontsize=32)
    ax.set_ylabel('Normalized cell adjacency degree',fontsize=32)
    if params.title: ax.set_title(params.title,fontsize=36)
    if params.fout: plt.savefig(params.fout)
    plt.show()


