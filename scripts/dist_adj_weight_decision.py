"""
dist_adj_weight_decision.py

Probability density of the log of surface area contacts for 
adjacencies that do (+) and do not (-) produce a synapse. 
Red line indicates a decision boundary. Right of the boundary 
an adjacency has a higher probability of producing a synapse.

created: Christopher Brittin
date: 01 November 2018

"""
import sys
sys.path.append(r'./volumetric_analysis')
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import digamma
import random

from connectome.load import from_db
from networks.stats import get_corresponding_edge_attr
from models.mutual_info import *
from figures.stats import *

db = 'N2U'
remove = ['VC01','VD01','VB01','VB02']
SCALE = 5*90*(1e-6)
KMAX = 100
krange = (3,KMAX)
THETA = 0
nbins=100

def get_data(G1,G2):
    #Get edge weights
    N = G2.ecount()
    data = np.zeros((N,2))
    for i in range(N):
        e = G2.es[i]
        data[i,0] = e['weight']
        if G1.are_connected(e.source,e.target):
            w = G1.es[G1.get_eid(e.source,e.target)]['weight']
            if w > THETA:
                data[i,1] = 1
        
    data[:,0] = np.log(data[:,0]*SCALE)
    return data
    

def run(fout=None):
    C = from_db(db,adjacency=True,chemical=True,electrical=True,remove=remove)
    C.C.reduce_to(C.A)
    C.E.reduce_to(C.A)
    N = C.C.ecount()

    C.C.to_undirected(combine_edges=sum)

    data = get_data(C.C,C.A)
    data = data[data[:,0].argsort()]
    pos = np.where(data[:,1] == 1)[0]
    neg = np.where(data[:,1] == 0)[0]


    h1,bins1 = np.histogram(data[neg,0],bins=nbins,
                            range=(-4,4),normed=True,density=True)
    dx = bins1[1] - bins1[0]
    h1 = np.cumsum(h1)*dx
    h2,bins2 = np.histogram(data[pos,0],bins=nbins,
                            range=(-4,4),normed=True,density=True)
    dx = bins2[1] - bins2[0]
    h2 = np.cumsum(h2[::-1])[::-1]*dx
    imin = np.argmin(abs(h1-h2)) + 1

    fig,ax = plt.subplots(1,1,figsize=(15,10))
    n,bins,_ = ax.hist(data[neg,0],bins=nbins,range=(-4,4),histtype='step',
                       density=True,cumulative=False,linewidth=4,
                       label='(-) synapse')
    n,bins,_ = ax.hist(data[pos,0],bins=nbins,range=(-4,4),histtype='step',
                       density=True,cumulative=False,linewidth=4,
                       label='(+) synapse')
    ax.axvline(bins[imin],color='r',linewidth=3)
    #print(bins[imin])
    #plot_skew_norm_fit(ax,data[neg,0],bins)
    #plot_skew_norm_fit(ax,data[pos,0],bins)
    plt.legend(loc='upper left',fontsize=24)
    ax.set_xlim([-4,4])
    ax.set_xlabel('log(surface area contact)',fontsize=32)
    ax.set_ylabel('Probability density',fontsize=32)
    plt.tight_layout()
    if fout: plt.savefig(fout)
    plt.show()


if __name__ == '__main__':
    run()
