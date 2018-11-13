"""
dist_adj_weight.py

Cumulative distribution of adult surface area contacts between cells. 
Inset: Normal probability plot of the log of surface area contacts. 
The middle ∼ 95% of data is lognormal.

created: Christopher Brittin
date: 01 November 2018 

"""
import sys
sys.path.append(r'./volumetric_analysis')
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


from connectome.load import from_db
from figures.stats import *

SCALE = 5*90*(1e-6)
db = 'N2U'
remove = ['VC01','VD01','VB01','VB02']

def run(fout=None):
    C = from_db(db,adjacency=True,remove=remove)
    td =  np.array(C.A.es['weight'])*SCALE


    fig,ax = plt.subplots(1,1,figsize=(15,10))
    plot_dist(ax,td,density=True,cumulative=True,xlim=[0,2],ylim=[0,1],
              hrange=(0,2),nbins=1000,
              xlabel='Surface area contact ($\mu$m$^2$)',
              ylabel='Cumulative distribution',fs=24)
    axins = inset_axes(ax,
                       width="50%",  # width = 30% of parent_bbox
                       height="50%",  # height : 1 inch
                       loc=5)
    plot_lognorm_probplot(axins,td,fout=None,fs=24,
                          ylabel='Ordered log($w$)')
    plt.tight_layout()
    if fout: plt.savefig(fout)
    plt.show()

if __name__=="__main__":
    run()
