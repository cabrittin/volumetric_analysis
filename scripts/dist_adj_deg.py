"""
dist_adj_deg.py

Adult adjacency degree distribution. Fit with a 3 component Gaussian 
Mixture Model.

created: Christopher Brittin
date: 01 November 2018 

"""

import sys
sys.path.append(r'./volumetric_analysis')
import numpy as np
import matplotlib.pyplot as plt

from connectome.load import from_db
from figures.stats import *



def run(fout=None):
    db = 'N2U'
    remove = ['VC01','VD01','VB01','VB02']
    C = from_db(db,adjacency=True,remove=remove)
    degree = C.A.degree()

    fig,ax = plt.subplots(1,1,figsize=(7,5))
    plot_adj_degree(ax,degree,density=True,fit_mode='GMM',ylim=[0,0.1],
                    xlabel='Adjacency degree, $d$',
                    ylabel='Probability',yfs=24,xfs=24)
    if fout: plt.savefig(fout)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    run()
