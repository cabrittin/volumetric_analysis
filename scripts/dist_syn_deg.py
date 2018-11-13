"""
dist_syn_deg.py

Adult postsynaptic, presynaptic and gap junction degree
distributions. Fit with a kernal density estimator.

created: Christopher Brittin
date: 01 November 2018 

"""
import sys
sys.path.append(r'./volumetric_analysis')
import numpy as np
import matplotlib.pyplot as plt

from connectome.load import from_db
from figures.stats import *

db = 'N2U'
remove = ['VC01','VD01','VB01','VB02']

def run(fout=None):
    C = from_db(db,adjacency=True,chemical=True,electrical=True,remove=remove)
    C.reduce_to_adjacency()
    degree_in = C.C.degree(mode='in')
    degree_out = C.C.degree(mode='out')
    degree_gap = C.E.degree()

    fig,ax = plt.subplots(1,3,figsize=(15,5),sharey=True)
    plot_syn_degree(ax[0],degree_out,density=True,fit_mode='KDE',ylim=[0,0.15],
                    xlabel='Postsynaptic degree, $d^\mathrm{out}$',
                    ylabel='Probability',xfs=24,yfs=24)
    plot_syn_degree(ax[1],degree_in,density=True,fit_mode='KDE',ylim=[0,0.15],
                    xlabel="Presynaptic degree, $d^\mathrm{in}$",xfs=24)
    plot_syn_degree(ax[2],degree_gap,density=True,fit_mode='KDE',ylim=[0,0.15],
                    xlabel='Gap junction degree, $d^\mathrm{gap}$',xfs=24)
    if fout: plt.savefig(fout)
    plt.show()

if __name__=="__main__":
    run()
