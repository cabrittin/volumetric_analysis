"""
modularity.py

Plots for modularity analysis

"""

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

mpl.rcParams['xtick.labelsize'] = 24
mpl.rcParams['ytick.labelsize'] = 24

from Connectome.format_igraph import *

ADULT_COL = '#FFC300'
L4_COL = '#3380FF'

def adjacency_matrix(fig,ax,vertex_order,G,fout=None):
    Z = G.get_numpy_array(vertex_order=vertex_order,edge_attr='weight')
    Z *= 5*90*(1e-6)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    cax = ax.pcolor(Z,vmin=0,vmax=2)
    cbar = fig.colorbar(cax)
    cbar.ax.get_yaxis().labelpad = 50
    cbar.set_label('Membrane contact, $\mu$m$^2$', rotation=270,fontsize=24)
    ax.set_title('Adjacency modularity',fontsize=32)
    plt.tight_layout()
    if fout: plt.savefig()

def synaptic_matrix(fig,ax,vertex_order,G,fout=None):
    Z = G.get_numpy_array(vertex_order=vertex_order,edge_attr='weight')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    cax = ax.pcolor(Z,vmin=0,vmax=10)
    cbar = fig.colorbar(cax)
    cbar.ax.get_yaxis().labelpad = 50
    cbar.set_label('#EM sections', rotation=270,fontsize=24)
    ax.set_title('Synaptic modularity',fontsize=32)
    plt.tight_layout()
    if fout: plt.savefig()


def adj_syn_norm_info(ax,n2u,jsh,fout=None):
    N = 3
    ind = np.arange(N)
    width = 0.35
    rects1 = ax.bar(ind,n2u,width,color=ADULT_COL)
    rects2 = ax.bar(ind+width,jsh,width,color=L4_COL)

    ax.set_ylabel('Normalized mutual information',fontsize=38)
    ax.set_title('Adjacency vs. synaptic modularity',fontsize=38)
    ax.set_xticks(ind + width / 2)
    ax.set_xticklabels(('Chem.','Elec.','Chem. + Elec.'))
    ax.set_ylim([0,1.0])
    ax.legend((rects1[0],rects2[0]),('Adult','L4'),fontsize=32)
    if fout: plt.savefig()


def mod_bar_graph(ax,n2u,jsh,width=0.15,N=3):
    ind = np.arange(N)

    rects1 = ax.bar(ind,n2u,width,color=ADULT_COL,label='Adult')
    rects2 = ax.bar(ind+width,jsh,width,color=L4_COL,label='L4')    
    ax.set_xticks(ind + width / 2)
    ax.set_ylim([0,1.0])
    ax.set_xticklabels(('','',''))
