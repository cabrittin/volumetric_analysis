"""
synapse_pos_correlation.py

Distribution of differences between homologous mean synapse positions.
Correlation of
mean synapse position for homologous gap junctions, presynaptic and postsynaptic contacts for neurons AIZL and AIZR in the adult.

Author: Christopher Brittin
Created: 07 February 2018

"""

import sys
sys.path.append(r'./volumetric_analysis/')
import matplotlib.pyplot as plt
import matplotlib as mpl

import db
from connectome.load import from_db
import connectome.synspecificity as synspec
import figures.stats as fstats
import aux

mpl.rcParams['xtick.labelsize'] = 24
mpl.rcParams['ytick.labelsize'] = 24

ADULT_COL = '#FFC300'
L4_COL = '#3380FF'
AL_COL = '#14FA29' 

PRESYN_COL = '#FA5252'
POSTSYN_COL = '#F9E530'
GAPJUNC_COL = '#47A705'

_db = 'N2U'
lr_pairs = './mat/lr_neurons.txt'
neuron1 = 'AIZL'
neuron2 = 'AIZR'


def run(fout=None):
    con = db.connect.default(_db)
    cur = con.cursor()

    lrd = aux.read.into_lr_dict(lr_pairs)
    
    S1 = synspec.synapse_positions(cur,neuron1)
    S2 = synspec.synapse_positions(cur,neuron2)
    S1,S2 = synspec.format_left_right_subcell(S1,S2,lrd)
    synspec.make_subcell_table(S1,S2)

    fig,ax = plt.subplots(1,1,figsize=(12,9))
    fstats.plot_corr(ax,S1['gap'],S2['gap'],
                    corr_label='$r^2_{\mathrm{gap}}$',
                    label='gap junction',col=GAPJUNC_COL,marker='^')
    fstats.plot_corr(ax,S1['pre'],S2['pre'],corr_label='$r^2_{\mathrm{pre}}$',
                    label='presynaptic',col=PRESYN_COL,_y=0.8)
    fstats.plot_corr(ax,S1['post'],S2['post'],_y=0.7,
                    corr_label='$r^2_{\mathrm{post}}$',
                    label='postsynaptic',col=POSTSYN_COL,marker='s')   
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])
    ax.set_xlabel('%s, mean synapse distance'%neuron1,fontsize=28)
    ax.set_ylabel('%s, mean synapse distance'%neuron2,fontsize=28)
    ax.legend(loc='lower right',fontsize=24,markerscale=2)
    plt.tight_layout()           
    if fout: plt.savefig(fout)
    plt.show()

if __name__=="__main__":
    run()
