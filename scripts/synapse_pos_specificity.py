"""
synaptic_specificity.py

Distribution of differences between homologous mean synapse positions.

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

lr_pairs = './mat/lr_neurons.txt'
lr_dict = './mat/lr_dict.txt'
homologs = './mat/homologs.txt'
left_nodes = './mat/left_nodes.txt'
right_nodes = './mat/right_nodes.txt'

def format_subcell(S,D,thresh=0.05):
    dgap,mgap = [],[]
    dpre,mpre = [],[]
    dpost,mpost = [],[]
    for n in S:
        if S[n][0] <= thresh and n in D:
            dgap += D[n][0]
            #mgap += D[n][0][1]
        if S[n][1] <= thresh and n in D:
            dpre += D[n][1]
            #mpre += D[n][1][1]
        if S[n][2] <= thresh and n in D:
            dpost += D[n][2]
            #mpost += D[n][2][1]

    
    # data =  [dgap,mgap,dpre,mpre,dpost,mpost]
    data = [dgap,dpre,dpost]
    #for i in xrange(len(data)):
    #    data[i] = [d for d in data[i] if d <= 1]
    return data

def run(fout=None):
    N2U = 'N2U'
    JSH = 'JSH'
    _remove = ['VC01','VD01','VB01','VB02']
    neurons = aux.read.into_list2(lr_pairs)
    lrd = aux.read.into_lr_dict(lr_dict)
    left = aux.read.into_list(left_nodes)
    left.remove('CEHDL')
    left.remove('CEHVL')
    left.remove('HSNL')
    left.remove('PVNL')
    left.remove('PLNL')

        
    N2U = from_db(N2U,adjacency=True,chemical=True,
                  electrical=True,remove=_remove,dataType='networkx')
    JSH = from_db(JSH,adjacency=True,chemical=True,
                  electrical=True,remove=_remove,dataType='networkx')

    n2ucon = db.connect.default('N2U')
    n2ucur = n2ucon.cursor()
    jshcon = db.connect.default('JSH')
    jshcur = jshcon.cursor()    
    
    both_nodes = set(N2U.A.nodes()) & set(JSH.A.nodes())
    both_nodes.remove('SABD')
    both_nodes.remove('FLPL')
    both_nodes.remove('FLPR')
    if 'VD01' in both_nodes: both_nodes.remove('VD01')
    
    S1 = synspec.get_bilateral_specificity(N2U,lrd,left)
    D1 = synspec.get_bilateral_subcell_specificity(n2ucur,neurons,lrd)
    B1 = format_subcell(S1,D1)

    S2 = synspec.get_bilateral_specificity(JSH,lrd,left)
    D2 = synspec.get_bilateral_subcell_specificity(jshcur,neurons,lrd)
    B2 = format_subcell(S2,D2)
        
    S3 = synspec.get_developmental_specificity(N2U,JSH,
                                       both_nodes=both_nodes)
    D3 = synspec.get_developmental_subcell_specificity(n2ucur,
                                                       jshcur,
                                                       both_nodes=both_nodes)
    B3 = format_subcell(S3,D3)
    n2ucon.close()
    jshcon.close()

    labels = None
    pos = [1.5,2,2.5,3.5,4,4.5,5.5,6,6.5]
    data = [B1[0],B2[0],B3[0],
            B1[1],B2[1],B3[1],
             B1[2],B2[2],B3[2]]

    colors = [ADULT_COL,L4_COL,AL_COL,
              ADULT_COL,L4_COL,AL_COL,
              ADULT_COL,L4_COL,AL_COL,]
    fig,ax = plt.subplots(1,1,figsize=(12,10))
    bp = fstats.plot_boxplots(ax,data,labels=labels,positions=pos,
                              ylim=[-5,5],
                              ylabel='Mean position difference',
                              title='Mean synapse position',
                              showfliers=True,width=0.2,colors=colors)

    ax.set_xticklabels(['gap j.',
                        'presyn.',
                        'postsyn.'])
    ax.set_xticks([2, 4, 6])
    ax.set_ylim([-1,1])
    ax.axvspan(0,3,facecolor='#C3C3C3')
    ax.axvspan(3,5,facecolor='#D8D7D7')
    ax.axvspan(5,8,facecolor='#C3C3C3')
    ax.axhline(0,color='r',linewidth=3,linestyle='--')
        
    _A, = ax.plot([1,1],ADULT_COL)
    _L, = ax.plot([1,1],L4_COL)
    _AL, = ax.plot([1,1],AL_COL)
    leg =ax.legend((_A, _L,_AL),
                   ('Adult L/R', 'L4 L/R','Adult/L4'),
                   fontsize=18)
    for legobj in leg.legendHandles:
        legobj.set_linewidth(4.0)
    _A.set_visible(False)
    _L.set_visible(False)
    _AL.set_visible(False)

    plt.tight_layout()
    if fout: plt.savefig(fout)
    plt.show()


    
if __name__=='__main__':
    run()
