"""
synaptic_specificity.py

Synapse probabilities are either below (+) or above (-) the 0.05 threshold. 
Left shows the different +/- combinations for p_s^gap , p_s^pre and p_s^post. 
Bar plots show fractions of neurons with the indicated combination of 
specificity probabilities. Bar color indicates either adult (yellow) or 
L4  (blue) left/right comparison or comparison between adult and L4 (green).

Author: Christopher Brittin
Created: 07 February 2018

"""

import sys
sys.path.append(r'./volumetric_analysis/')
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as patches
from matplotlib import lines
import numpy as np


from connectome.load import from_db
import connectome.synspecificity as synspec
import aux

ADULT_COL = '#FFC300'
L4_COL = '#3380FF'
AL_COL = '#14FA29' 

lr_dict = './mat/lr_dict.txt'
homologs = './mat/homologs.txt'
left_nodes = './mat/left_nodes.txt'
right_nodes = './mat/right_nodes.txt'

def get_bilateral_specificity(C):
    db = C.db
    lrd = aux.read.into_lr_dict(lr_dict)
    nclass = aux.read.into_list2(homologs)
    left = aux.read.into_list(left_nodes)
    
    ndict = {}    
    for n in nclass:
        for _n in n[1:]:ndict[_n] = n[0]
    
    left.remove('CEHDL')
    left.remove('CEHVL')
    left.remove('HSNL')
    left.remove('PVNL')
    left.remove('PLNL')
    _pre = synspec.bilateral_specificity(C,left,lrd)
    _post = synspec.bilateral_specificity(C,left,lrd,mode='post')
    _gap = synspec.bilateral_specificity(C,left,lrd,mode='gap') 

    spec = {}
    for nl in left:
        spec[nl] = [_gap[nl].p,_pre[nl].p,_post[nl].p]
    return spec

def get_developmental_specificity(C1,C2):
    both_nodes = set(C1.A.nodes()) & set(C2.A.nodes())
    both_nodes.remove('SABD')
    if 'VD01' in both_nodes: both_nodes.remove('VD01')    
    
    _pre = synspec.developmental_specificity(C1,C2,both_nodes)
    _post = synspec.developmental_specificity(C1,C2,both_nodes,mode='post')
    _gap = synspec.developmental_specificity(C1,C2,both_nodes,mode='gap')

    #for n in sorted(both_nodes):
    #    print(n,_pre[n].p,_post[n].p,_gap[n].p)
    
    spec = {}
    for nl in both_nodes:
        spec[nl] = [_gap[nl].p,_pre[nl].p,_post[nl].p]
    return spec   

def plot_specificity(ax,db1,db2,fout=None):
    S1 = get_bilateral_specificity(db1)
    S2 = get_bilateral_specificity(db2)
    S3 = get_developmental_specificity(db1,db2)

    B1,O1 = format_bar_specificity(S1)
    B2,O2 = format_bar_specificity(S2)
    B3,O3 = format_bar_specificity(S3)
    
    fs = 24
    width = 0.25
    y_pos = np.arange(8)
    ax.grid(zorder=0)
    ax.barh(y_pos-width,B1,width,label='Adult L/R',zorder=3,color=ADULT_COL)
    ax.barh(y_pos,B2,width,label='L4 L/R',zorder=3,color=L4_COL)
    ax.barh(y_pos+width,B3,width,label='Adult/L4',zorder=3,color=AL_COL)
    ax.set_xlim([0,0.5])
    ax.set_yticks(y_pos)
    ax.set_yticklabels(('','','','','','','',''))#gpp)
    ax.invert_yaxis()
    ax.text(0.01, 0.9, '$p^{\mathrm{gap}}_s$',rotation=45,ha='left',
            fontsize=fs, transform=plt.gcf().transFigure)
    ax.text(0.055, 0.9, '$p^{\mathrm{pre}}_s$',rotation=45,ha='left',
            fontsize=fs, transform=plt.gcf().transFigure)
    ax.text(0.09, 0.9, '$p^{\mathrm{post}}_s$',rotation=45,ha='left',
            fontsize=fs, transform=plt.gcf().transFigure)
    ptxt = [(0.01,0.8),(0.01,0.71),(0.01,0.62),(0.01,0.53),
            (0.055,0.8),(0.055,0.71),(0.055,0.44),(0.055,0.35),
            (0.09,0.8),(0.09,0.62),(0.09,0.44),(0.09,0.26)]
    for c in ptxt:
        ax.text(c[0],c[1],'+',fontsize=fs,transform=plt.gcf().transFigure)
    ntxt = [(0.015,0.44),(0.015,0.35),(0.015,0.26),(0.015,0.17),
            (0.06,0.62),(0.06,0.53),(0.06,0.26),(0.06,0.17),
            (0.095,0.71),(0.095,0.53),(0.095,0.35),(0.095,0.17)]
    for c in ntxt:
        ax.text(c[0],c[1],'-',fontsize=fs,transform=plt.gcf().transFigure) 
    ax2 = plt.axes([0,0,1,1], facecolor=(1,1,1,0))
    x,y = np.array([[0.045,0.045], [0.1,0.85]])
    line = lines.Line2D(x, y, lw=3., color='k')
    ax2.add_line(line)
    x,y = np.array([[0.085,0.085], [0.1,0.85]])
    line = lines.Line2D(x, y, lw=3., color='k')
    ax2.add_line(line)   
    ax.legend(loc='lower right',fontsize=fs)
    ax.set_xlabel('Fraction of neurons below specificity threshold ($p_s<0.05$)',
                  fontsize=32)
    ax.set_title('Synaptic specificity',fontsize=32)
    if fout: plt.savefig(fout)

    
def format_bar_specificity(S,thresh = 0.05):
    N = float(len(S))
    data = np.zeros(8)
    outliers = []
    for s in sorted(S):
        [gap,pre,post] = S[s]
        if gap <= thresh and pre <= thresh and post <= thresh:
            data[0] += 1
        elif gap <= thresh and pre <= thresh and post > thresh:
            data[1] += 1
        elif gap <= thresh and pre > thresh and post <= thresh:
            data[2] += 1
        elif gap <= thresh and pre > thresh and post > thresh:
            data[3] += 1
        elif gap > thresh and pre <= thresh and post <= thresh:
            data[4] += 1        
        elif gap > thresh and pre <= thresh and post > thresh:
            data[5] += 1
        elif gap > thresh and pre > thresh and post <= thresh:
            data[6] += 1
        elif gap > thresh and pre > thresh and post > thresh:
            data[7] += 1
            outliers.append(s)
    return data/N,sorted(outliers)

def run(fout=None):
    N2U = 'N2U'
    JSH = 'JSH'
    _remove = ['VC01','VD01','VB01','VB02']

    N2U = from_db(N2U,adjacency=True,chemical=True,
                  electrical=True,remove=_remove,dataType='networkx')
    JSH = from_db(JSH,adjacency=True,chemical=True,
                  electrical=True,remove=_remove,dataType='networkx')
    fig,ax = plt.subplots(1,1,figsize=(18,10))
    plot_specificity(ax,N2U,JSH,fout=fout)
    plt.show()

if __name__=='__main__':
    run()
