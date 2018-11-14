"""
count_polyads_repeat.py

Fraction of polyadic synapses that are conserved between homologous cells.

created: Christopher Brittin
date: 01 November 2018

"""
import sys
sys.path.append(r'./volumetric_analysis')
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from itertools import combinations

import db
import aux

SCREEN = ['old','duplicate','Frag','error','unk','sheath']
mpl.rcParams['xtick.labelsize'] = 32
mpl.rcParams['ytick.labelsize'] = 32

ADULT_COL = '#FFC300'
L4_COL = '#3380FF'
AL_COL = '#14FA29'

lr_dict = './mat/lr_dict.txt'
_group = './mat/musgrp.txt'

def scrub_neurons(_neurons):
    neurons = []
    for n in _neurons:
        remove = False 
        for s in SCREEN:
            if s in n:
                remove = True
                break
        if not remove: neurons.append(n)
    return neurons

def get_poly_synapse(_db):
    if _db == 'N2U':
        end = 325
    else:
        end = 500    

    group = aux.read.into_map(_group)
    lrd = aux.read.into_dict(lr_dict)
    syn = []
    
    con = db.connect.default(_db)
    cur = con.cursor()
    
    neurons = scrub_neurons(db.mine.get_neurons(cur))
    synapses = db.mine.get_synapse_data(cur,'chemical',end=end)
    
    for s in synapses:
        pre = s[0]
        _post = [p for p in s[1].split(',') if p in neurons]
        post = []
        for _p in _post:
            if _p in group: _p = group[_p]
            if _p not in post: post.append(_p)
        post = sorted(post)
        N = len(post)
        if N > 1:
            comb = combinations(post,2)
            for (c1,c2) in comb:
                tmp = [pre,c1,c2]
                if tmp not in syn: syn.append([pre,c1,c2])

    return syn
                
def count_repeats(syn1,syn2):
    count = 0
    for s in syn2:
        if s in syn1:
            count += 1
            idx = syn1.index(s)
    return count


def run(fout=None):
    lrd = aux.read.into_lr_dict(lr_dict)
    n2u = get_poly_synapse('N2U')
    jsh = get_poly_synapse('JSH')

    rn2u = []
    for [s1,s2,s3] in n2u:
        try:
            tmp = [lrd[s1],lrd[s2],lrd[s3]]
            if tmp not in rn2u: rn2u.append(tmp)
        except:
            pass

    rjsh = []
    for [s1,s2,s3] in jsh:
        try:
            tmp = [lrd[s1],lrd[s2],lrd[s3]]
            if tmp not in rn2u: rjsh.append(tmp)        
        except:
            pass

    ncount = count_repeats(n2u,rn2u)
    jcount = count_repeats(jsh,rjsh)
    bcount = count_repeats(n2u,jsh)

    print(ncount,jcount,bcount)
    print(len(rn2u),len(rjsh),len(jsh))

    nfrac = ncount / float(len(rn2u))
    jfrac = jcount / float(len(rjsh))
    bfrac = bcount / float(len(jsh))

    frac = [nfrac,jfrac,bfrac]

    n_groups = 3
    bar_width = 0.25
    opacity = 0.5
    index = np.array([1,2,3])
    
    fig,ax = plt.subplots(1,1,figsize=(10,10))
    rects1 = ax.bar(index, frac, bar_width,
                    color=[ADULT_COL,L4_COL,AL_COL])
    
    ax.set_ylabel('Fraction of polyads',fontsize=32)
    ax.set_title('Polyads conserved among homologs',fontsize=32)
    ax.set_xticks(index)
    ax.set_xticklabels(('Adult L/R', 'L4 L/R', 'Adult/L4'))
    ax.set_xlim([0.5,3.5])
    ax.set_ylim([0,1])
    fig.tight_layout()
    if fout: plt.savefig(fout)
    plt.show()

    
if __name__=="__main__":
    run()
