"""
count_polyads.py

Fraction of synapses that are polyadic.

created: Christopher Brittin
date: 01 November 2018

"""
import sys
sys.path.append(r'./volumetric_analysis')
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

import db

SCREEN = ['old','duplicate','Frag','error','unk','sheath']
mpl.rcParams['xtick.labelsize'] = 32
mpl.rcParams['ytick.labelsize'] = 32

ADULT_COL = '#FFC300'
L4_COL = '#3380FF'

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

def count_mon_poly(_db):
    if db == 'N2U':
        end = 325
    else:
        end = 500    

    con = db.connect.default(_db)
    cur = con.cursor()
        
    neurons = scrub_neurons(db.mine.get_neurons(cur))
    synapses = db.mine.get_synapse_data(cur,'chemical',end=end)

    mon,poly = 0,0
    for s in synapses:
        pre = s[0]
        post = [p for p in s[1].split(',') if p in neurons]
        N = len(post)
        if N == 1:
            mon += 1
        elif N > 1:
            poly += 1

    con.close()
    return [mon,poly]

def run(fout=None):
    n2u = count_mon_poly('N2U')
    jsh = count_mon_poly('JSH')
    print(n2u,jsh)
    
    mon = [n2u[0],jsh[0]]
    poly = [n2u[1],jsh[1]]

    n_groups = 2
    bar_width = 0.25
    opacity = 0.5
    index = np.array([1,2])

    fig,ax = plt.subplots(1,1,figsize=(7,10))
    rects1 = ax.bar(index, mon, bar_width,color='r',
                    alpha=opacity,label='Monadic')

    rects2 = ax.bar(index+bar_width, poly, bar_width,color='b',alpha=opacity,
                    label='Polyadic')
    
    ax.set_ylabel('# synapses',fontsize=32)
    ax.set_title('Polyad frequency',fontsize=32)
    ax.set_xticks(index + bar_width / 2)
    ax.set_xticklabels(('Adult', 'L4'))
    ax.legend(fontsize=24)
    ax.set_xlim([0.75,2.5])
    fig.tight_layout()
    if fout: plt.savefig(fout)
    plt.show()

if __name__=="__main__":
    run()


    
