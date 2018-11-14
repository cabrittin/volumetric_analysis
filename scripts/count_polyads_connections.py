"""
count_polyads_connection.py

Fraction synaptic contacts that have at least one polyadic synapse.

created: Christopher Brittin
date: 01 November 2018

"""
import sys
sys.path.append(r'./volumetric_analysis')
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

from connectome.load import from_db

SCREEN = ['old','duplicate','Frag','error','unk','sheath']
mpl.rcParams['xtick.labelsize'] = 32
mpl.rcParams['ytick.labelsize'] = 32

fout = '/home/cabrittin/Dropbox/PhD/sr_vol/figs2/polyad_frequency_syn.png'

ADULT_COL = '#FFC300'
L4_COL = '#3380FF'

def count_mon_poly(_db):
    _remove = ['VC01','VD01','VB01','VB02']
    C = from_db(_db,adjacency=True,chemical=True,add_poly=True,remove=_remove)
    C.reduce_to_adjacency()

    mon,poly = 0,0
    for e in C.C.es:
        if e['S'] > 0 and e['Sp'] == 0:
            mon += 1
        elif e['Sp'] > 0:
            poly += 1

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

    ax.set_ylabel('# synaptic connections',fontsize=32)
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

    
