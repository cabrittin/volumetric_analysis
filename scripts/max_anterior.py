"""
max_anterior.py

The projection distance of cells into the nerve ring. Posterior nuclei 
(blue) project anteriorly. Anterior nuclei (red) project posteriorly. 
Cells grouped based on function and length of projection (see main text). 
Vertical bars are the range of projections for given cell group. Width of 
bar is the fraction of cells at given length. Middle ticks are the mean 
and median.

created: Christopher Brittin
date: 01 November 2018 

"""

import sys
sys.path.append(r'./volumetric_analysis')
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches

mpl.rcParams['xtick.labelsize'] = 28
mpl.rcParams['ytick.labelsize'] = 28
mpl.rc('font', weight='bold')

#Brittin modules
import aux

MAXANT = 200
SCALE = 80/1000.
SCREEN = ['SABD']

n2u = './data/N2U_maxAnterior.csv'   
jsha = './data/JSH_maxAnterior.csv'
jshp = './data/JSH_maxPosterior.csv'

neuron_class = './mat/nerve_ring_classes.txt'
class_order_a = ['Sp1','Sp2','I1','I2','SMN','HMNp']
class_order_p = ['Sa','HMNa']

def maxAnt(pos):
    return (MAXANT - int(pos))*SCALE

def maxPost(pos):
    return int(pos)*SCALE

def format_data(fin,func,class_order,nclass,screen=[]):
    pos = aux.read.into_dict(fin)
    data = dict([(c,[]) for c in nclass.values()])
    for n in pos:
        if n in screen: continue
        data[nclass[n]].append(func(pos[n]))
    return [data[c] for c in class_order]

def print_data(fin,class_order,nclass):
    pos = aux.read.into_dict(fin)
    rev_nclass = dict([(c,[]) for c in set(nclass.values())])
    for n in nclass:
        rev_nclass[nclass[n]].append(n)
        
    for c in class_order:
        print('Class: ' + c)
        for n in sorted(rev_nclass[c]):
            if n in pos:
                print('\t* %s -- %s' %(n,pos[n]))


def run(fout=None):
    nclass = aux.read.into_dict(neuron_class)

    data_a = format_data(jsha,maxAnt,class_order_a,nclass,SCREEN)
    print_data(jsha,class_order_a,nclass)

    data_p = format_data(jshp,maxPost,class_order_p,nclass,SCREEN)
    print_data(jshp,class_order_p,nclass)

    N = len(class_order_a)
    pos_a = [i for i in range(N)]
    pos_p = [N, N+1]


    fig = plt.figure(figsize=(15,10))
    ax1 = fig.add_subplot(111)
    v1 = ax1.violinplot(data_a,pos_a,points=10,widths=0.5,
                        showmeans=True, showextrema=True, showmedians=True)
    for pc in v1['bodies']:
        pc.set_color('#6CC7FF')
        pc.set_alpha(0.3)


    ax2 = fig.add_subplot(111,sharex=ax1, frameon=False)
    v2 = ax2.violinplot(data_p,pos_p,points=10,widths=0.5,
                        showmeans=True, showextrema=True, showmedians=True)
    for pc in v2['bodies']:
        pc.set_color('#D43F3A')
        pc.set_alpha(0.3)

    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    ax2.set_ylabel("Posterior projection, $\mu$m",
                   fontsize=32,rotation=-90,
                   labelpad=40,weight='bold')
    ax2.set_ylim([0,16])
    ax2.invert_yaxis()

    ax1.set_xticks(pos_a + pos_p)
    ax1.set_xticklabels(class_order_a + class_order_p)
    ax1.set_ylabel('Anterior projection, $\mu$m',
                   fontsize=32,weight='bold')
    ax1.set_title('Projection distance into nerve ring'
                  ,fontsize=32,weight='bold')
    ax1.set_ylim([0,16])
    ax1.yaxis.grid(True)

    post = mpatches.Patch(color='#6CC7FF', alpha=0.3, label='Posterior cells')
    ant = mpatches.Patch(color='#D43F3A', alpha=0.3, label='Anterior cells')
    
    
    plt.legend(handles=[post,ant],fontsize=18,loc='best')
    if fout: plt.savefig(fout)
    plt.show()

if __name__=="__main__":
    run()

    
