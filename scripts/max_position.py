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

mpl.rcParams['xtick.labelsize'] = 24
mpl.rcParams['ytick.labelsize'] = 28
mpl.rc('font', weight='bold')

#Brittin modules
import aux

SCREEN = ['SABD']

N2UPOSFILE = './data/N2U_maxposition.csv'
JSHPOSFILE = './data/JSH_maxposition.csv'

N2UCONTINSCREEN = './mat/n2u_sp_contins.txt'

MAXPOS = 200
SCALE = 0.08

class_order_a = ['Sp1','Sp2','I1','I2','SMN','HMNp']
class_order_p = ['Sa','HMNa']

def format_data(fin,ant_class_order,post_class_order,screen=[],contin_screen=None):
    pos = aux.read.into_list2(fin)
    #First pass organize cells
    data1 = {}
    for (cell,nclass,contin,apos,ppos) in pos:
        if cell not in data1:
            data1[cell] = {
                'nclass':nclass,
                'apos' : [],
                'ppos' : []
            }
        if contin_screen and cell in contin_screen:
            if int(contin) != int(contin_screen[cell]): continue
        _apos,_ppos = int(apos),int(ppos)
        if _apos > -1: data1[cell]['apos'].append(_apos)
        if _ppos > -1: data1[cell]['ppos'].append(_ppos)

    data2 = {}
    for cell in data1:
        if not data1[cell]['apos'] or not data1[cell]['ppos']:
            continue
        data2[cell] = {
            'nclass':data1[cell]['nclass'],
            'apos' : min(data1[cell]['apos']),
            'ppos' : max(data1[cell]['ppos'])
            }
            
    data = {}
    for cell in sorted(data2):
        if cell in screen: continue
        nclass = data2[cell]['nclass']
        _apos = data2[cell]['apos']
        _ppos = data2[cell]['ppos']
        if nclass not in data: data[nclass] = []
        if nclass in ant_class_order:
            apos = float(_apos)
            if apos > -1:
                apos = (MAXPOS - apos)*SCALE
                data[nclass].append(apos)
                if nclass in ['Sp1','Sp2']:
                    print(cell,nclass,apos)
        if nclass in post_class_order:
            ppos = float(_ppos)
            ppos = float(min(MAXPOS,ppos))*SCALE
            if ppos > -1: data[nclass].append(ppos)
    adata = [data[c] for c in ant_class_order]
    pdata = [data[c] for c in post_class_order]
    return adata,pdata

def run(_db='N2U',fout=None):
    if _db == 'JSH':
        fin = JSHPOSFILE
        contin_screen = None
    elif _db == 'N2U':
        fin = N2UPOSFILE
        contin_screen = aux.read.into_dict(N2UCONTINSCREEN)
    
    data_a,data_p = format_data(fin,class_order_a,class_order_p,
                                screen=SCREEN,contin_screen=contin_screen)


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

    ticklabels = []
    for i in range(len(class_order_a)):
        ticklabels.append(class_order_a[i] + '\n(n=%d)'%len(data_a[i]))
    for i in range(len(class_order_p)):
        ticklabels.append(class_order_p[i] + '\n(n=%d)'%len(data_p[i]))    
    ax1.set_xticks(pos_a + pos_p)
    ax1.set_xticklabels(ticklabels)
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

    
