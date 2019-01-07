"""
pre_post_specificity.py

Plot of p_s^pre vs. p_s^post in the adult. Cells fall into one of 
four categories: (+,+), (+,-), (+,-) or (-,-) indicated by red, yellow, 
orange and blue, respectively. Table gives the fraction of neurons in each 
category. Outlier neurons in the last category (-,-) are labeled. Homologous 
neurons are considered outliers if both p_s^pre and p_s^post are greater than α = 0.05. Red dashed line marks where the probability is 0.05.

Author: Christopher Brittin
Created: 07 February 2018

"""

import sys
sys.path.append(r'./volumetric_analysis/')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl


from connectome.load import from_db
import connectome.synspecificity as synspec
import aux

mpl.rcParams['xtick.labelsize'] = 32
mpl.rcParams['ytick.labelsize'] = 32

lr_dict = './mat/lr_dict.txt'
homologs = './mat/homologs.txt'
left_nodes = './mat/left_nodes.txt'
db = 'N2U'
_remove = ['VC01','VD01','VB01','VB02']

def get_outliers(pre,post,ndict):
    prob_both = []
    outliers = []
    count1,count2,count3,count4 = 0,0,0,0
    for n in sorted(pre):
        prob_both.append((pre[n].p,post[n].p))
        if (pre[n].p > 0.05 and post[n].p > 0.05):
            count4 += 1
            print(n,
                  pre[n].p,post[n].p,
                  pre[n].c1,post[n].c1,
                  pre[n].k,post[n].k,
                  pre[n].M,post[n].M)
            outliers.append([ndict[n],pre[n].p,post[n].p])
        elif pre[n].p <= 0.05 and post[n].p <= 0.05:
            count1 += 1
        elif pre[n].p <= 0.05 and post[n].p > 0.05:
            count2 += 1
        elif pre[n].p > 0.05 and post[n].p <= 0.05:
            count3 += 1
    N = float(len(pre))
    print('Type1: %1.4f, Type2: %1.4f, Type3: %1.4f, Type4: %1.4f'
          %(count1/N,count2/N,count3/N,count4/N))
    print('Type1: %d, Type2: %d, Type3: %d, Type4: %d, Total: %d'
          %(count1,count2,count3,count4,N))
    return prob_both,outliers,[count1/N,count2/N,count3/N,count4/N]

def plot_pre_post_specificity(prob_both,outliers,probs,fout=None):
    plt.figure(figsize=(18,10))
    x,y = zip(*prob_both)
    plt.plot(x,y,'ko')
    plt.axhspan(ymin=0,ymax=0.05,xmin=0.0,xmax=0.05,facecolor='#FA4E29',alpha=0.5)
    plt.axhspan(ymin=0,ymax=0.05,xmin=0.05,xmax=1.0,facecolor='#F7CF64',alpha=0.5)
    plt.axhspan(ymin=0.05,ymax=1.0,xmin=0.0,xmax=0.05,facecolor='#F7F364',alpha=0.5)
    plt.axhspan(ymin=0.05,ymax=1.0,xmin=0.05,xmax=0.7,facecolor='#648EF7',alpha=0.5,linewidth=0)
    plt.axhspan(ymin=0.05,ymax=1.0,xmin=0.9,xmax=1.0,facecolor='#648EF7',alpha=0.5,linewidth=0)
    plt.axhspan(ymin=0.05,ymax=0.5,xmin=0.7,xmax=0.9,facecolor='#648EF7',alpha=0.5,linewidth=0)
    plt.axhspan(ymin=0.7,ymax=1.0,xmin=0.7,xmax=0.9,facecolor='#648EF7',alpha=0.5,linewidth=0)
    plt.text(0.74,0.72,'+',fontsize=32)
    plt.text(0.84,0.72,'-',fontsize=32)
    plt.text(0.67,0.64,'+',fontsize=32)
    plt.text(0.67,0.54,'-',fontsize=32)
    plt.text(0.78,0.78,'$p^{\mathrm{pre}}_s$',fontsize=32)
    plt.text(0.63,0.63,'$p^{\mathrm{post}}_s$',fontsize=32,rotation=90)
    plt.text(0.72,0.63,'%1.2f'%probs[0],fontsize=32)
    plt.text(0.72,0.53,'%1.2f'%probs[1],fontsize=32)
    plt.text(0.82,0.63,'%1.2f'%probs[2],fontsize=32)
    plt.text(0.82,0.53,'%1.2f'%probs[3],fontsize=32)
    plt.axvline(x=0.05,ymin=0.05,ymax=1.0,linewidth=4,linestyle='--',color='r')
    plt.axhline(y=0.05,xmin=0.05,xmax=1.0,linewidth=4,linestyle='--',color='r')
    rect2 = patches.Rectangle((0.7,0.5),0.1,0.1,linewidth=2,
                              edgecolor='k',facecolor='#F7CF64',alpha=0.5)
    rect1 = patches.Rectangle((0.7,0.6),0.1,0.1,linewidth=4,
                              edgecolor='k',facecolor='#FA4E29',alpha=0.5)
    rect3 = patches.Rectangle((0.8,0.6),0.1,0.1,linewidth=4,
                              edgecolor='k',facecolor='#F7F364',alpha=0.5)
    rect4 = patches.Rectangle((0.8,0.5),0.1,0.1,linewidth=4,
                              edgecolor='k',facecolor='#648EF7',alpha=0.5)
    plt.text(0.63,0.43,'Fraction of neurons',fontsize=32)
    plt.gca().add_patch(rect1)
    plt.gca().add_patch(rect2)
    plt.gca().add_patch(rect3)
    plt.gca().add_patch(rect4)
    
    plt.xlabel('$p^{\mathrm{pre}}_s$',fontsize=32)
    plt.ylabel('$p^{\mathrm{post}}_s$',fontsize=32)
    plt.xticks([0,0.2,0.4,0.6,0.8,1.0],['',0.2,0.4,0.6,0.8,1.0])
    for o in outliers:
        _x = o[1]+0.01
        _y = o[2]+0.01
        if o[1] == 1.0: _x = o[1]-0.05
        if o[2] == 1.0: _y = o[2]-0.05
        plt.text(_x,_y,o[0],fontsize=24)
    plt.xlim([0,1])
    plt.ylim([0,1])
    plt.tight_layout()
    if fout: plt.savefig(fout)
    
def run(fout=None):
    C = from_db(db,adjacency=True,chemical=True,electrical=True,
                remove=_remove,dataType='networkx')

    #C.remove_self_loops()
    C.reduce_to_adjacency()    

    lrd = aux.read.into_lr_dict(lr_dict)

    nclass = aux.read.into_list2(homologs)
    ndict = {}    
    
    for n in nclass:
        for _n in n[1:]:ndict[_n] = n[0]
    
    left = aux.read.into_list(left_nodes)
    
    left.remove('CEHDL')
    left.remove('CEHVL')
    left.remove('HSNL')
    left.remove('PVNL')
    left.remove('PLNL')

    _pre = synspec.bilateral_specificity(C,left,lrd)
    _post = synspec.bilateral_specificity(C,left,lrd,mode='post')

    prob_both,outliers,probs = get_outliers(_pre,_post,ndict)
    plot_pre_post_specificity(prob_both,outliers,probs,fout=fout)
    plt.show()


if __name__=='__main__':
    run()
