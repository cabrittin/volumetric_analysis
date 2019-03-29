"""
tpair_sa.py

Plots distribution of the differences in surface areas.

created: Christopher Brittin
date: 01 November 2018

"""
import sys
sys.path.append(r'./volumetric_analysis')
import argparse
import matplotlib.pyplot as plt
from scipy.stats import ttest_rel
from scipy.stats import ttest_ind
from scipy.stats import wilcoxon
from scipy.stats import mannwhitneyu
from scipy.stats import ttest_1samp

#Brittin modules
from connectome.load import from_db
from networks.stats import *
from figures.stats import *
import aux

SCALE = 450e-6

def write_out(nodes,score1,score2,_fout):
    _out = []
    for (n1,n2) in nodes:
        _out.append([n1,n2,score1[n1],score2[n2],abs(score1[n1]-score2[n2])])

    aux.write.from_list(_fout,_out)

if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('n2u',
                        action="store",
                        help="Path to adult size stats file")
   
    parser.add_argument('jsh',
                        action="store",
                        help="Path to l4 size stats file")
   

    parser.add_argument('lr_pairs',
                        action="store",
                        help="Path to L/R pairs file")

    parser.add_argument('-o','--output',
                        dest = 'fout',
                        action="store",
                        required= False,
                        default = None,
                        help="Output file, if you wish to save the output."
                        )

    params = parser.parse_args()
    

    _lr = aux.read.into_list2(params.lr_pairs) 
    lr = list(zip(*_lr))
    nodes = lr[0]
    _remove = ['VC01','VD01','VB01','VB02','Pharynx']

    n2u = aux.read.into_list2(params.n2u)
    lnd = dict([(d[0],int(d[1]) * SCALE) for d in n2u if d[0] in lr[0]])
    rnd = dict([(d[0],int(d[1]) * SCALE) for d in n2u if d[0] in lr[1]])
    
    jsh = aux.read.into_list2(params.jsh)
    ljd = dict([(d[0],int(d[1]) * SCALE) for d in jsh if d[0] in lr[0]])
    rjd = dict([(d[0],int(d[1]) * SCALE) for d in jsh if d[0] in lr[1]])
    
    ncells = set([d[0] for d in n2u]) -set(_remove)
    jcells = set([d[0] for d in jsh]) -set(_remove)
    cells = ncells & jcells
    bnd = dict([(d[0],int(d[1]) * SCALE) for d in n2u if d[0] in cells])
    bjd = dict([(d[0],int(d[1]) * SCALE) for d in jsh if d[0] in cells])


    ndelta = [2 * (lnd[k1] - rnd[k2]) / (lnd[k1] + rnd[k2])  for (k1,k2) in _lr]
    jdelta = [2 * (ljd[k1] - rjd[k2]) / (lnd[k1] + rnd[k2]) for (k1,k2) in _lr]
    bdelta = [2 * (bnd[k] - bjd[k]) / (bnd[k] + bjd[k]) for k in bnd.keys() if k in bjd.keys()]

    data = [ndelta,jdelta,bdelta]
    
    print('Stats:')
    print_wilcoxon(data[0],'Adult L/R')
    print_wilcoxon(data[1],'L4 L/R')
    print_wilcoxon(data[2],'Adult/L4')
    
    #tval1,pval1 = ttest_ind(data[0],data[1])
    #tval2,pval2 = ttest_ind(data[1],data[2])
    #tval3,pval3 = ttest_ind(data[0],data[2])
    
    tval1,pval1 = mannwhitneyu(data[0],data[1])
    tval2,pval2 = mannwhitneyu(data[1],data[2])
    tval3,pval3 = mannwhitneyu(data[0],data[2],alternative='less')

    pval = [(0,1,pval1),(1,2,pval3)]
    fig,ax = plt.subplots(1,1,figsize=(12,10))
    tpair_adj_deg(ax,data,pval)
    ax.set_ylim([-1,1])
    ax.set_ylabel('Normalized surface area difference',fontsize=32)
    ax.set_title('Homologous surface area difference',fontsize=38)
    
    ax.text(1.4,0.42,'n.s.',fontsize=18)
    ax.text(2.4,0.77,'****',fontsize=18)
    if params.fout: plt.savefig(params.fout)
    plt.show()

