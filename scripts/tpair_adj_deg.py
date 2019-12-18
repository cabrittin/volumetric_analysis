"""
tpair_adj_deg.py

Plots distribution of the differences in homologous adjacency degrees.

created: Christopher Brittin
date: 01 November 2018

"""
import sys
sys.path.append(r'./volumetric_analysis')
import matplotlib.pyplot as plt
from scipy.stats import ttest_rel
from scipy.stats import ttest_ind
from scipy.stats import wilcoxon
from scipy.stats import ttest_1samp

#Brittin modules
from connectome.load import from_db
from networks.stats import *
from figures.stats import *
import aux


lr_pairs = './mat/lr_neurons.txt'

def write_out(nodes,score1,score2,_fout):
    _out = []
    for (n1,n2) in nodes:
        _out.append([n1,n2,score1[n1],score2[n2],abs(score1[n1]-score2[n2])])

    aux.write.from_list(_fout,_out)

def run(fout=None,source_data=None):
    
    N2U = 'N2U'
    JSH = 'JSH'

    _lr = [l for l in aux.read.into_list2(lr_pairs) if '#' not in l[0]]
    lr = list(zip(*_lr))
    nodes = lr[0]
    _remove = ['VC01','VD01','VB01','VB02']

    n2u = from_db(N2U,adjacency=True,remove=_remove)
    lnd = get_adj_deg(n2u,vertices = lr[0])
    rnd = get_adj_deg(n2u,vertices = lr[1])
    if source_data:
        fsplit = source_data.split('.')
        dout = fsplit[0] + '_adult_contralateral.' + fsplit[1]
        write_out(_lr,lnd,rnd,dout)

    jsh = from_db(JSH,adjacency=True,remove=_remove)
    ljd = get_adj_deg(jsh,vertices = lr[0])
    rjd = get_adj_deg(jsh,vertices = lr[1])
    if source_data:
        fsplit = source_data.split('.')
        dout = fsplit[0] + '_l4_contralateral.' + fsplit[1]
        write_out(_lr,ljd,rjd,dout)

    cells = sorted((set(n2u.neurons)&set(jsh.neurons))-set(_remove))
    bnd = get_adj_deg(n2u,vertices = cells)
    bjd = get_adj_deg(jsh,vertices = cells)
    if source_data:
        fsplit = source_data.split('.')
        dout = fsplit[0] + '_adult_l4_homologous.' + fsplit[1]
        write_out(zip(cells,cells),bnd,bjd,dout)    
    
    
    ndelta = [lnd[k1] - rnd[k2] for (k1,k2) in _lr]
    jdelta = [ljd[k1] - rjd[k2] for (k1,k2) in _lr]
    bdelta = [bnd[k] - bjd[k] for k in bnd.keys() if k in bjd.keys()]

    data = [ndelta,jdelta,bdelta]
    
    print(bdelta)
    print('Stats:')
    print_wilcoxon(ndelta,'Adult L/R')
    print_wilcoxon(jdelta,'L4 L/R')
    print_wilcoxon(bdelta,'Adult/L4',alternative="greater")
    
    tval1,pval1 = ttest_ind(ndelta,jdelta)
    tval2,pval2 = ttest_ind(jdelta,bdelta)
    tval3,pval3 = ttest_ind(ndelta,bdelta)
    
    pval = [(0,1,pval1),(1,2,pval3)]
    fig,ax = plt.subplots(1,1,figsize=(12,10))
    tpair_adj_deg(ax,data,pval,fout=fout)
    plt.show()


if __name__ == '__main__':
    run()
