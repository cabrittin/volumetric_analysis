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

#Brittin modules
from connectome.load import from_db
from networks.stats import *
from figures.stats import *
import aux


lr_pairs = './mat/lr_neurons.txt'
dout1 = '../outputs/n2u_adj_def_diff.csv'
dout2 = '../outputs/jsh_adj_def_diff.csv'
dout3 = '../outputs/both_adj_def_diff.csv'

def write_out(node,score1,score2,_fout):
    _out = []
    for i in range(len(node)):
        _out.append([i,node[i],score1[i],score2[i],abs(score1[i]-score2[i])])

    aux.write.from_list(_fout,_out)

def run(fout=None):
    
    N2U = 'N2U'
    JSH = 'JSH'

    _lr = aux.read.into_list2(lr_pairs) 
    lr = list(zip(*_lr))
    nodes = lr[0]
    _remove = ['VC01','VD01','VB01','VB02']

    n2u = from_db(N2U,adjacency=True,remove=_remove)
    lnd = get_adj_deg(n2u,vertices = lr[0])
    rnd = get_adj_deg(n2u,vertices = lr[1])
    #write_out(nodes,lnd,rnd,dout1)

    jsh = from_db(JSH,adjacency=True,remove=_remove)
    ljd = get_adj_deg(jsh,vertices = lr[0])
    rjd = get_adj_deg(jsh,vertices = lr[1])
    #write_out(nodes,ljd,rjd,dout2)

    cells = sorted((set(n2u.neurons)&set(jsh.neurons))-set(_remove))
    bnd = get_adj_deg(n2u,vertices = cells)
    bjd = get_adj_deg(jsh,vertices = cells)
    #write_out(nodes,bnd,bjd,dout3)

    ndelta = [lnd[k1] - rnd[k2] for (k1,k2) in _lr]
    jdelta = [ljd[k1] - rjd[k2] for (k1,k2) in _lr]
    bdelta = [bnd[k] - bjd[k] for k in bnd.keys() if k in bjd.keys()]

    data = [ndelta,jdelta,bdelta]
    #print(np.mean(ndelta),np.std(ndelta))
    #print(np.mean(jdelta),np.std(jdelta))
    #print(np.mean(bdelta),np.std(bdelta))

    tval1,pval1 = ttest_ind(ndelta,jdelta)
    tval2,pval2 = ttest_ind(jdelta,bdelta)
    tval3,pval3 = ttest_ind(ndelta,bdelta)
    
    pval = [(0,1,pval1),(1,2,pval3)]
    fig,ax = plt.subplots(1,1,figsize=(12,10))
    tpair_adj_deg(ax,data,pval,fout=fout)
    plt.show()


if __name__ == '__main__':
    run()
