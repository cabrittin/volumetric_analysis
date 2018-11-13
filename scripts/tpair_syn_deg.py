"""
tpair_syn_deg.py

Distribution of synaptic degree differences for 
homologous neurons

created: Christopher Brittin
date: 01 November 2018 

"""
import sys
sys.path.append(r'./volumetric_analysis')
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import ttest_rel
from scipy.stats import ttest_ind

#Brittin modules
from Connectome.load import from_db
from Networks.stats import *
from Figures.stats import *
import aux

mpl.rcParams['xtick.labelsize'] = 32
mpl.rcParams['ytick.labelsize'] = 32
lr_pairs = './mat/lr_neurons.txt'

def get_deg_data(C,lr):
    ld_out = get_deg(C.C,vertices=lr[0],mode='OUT')
    rd_out = get_deg(C.C,vertices=lr[1],mode='OUT')
    ld_in = get_deg(C.C,vertices=lr[0],mode='IN')
    rd_in = get_deg(C.C,vertices=lr[1],mode='IN')
    ld_gap = get_deg(C.E,vertices=lr[0],mode='ALL')
    rd_gap = get_deg(C.E,vertices=lr[1],mode='ALL')
    return [ld_gap,ld_out,ld_in],[rd_gap,rd_out,rd_in]

def get_dev_deg_data(C1,C2,cells):
    d1_out = get_deg(C1.C,vertices=cells,mode='OUT')
    d2_out = get_deg(C2.C,vertices=cells,mode='OUT')
    d1_in = get_deg(C1.C,vertices=cells,mode='IN')
    d2_in = get_deg(C2.C,vertices=cells,mode='IN')
    d1_gap = get_deg(C1.E,vertices=cells,mode='ALL')
    d2_gap = get_deg(C2.E,vertices=cells,mode='ALL')
    return [d1_gap,d1_out,d1_in],[d2_gap,d2_out,d2_in]    



def run(fout=None):
    N2U = 'N2U'
    JSH = 'JSH'
    _lr = aux.read.into_list2(lr_pairs)
    lr = list(zip(*_lr))
    _remove = ['VC01','VD01','VB01','VB02']
    n2u = from_db(N2U,adjacency=True,chemical=True,electrical=True,
                  remove=_remove)
    ndeg_l,ndeg_r = get_deg_data(n2u,lr)

    jsh = from_db(JSH,adjacency=True,chemical=True,electrical=True,
                  remove=_remove)
    jdeg_l,jdeg_r = get_deg_data(jsh,lr)

    cells = sorted((set(n2u.neurons)&set(jsh.neurons))-set(_remove))
    bdeg_l,bdeg_r = get_dev_deg_data(n2u,jsh,cells)

    data = [[ndeg_l[0][k1] - ndeg_r[0][k2] for (k1,k2) in _lr],
            [jdeg_l[0][k1] - jdeg_r[0][k2] for (k1,k2) in _lr],
            [bdeg_l[0][k] - bdeg_r[0][k] for k in cells],
            [ndeg_l[1][k1] - ndeg_r[1][k2] for (k1,k2) in _lr],
            [jdeg_l[1][k1] - jdeg_r[1][k2] for (k1,k2) in _lr],
            [bdeg_l[1][k] - bdeg_r[1][k] for k in cells],
            [ndeg_l[2][k1] - ndeg_r[2][k2] for (k1,k2) in _lr],
            [jdeg_l[2][k1] - jdeg_r[2][k2] for (k1,k2) in _lr],
            [bdeg_l[2][k] - bdeg_r[2][k] for k in cells]]

            
    #for i in [0,1,2]:
    #    print(ttest_rel(ndeg_l[i],ndeg_r[i]))
    #    print(ttest_rel(jdeg_l[i],jdeg_r[i]))
    #    print(ttest_rel(bdeg_l[i],bdeg_r[i]))

    tval0,pval0 = ttest_ind(data[0],data[2])
    tval1,pval1 = ttest_ind(data[0],data[1])
    tval1,pval2 = ttest_ind(data[1],data[2])
    tval1,pval3 = ttest_ind(data[3],data[4])
    tval1,pval4 = ttest_ind(data[4],data[5])
    tval1,pval5 = ttest_ind(data[6],data[7])
    tval1,pval6 = ttest_ind(data[7],data[8])
        

    pval = [
        (0,1,pval1),(1,2,pval2),
        (3,4,pval3),(4,5,pval4),
        (6,7,pval5),(7,8,pval6)]

    fig,ax = plt.subplots(1,1,figsize=(12,10))
    tpair_syn(ax,data,pval,
              fout=fout,
              ylabel='synaptic degree difference',
              title = 'Homologous synaptic degrees',
              xticklabels = ['$d_{\mathrm{gap}}$',
                             '$d_{\mathrm{pre}}$',
                             '$d_{\mathrm{post}}$'],
              ylim=[-20,20])
    plt.tight_layout()
    plt.show()

if __name__=='__main__':
    run()
