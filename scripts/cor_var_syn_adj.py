"""
cor_var_syn_adj.py

Correlation between synaptic and adjacency degree differences 
betweens homologous neurons.

created: Christopher Brittin
date: 01 November 2018

"""
import sys
sys.path.append(r'./volumetric_analysis')
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import ttest_rel
from scipy.stats import ttest_ind
from scipy.stats import pearsonr
from scipy.stats import spearmanr

#Brittin modules
from connectome.load import from_db
from networks.stats import *
from figures.stats import *
import aux

mpl.rcParams['xtick.labelsize'] = 32
mpl.rcParams['ytick.labelsize'] = 32

def get_deg_data(C,lr):
    L = [None]*4
    R = [None]*4

    L[0] = get_deg(C.A,vertices=lr[0],mode='ALL')
    L[1] = get_deg(C.E,vertices=lr[0],mode='ALL')
    L[2] = get_deg(C.C,vertices=lr[0],mode='OUT')
    L[3] = get_deg(C.C,vertices=lr[0],mode='IN')

    R[0] = get_deg(C.A,vertices=lr[1],mode='ALL')
    R[1] = get_deg(C.E,vertices=lr[1],mode='ALL')
    R[2] = get_deg(C.C,vertices=lr[1],mode='OUT')
    R[3] = get_deg(C.C,vertices=lr[1],mode='IN')    

    return L,R

def get_dev_deg_data(C1,C2,cells):
    D1 = [None]*4
    D2 = [None]*4

    D1[0] = get_deg(C1.A,vertices=cells,mode='ALL')
    D1[1] = get_deg(C1.E,vertices=cells,mode='ALL')
    D1[2] = get_deg(C1.C,vertices=cells,mode='OUT')
    D1[3] = get_deg(C1.C,vertices=cells,mode='IN')

    D2[0] = get_deg(C2.A,vertices=cells,mode='ALL')
    D2[1] = get_deg(C2.E,vertices=cells,mode='ALL')
    D2[2] = get_deg(C2.C,vertices=cells,mode='OUT')
    D2[3] = get_deg(C2.C,vertices=cells,mode='IN')    
    
    return D1,D2

ADULT_COL = '#FFC300'
L4_COL = '#3380FF'
AL_COL = '#14FA29'
lr_pairs = './mat/lr_neurons.txt'

def run(fout=None):
    N2U = 'N2U'
    JSH = 'JSH'
    _lr = aux.read.into_list2(lr_pairs)
    lr = list(zip(*_lr))
    _remove = ['VC01','VD01','VB01','VB02']
    n2u = from_db(N2U,adjacency=True,chemical=True,electrical=True,remove=_remove)
    ndeg_l,ndeg_r = get_deg_data(n2u,lr)

    jsh = from_db(JSH,adjacency=True,chemical=True,electrical=True,remove=_remove)
    jdeg_l,jdeg_r = get_deg_data(jsh,lr)

    cells = sorted((set(n2u.neurons)&set(jsh.neurons))-set(_remove))
    bdeg_1,bdeg_2 = get_dev_deg_data(n2u,jsh,cells)
    
    ndeg,jdeg,bdeg = [None]*4,[None]*4,[None]*4
    for i in range(4):
        ndeg[i] = [ndeg_l[i][k1] - ndeg_r[i][k2] for (k1,k2) in _lr]
        jdeg[i] = [jdeg_l[i][k1] - jdeg_r[i][k2] for (k1,k2) in _lr]
        bdeg[i] = [bdeg_1[i][k] - bdeg_2[i][k] for k in cells]

    ndeg = np.array(ndeg)
    jdeg = np.array(jdeg)
    bdeg = np.array(bdeg)
    data = []
    for A in [ndeg,jdeg,bdeg]:
        tmp = []
        for j in [1,2,3]:
            rho,p = spearmanr(A[0,:],A[j,:]) 
            tmp.append(rho**2)
        data.append(tmp)
        
    ind = np.arange(3)
    width = 0.15

    fig,ax = plt.subplots(1,1,figsize=(10,10))
    ax.axvspan(-0.5,0.5,facecolor='#C3C3C3')
    ax.axvspan(0.5,1.5,facecolor='#D8D7D7')
    ax.axvspan(1.5,2.5,facecolor='#C3C3C3')
    ax.grid(axis='y',color='0.9',linestyle='-',linewidth=1)
    rects0 = ax.bar(ind-width, data[0], width, color=ADULT_COL)
    rects1 = ax.bar(ind, data[1], width, color=L4_COL)
    rects2 = ax.bar(ind+width, data[2], width, color=AL_COL) 
    ax.set_xticks(ind)
    ax.set_xticklabels(('$r^2_{\mathrm{gap}}$',
                        '$r^2_{\mathrm{pre}}$',
                        '$r^2_{\mathrm{post}}$'))
    ax.set_yticks([0,0.05,0.1,0.15,0.2])
    ax.set_ylim([0,0.2])
    ax.set_xlim([-.5,2.5])
    ax.legend((rects0[0], rects1[0], rects2[0]),
              ('Adult L/R', 'L4 L/R', 'Adult/L4'),
              fontsize=18)
    ax.set_ylabel('Coefficient of determination',fontsize=32)
    ax.set_title('Correlation between synaptic and\nadjacency degree differences',
                 fontsize=32)
    plt.tight_layout()
    if fout:plt.savefig(fout)
    plt.show()

if __name__=='__main__':
    run()


