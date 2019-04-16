"""
cam_lus.py

The LUS scores for the combinatoric CAM models WBE, SBE and IE. 
Bar heights are the average LUS across neurons. Error bars represent 
the standard error.

created: Christopher Brittin
date: 01 November 2018
"""

import sys
sys.path.append(r'./volumetric_analysis')
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

import cam.cam_plots as expplt
import aux

"""
SBE_SCREEN : Removing these cells from the SBE results
because they are not in the WBE results. This is due to
a discrepancy between the adjacency data and the synaptic
connectivity data. The synaptic connectivity data says that
PLN makes synapses with SAA and SMB, but in the adjacency data
these cells do not make contact. There are 70/2319 (3%) synaptic
edges in the adult synaptic connectivity data that are not accounted
for in the adjacency data. There are 20/1571 (1.2%) synaptic 
edges in the l4 synaptic connectivity data that are not accounted for
in the adjacency data.
"""
SBE_SCREEN = ['PLNL','PLNR']


WBE = './cam_analysis/results/conserved_deg%d_wbe.csv'
SBE = './cam_analysis/results/conserved_deg%d_sbe.csv'
IE =  './cam_analysis/results/conserved_deg%d_ie.csv'

DEG = [1,2,3,4]

FOUT = '/home/cabrittin/Dropbox/PhD/sr_vol/figs2/fig9/cam_lus.png'


def plot_lus(ax,mu,yerr):
    error_kw = {'capsize': 10, 'capthick': 2,
                'ecolor': 'black','elinewidth':2}
    
    ax.axvspan(-0.5,0.4,facecolor='#C3C3C3')
    ax.axvspan(0.4,1.4,facecolor='#D8D7D7')
    ax.axvspan(1.4,2.4,facecolor='#C3C3C3')
    ax.axvspan(2.4,3.4,facecolor='#D8D7D7')
    ax.axvspan(3.4,4.4,facecolor='#C3C3C3')
    ax.axvspan(4.4,5.4,facecolor='#D8D7D7')

    (n,m) = mu.shape
    ind = np.arange(n)
    width = 0.2
   
    x = [ind-2*width,ind-width,ind,ind+width]
    rects = []
    col = ['r','g','b','y']
    for i in range(m):
        _ax = ax.bar(x[i],mu[:,i],width,color=col[i],
                    yerr=yerr[:,i],error_kw=error_kw)
        rects.append([_ax])
    
    ax.set_ylabel('LUS',fontsize=28)
    ax.set_xticks(ind-(0.5*width))
    ax.set_xticklabels(('WBE chem.','WBE gap','SBE chem.','SBE gap','IE chem.','IE gap'))
    ax.set_ylim([0,1.])
    ax.set_xlim([-0.5,n-0.5])
    ax.legend((rects[0][0],rects[1][0],rects[2][0],rects[3][0]),
            ('Conserved 1','Conserved 2','Conserved 3','Conserved 4'),
            fontsize=24)
    ax.set_title('Combinatorial CAM models', fontsize=28)
    


def run(fout=None):
   
    n = len(DEG)
    size = np.zeros((6,n))
    mu = np.zeros((6,n))
    std = np.zeros((6,n))
    
    idx = 0
    for deg in DEG:
        wbe = aux.read.into_list2(WBE%deg )
        sbe = aux.read.into_list2(SBE%deg)
        ie = aux.read.into_list2(IE%deg)

        sbe = [s for s in sbe if s[0] not in SBE_SCREEN]
        wbe_pre = [float(d[1]) for d in wbe if float(d[1]) >= 0]
        sbe_pre = [float(d[1]) for d in sbe if float(d[1]) >= 0]
        ie_pre  = [float(d[0]) for d in ie if float(d[0]) >= 0]
        wbe_gap = [float(d[2]) for d in wbe if float(d[2]) >= 0]
        sbe_gap = [float(d[2]) for d in sbe if float(d[2]) >= 0]
        ie_gap = [float(d[1]) for d in ie if float(d[1]) >= 0]
        
        _data = [wbe_pre,wbe_gap,sbe_pre,sbe_gap,ie_pre,ie_gap]
        size[:,idx] = np.array([len(d) for d in _data])
        mu[:,idx] = np.array([np.mean(d) for d in _data])
        std[:,idx] = np.array([np.std(d) for d in _data])
        
        
        idx += 1

    print(mu)

    fig,ax = plt.subplots(1,1,figsize=(16,8))
    plot_lus(ax,mu,std)
    #plt.tight_layout()
    plt.savefig(FOUT)
    plt.show()


if __name__=='__main__':
    run()
