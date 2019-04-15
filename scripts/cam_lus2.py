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

def run(fout=None):
    db = 'N2U'
    wbe = aux.read.into_list2('./cam_analysis/results/conserved_deg3_wbe.csv' )
    sbe = aux.read.into_list2('./cam_analysis/results/conserved_deg3_sbe.csv')
    ie = aux.read.into_list2('./cam_analysis/results/conserved_deg3_ie.csv')

    sbe = [s for s in sbe if s[0] not in SBE_SCREEN]
    wbe_pre = [float(d[1]) for d in wbe if float(d[1]) >= 0]
    sbe_pre = [float(d[1]) for d in sbe if float(d[1]) >= 0]
    ie_pre  = [float(d[0]) for d in ie if float(d[0]) >= 0]
    wbe_gap = [float(d[2]) for d in wbe if float(d[2]) >= 0]
    sbe_gap = [float(d[2]) for d in sbe if float(d[2]) >= 0]
    ie_gap = [float(d[1]) for d in ie if float(d[1]) >= 0]

    pre_size = [len(wbe_pre),len(sbe_pre),len(ie_pre)]
    gap_size = [len(wbe_gap),len(sbe_gap),len(ie_gap)]
    pre_mean = [np.mean(wbe_pre),np.mean(sbe_pre),np.mean(ie_pre)]
    gap_mean = [np.mean(wbe_gap),np.mean(sbe_gap),np.mean(ie_gap)]
    pre_std =  [np.std(wbe_pre),np.std(sbe_pre),np.std(ie_pre)]
    gap_std =  [np.std(wbe_gap),np.std(sbe_gap),np.std(ie_gap)]

    fig,ax = plt.subplots(1,1,figsize=(12,8))
    expplt.plot_cam_lus(ax,[pre_mean,gap_mean],
                        [pre_std,gap_std],
                        fout = fout)
    print(pre_size,gap_size)
    ax.set_xticklabels(('WBE\n($n=%d,%d$)'%(pre_size[0],gap_size[0]),
                        'SBE\n($n=%d,%d$)'%(pre_size[1],gap_size[1]),
                        'IE\n($n=%d,%d$)'%(pre_size[2],gap_size[2])))
    #ax.xaxis.set_tick_params(labelsize=28)
    plt.tight_layout()
    if fout: plt.savefig(fout)
    plt.show()    


if __name__=='__main__':
    run()
