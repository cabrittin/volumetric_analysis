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

def run(fout=None):
    db = 'N2U'
    wbe = aux.read.into_list2('./data/'+db+'_wbe.csv')
    sbe = aux.read.into_list2('./data/'+db+'_sbe.csv')
    ie = aux.read.into_list2('./data/'+db+'_ie.csv')

    pre_mean = [np.mean([float(d[1]) for d in wbe if float(d[1]) >= 0]),
                np.mean([float(d[1]) for d in sbe if float(d[1]) >= 0]),
                np.mean([float(d[0]) for d in ie if float(d[0]) >=0 ])]
    gap_mean = [np.mean([float(d[2]) for d in wbe if float(d[2]) >= 0]),
                np.mean([float(d[2]) for d in sbe if float(d[2]) >= 0]),
                np.mean([float(d[1]) for d in ie if float(d[1]) >= 0])]

        
    pre_std = [np.std([float(d[1]) for d in wbe if float(d[1]) >= 0]),
               np.std([float(d[1]) for d in sbe if float(d[1]) >= 0]),
               np.std([float(d[0]) for d in ie if float(d[0]) >=0 ])]
    gap_std = [np.std([float(d[2]) for d in wbe if float(d[2]) >= 0]),
               np.std([float(d[2]) for d in sbe if float(d[2]) >= 0]),
               np.std([float(d[1]) for d in ie if float(d[1]) >= 0])]
        
    fig,ax = plt.subplots(1,1,figsize=(10,8))
    expplt.plot_cam_lus(ax,[pre_mean,gap_mean],
                        [pre_std,gap_std],
                        fout = fout)
    plt.show()    


if __name__=='__main__':
    run()
