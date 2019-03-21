"""
dist_sa_subgrp2.py

Plots the distributions of surface areas

@author Christopher Brittin
@date 2019-03-17
"""

import sys
sys.path.append('./volumetric_analysis')
import argparse

from connectome.load import from_db
import aux
from figures.stats import *


SCALE = 450e-6

def get_sa(stats,nclass):
    sa = {}
    i = 0
    for d in stats:
        if d[0] not in nclass: continue
        nc = nclass[d[0]]
        if nc not in sa: sa[nc] = []
        sa[nc].append(int(d[1])*SCALE)
    return sa  

if __name__== '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('n2u',
                        action="store",
                        help="Path to adult size stats file")
   
    parser.add_argument('jsh',
                        action="store",
                        help="Path to l4 size stats file")
   

    parser.add_argument('nclass',
                        action="store",
                        help="Path the class file")

    parser.add_argument('-o','--output',
                        dest = 'fout',
                        action="store",
                        required= False,
                        default = None,
                        help="Output file, if you wish to save the output."
                        )

    params = parser.parse_args()
    
    nclass = aux.read.into_dict(params.nclass)
 
    n2u = get_sa(aux.read.into_list2(params.n2u),nclass)
    jsh = get_sa(aux.read.into_list2(params.jsh),nclass)

  
    NSA = [n2u['Sp1'] + n2u['Sp2'],
            n2u['Sa'],n2u['I1'],n2u['I2'],n2u['SMN'],
            n2u['HMNp'],n2u['HMNa']]
    
    JSA = [jsh['Sp1'] + jsh['Sp2'],
            jsh['Sa'],jsh['I1'],jsh['I2'],jsh['SMN'],
            jsh['HMNp'],jsh['HMNa']]

    data = []
    for i in range(len(NSA)):
        data.append(NSA[i])
        data.append(JSA[i])
    
    fig,ax = plt.subplots(1,1,figsize=(15,10))
    dist_adj_subgroups2(ax,data)
    ax.set_ylim([0,500])    
    ax.set_ylabel(r'Surface area ($\mu$m$^2$)',fontsize=32)
    ax.set_title('Cell surface area by group',fontsize=36)
    if params.fout: plt.savefig(params.fout)
    plt.show()




