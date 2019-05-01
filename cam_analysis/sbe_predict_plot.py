"""
sbe_predict_plot.py

Plots sbe predict results

@author Christopher Brittin
@data 23 April 2019
"""

import sys
sys.path.append('./volumetric_analysis')
import argparse
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

import aux

mpl.rcParams['xtick.labelsize'] = 24
mpl.rcParams['ytick.labelsize'] = 28


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('fin',
                        action='store',
                        help = 'Path to file')

    parser.add_argument('-o','--output',
                        action ='store',
                        dest = 'fout',
                        required = False,
                        default = None,
                        help = "Output file")

    params = parser.parse_args()


    data = aux.read.into_dict2(params.fin)

    pred = []
    for (k,d) in data.items():
        if float(d[0]) < 0: continue
        pred.append([float(d[2]),float(d[0]),float(d[1])])

    pred = np.array(pred)
    
    error_kw = {'capsize': 10, 'capthick': 2,
                'ecolor': 'black','elinewidth':2}
    
    mu = np.mean(pred,axis=0)
    std = np.std(pred,axis=0)
    ind = range(3)
    
    fig,ax = plt.subplots(1,1,figsize=(10,10))
    ax.bar(ind,mu,yerr=std, error_kw=error_kw)

    ax.set_ylabel('Fraction of correctly predicted synapses',fontsize=24)
    ax.set_xticks(ind)
    ax.set_xticklabels(('Rand.','Gene match','Gene match +\nspace loc.'))
    ax.set_title(('Postsynaptic exression patterns when combined\n with ' 
                  'spatial localization predicts synaptic partners '),
                fontsize=24)

    if params.fout: plt.savefig(params.fout)
    plt.show()


    
