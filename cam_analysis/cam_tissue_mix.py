"""
cam_tissue_mix.py

calculate the average diversity of cell clusters

"""

import sys
sys.path.append('./volumetric_analysis')
sys.path.append('.')
import argparse
import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt
import matplotlib as mpl
import networkx as nx
import seaborn as sns
from collections import defaultdict

from cam.expression import Matrix
from mat_loader import MatLoader
from connectome.load import from_db
import aux

mpl.rcParams['xtick.labelsize'] = 16 
mpl.rcParams['xtick.labelsize'] = 14 

cam_class = ['cad','igsf','lrr','nrx','all']
METRIC = 'pearsonr'
NODE_SCREEN = ['NSM','MC']

if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
 
    params = parser.parse_args()
    
    ML = MatLoader()
    ML.load_lrmap()
    #nodes = sorted(ML.load_reduced_nodes())
    nodes = sorted(ML.load_all_tissue())
    neurons = sorted(ML.load_reduced_nodes()) + NODE_SCREEN
    
    for c in cam_class:
        camclass = ML.load_cam_class(METRIC,c)
        k = max(camclass.values()) + 1
        num = np.zeros(k)
        den = np.zeros(k)
        for (i,j) in camclass.items():
            den[j] += 1
            if i in neurons: num[j] += 1

        p = num/den
        plog = np.log2(p)
        plog[plog==-np.inf] = 0
        q = 1 - p
        qlog = np.log2(q)
        qlog[qlog==-np.inf] = 0
        H = -p*plog -q*qlog
        hmu = np.mean(H)
        hstd = np.std(H)
        print(H)
        print(hmu,hstd)

    

    
