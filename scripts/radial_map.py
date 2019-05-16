"""
radial_map.py

Plots (r, z) map of neurites. Colors represent cell classes.

created: Christopher Brittin
date: 01 November 2018 

"""

import sys
sys.path.append('.')
sys.path.append(r'./volumetric_analysis')
import matplotlib.pyplot as plt
from collections import defaultdict
from tqdm import tqdm
import pandas as pd
import numpy as np
import seaborn as sns

import db
import figures.spatialMap as smplt
from mat_loader import MatLoader

_db = 'JSH'

remap = {'I1':'I1','Sp1':'Sp','Sp2':'Sp','Sp':'Sp',
        'HMNp':'HMN','HMNa':'HMN','HMN':'HMN','I2':'I2',
        'Sa':'Sa','SMN':'SMN'}

HUE_ORDER = ['Sp','Sa','I1','I2','HMN','SMN']
PALETTE = {'Sp':'#910000','Sa':'#FE6F00','I1':'#FF00FF',
        'I2':'#1700FF','HMN':'#1D6500','SMN':'#2aff00'}

XSCALE = 0.09
RSCALE = 0.005

if __name__=='__main__':
    M = MatLoader()
    M.load_left()
    M.load_right()
    nclass = M.load_nerve_ring_classes()
    con = db.connect.default(_db)
    cur = con.cursor()
    #radial = defaultdict(lambda:defaultdict(list))
    data = []
    for (k,v) in tqdm(nclass.items(),desc="Cells:"):
        if k not in M.left: continue
        if v not in remap: continue
        rv = remap[v]
        loc = db.mine.neuron_cylinder(cur,k)
        if not loc: continue
        #for l in loc: radial[v][l[2]].append(l[0])
        for l in loc: 
            if l[2] > 160: continue
            data.append([k,rv,l[2]*XSCALE,l[0]*RSCALE])

    #data = []    
    #for (cell,rad) in radial.items():
    #    for (sect,dist) in rad.items():
    #        mu = np.mean(dist)
    #        std = np.std(dist)
    #        rmin = np.min(dist)
    #        rmax = np.max(dist)
    #        pct = np.percentile(dist,q=np.array([10,90]))
    #        data.append([cell,sect,mu,std,rmin,rmax,pct[0],pct[1]])


    #df = pd.DataFrame(data, columns = ['cell_class', 'sect_num','mu','std','rmin','rmax','p10','p90']) 
    df = pd.DataFrame(data, columns = ['Cell','Cell class','Axial position','Radial distance'])

    sns.lineplot(x='Axial position',y='Radial distance',hue='Cell class',data=df,
            hue_order=HUE_ORDER,palette=PALETTE,ci='sd')
    plt.show()

    
