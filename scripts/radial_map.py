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

import db
import figures.spatialMap as smplt
from mat_loader import MatLoader

_db = 'JSH'


if __name__=='__main__':
    M = MatLoader()
    nclass = M.load_nerve_ring_classes()
    con = db.connect.default(_db)
    cur = con.cursor()
    radial = defaultdict(lambda:defaultdict(list))
    print(nclass)
    for (k,v) in tqdm(nclass.items(),desc="Cells:"):
        loc = db.mine.neuron_cylinder(cur,k)
        if not loc: continue
        for l in loc: radial[v][l[2]].append(l[0])


    
    #for (k,v) in radial.items():
    #    print(k,v)

    
