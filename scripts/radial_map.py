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


import db
import figures.spatialMap as smplt
from mat_loader import MatLoader

_db = 'JSH'

def run(fout=None):
    M = MatLoader()
    nclass = M.load_nerve_ring_classes()
    con = db.connect.default(_db)
    cur = con.cursor()
    loc = db.mine.neuron_cylinder(cur)
    print(loc)
    for l in loc:
        print(l)

if __name__=='__main__':
    run()
