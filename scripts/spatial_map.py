"""
spatial_map.py

The (r, z) map of synapses. Colors represent the presynaptic (left),
postsynaptic (middle) and gap junction (right) partners. Maps are 
split to show the left and right side of nerve ring. (f) The (φ , z) 
map of synapses.

created: Christopher Brittin
date: 01 November 2018 

"""

import sys
sys.path.append(r'./volumetric_analysis')

import db
import figures.spatialMap as smplt

display = 'Sp1,Sp2,I1,I2,SMN,HMNp,Sa,HMNa'
color_code = './mat/color_code.txt'
neuron_class = './mat/nerve_ring_classes.txt'
_db = 'JSH'
left = './mat/left_nodes.txt'
right = './mat/right_nodes.txt'


def run(fout=None):
    con = db.connect.default(_db)
    cur = con.cursor()
    pre,post = db.mine.syn_cylinder(cur)
    gap = db.mine.gap_cylinder(cur)
    smplt.plot_syn_phiz(_db,pre,post,gap,display,
                        neuron_class,color_code,left,right,fout=fout)


if __name__=='__main__':
    run()
