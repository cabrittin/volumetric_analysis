"""
bar_alt_spliced.py

The number of neurons expressing alternatively spliced genes.

created: Christopher Brittin
date: 01 November 2018
"""

import sys
sys.path.append(r'./volumetric_analysis')
import matplotlib.pyplot as plt
import matplotlib as mpl

import db
from cam.expression import Expression
import cam.cam_plots as expplt
import aux

cam = 'mat/cam_nr_pre_post.csv'
homologs = './mat/homologs.txt'

def run(fout=None):
    nclass = aux.read.into_list2(homologs)    
    mode = 'post'
    con = db.connect.default('N2U')
    cur = con.cursor()
    nodes = sorted(db.mine.get_adjacency_cells(cur))
    e = Expression(cur,cam,nodes)
    e.assign_expression_patterns(mode = mode)
    fig,ax = plt.subplots(1,1,figsize=(10,7))
    expplt.plot_alt_spliced_genes(ax,e,nclass,fout=fout)
    plt.xticks(fontsize=14)
    plt.show()           


if __name__=="__main__":
    run()
