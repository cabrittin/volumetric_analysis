"""
dist_neurons_expr_gene.py

Cumulative distribution of the number of NR neurons that 
express a given gene.

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

cam = 'mat/cam_nr_pre_post.csv'


def run(fout=None):
    mode = 'post'
    con = db.connect.default('N2U')
    cur = con.cursor()
    nodes = sorted(db.mine.get_adjacency_cells(cur))
    e = Expression(cur,cam,nodes)
    e.assign_expression_patterns(mode = mode)
    fig,ax = plt.subplots(1,1,figsize=(10,7))
    expplt.plot_neurons_expressing_genes(ax,e,fout)
    plt.show()           


if __name__=="__main__":
    run()


