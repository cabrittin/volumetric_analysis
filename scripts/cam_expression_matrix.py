"""
cam_expression_matrix.py

Expression matrix of 35 CAM genes used in this study.

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
    fig,ax = plt.subplots(1,1,figsize=(5,10))
    expplt.plot_expression_patterns(ax,e,fout=fout)
    plt.show()


if __name__ == '__main__':
    run()
