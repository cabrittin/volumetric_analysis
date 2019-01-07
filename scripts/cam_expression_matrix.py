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
import aux

cam = 'mat/cam_nr_pre_post.csv'
_db = 'JSH'
def write_source(fout,exp):
    (n,m) = exp.E.shape
    genes = dict((exp.genes[g].idx,g) for g in exp.genes)
    cells = dict((exp.nodes[c],c) for c in exp.nodes)
    data = [ ','.join([""] + [genes[i] for i in range(m)])]
    for i in range(n):
        data.append(','.join([cells[i]] + list(map(str,exp.E[i,:].tolist()))))
    aux.write.from_list(fout,data)

def write_source_identical(fout,exp):
    exp.create_links()
    exptag = {}
    for n,tag in exp.exp_link.items():
        if tag not in exptag: exptag[tag] = []
        exptag[tag].append(n)
    aux.write.from_dict(fout,exptag)
    
def run(fout=None,source_data=None):
    mode = 'post'
    con = db.connect.default(_db)
    cur = con.cursor()
    nodes = sorted(db.mine.get_adjacency_cells(cur))
    e = Expression(cur,cam,nodes)
    e.assign_expression_patterns(mode = mode)
    if source_data:
        fsplit = source_data.split('.')
        fident = fsplit[0] + '_' + _db + '.' + fsplit[1]
        write_source(fident,e)
        fident = fsplit[0] + '_' + _db + '_identical.' + fsplit[1]
        write_source_identical(fident,e)
    fig,ax = plt.subplots(1,1,figsize=(5,10))
    expplt.plot_expression_patterns(ax,e,fout=fout)
    plt.show()


if __name__ == '__main__':
    run()
