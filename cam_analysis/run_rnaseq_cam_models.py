"""
run_cam_models.py

Runs the CAM models:
  WBE: Whole-cell binary expression
  SBE: Sub-cell binary expression
  IE: Isoform expression

created: Christopher Brittin
date: 01 November 2018 

Synopsis:
  python run_cam_models db_name 


Parameters:
  db (str) : database name

"""

import sys
sys.path.append('./volumetric_analysis')
import argparse
import numpy as np

import db
from cam.expression import Matrix 
from connectome.load import from_db
import cam.cam_lus as cam_lus
import aux

cam = 'mat/cam_nr_pre_post.csv'
ie_iter = 1000

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('db',
                        action = 'store',
                        help = 'Database name')
    
    parser.add_argument('matrix',
                        action = 'store',
                        help = 'Path to matrix file')

    parser.add_argument('fout',
                        action='store',
                        help = 'Path to output file')

    params = parser.parse_args()

    _end = 500
    if params.db == 'N2U': _end = 325
    
    mode = 'post'
    con = db.connect.default(params.db)
    cur = con.cursor()
    nodes = sorted(db.mine.get_adjacency_cells(cur))

    e = Matrix(params.matrix)
    e.binarize()
    e.difference_matrix()

    e.cur = cur
    e.nodes = e.cells.keys()
    remove = set(nodes) - set(e.cells.keys())

    #e = Expression(cur,cam,nodes)
    #e.assign_expression_patterns(mode = mode)   

    C = from_db(params.db,adjacency=True,chemical=True,
                remove=remove,electrical=True,dataType='networkx')
    C.remove_self_loops()
    C.reduce_to_adjacency()    

    
    wbe = cam_lus.wbe(e,C)
    wbe_data = wbe.get_data()
    tmp = params.fout.replace('.','_wbe.')
    aux.write.from_list(tmp,wbe_data) 

    
    sbe = cam_lus.sbe(e,end=_end)       
    sbe_data = sbe.get_data()
    tmp = params.fout.replace('.','_sbe.')
    aux.write.from_list(tmp,sbe_data)
        
    #e.splice = True
    #e.load_genes()
    #ie = cam_lus.ie(e,C,iters=ie_iter,mode=mode)
    #np.savetxt('results/'+params.db+'_ie.csv',ie,delimiter=',',fmt='%1.4f')
