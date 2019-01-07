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

import argparse
import numpy as np

import db
from cam.expression import Expression
from connectome.load import from_db
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

    params = parser.parse_args()

    _end = 500
    if params.db == 'N2U': _end = 325
    
    mode = 'post'
    con = db.connect.default(params.db)
    cur = con.cursor()
    nodes = sorted(db.mine.get_adjacency_cells(cur))
    e = Expression(cur,cam,nodes)
    e.assign_expression_patterns(mode = mode)   

    C = from_db(params.db,adjacency=True,chemical=True,
                electrical=True,dataType='networkx')
    C.remove_self_loops()
    C.reduce_to_adjacency()    
    
    wbe = cam_lus.wbe(e,C)
    wbe_data = wbe.get_data()
    aux.write.from_list('results/'+params.db+'_wbe.csv',wbe_data) 
    
    sbe = cam_lus.sbe(e,end=_end)       
    sbe_data = sbe.get_data()
    aux.write.from_list('results/'+params.db+'_sbe.csv',sbe_data)
        
    e.splice = True
    e.load_genes()
    ie = cam_lus.ie(e,C,iters=ie_iter,mode=mode)
    np.savetxt('results/'+params.db+'_ie.csv',ie,delimiter=',',fmt='%1.4f')
