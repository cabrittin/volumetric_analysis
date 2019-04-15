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
sys.path.append('.')
import argparse
import numpy as np

import db
from cam.expression import Matrix 
from connectome.load import from_db
from connectome import consensus
import cam.cam_lus as cam_lus
import aux
from mat_loader import MatLoader

cam = 'mat/cam_isoforms.txt'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('matrix',
                        action = 'store',
                        help = 'Path to matrix file')

    parser.add_argument('deg',
                        action = 'store',
                        type= int,
                        help = 'Conserved degree')

    parser.add_argument('fout',
                        action='store',
                        help = 'Path to output file')

    parser.add_argument('--ie_iter',
                        dest = 'ie_iter',
                        action = 'store',
                        required = False,
                        default = 1000,
                        type = int,
                        help = "Number IE iterations")
    
    params = parser.parse_args()
    
    M = MatLoader()
    M.load_left()
    D = M.load_consensus_graphs(params.deg)
    S = M.load_consensus_chemical_synapse(params.deg)
    G = M.load_consensus_gap_junctions(params.deg)
    
    nodes = sorted(D.A.nodes()) 

    e = Matrix(cam,params.matrix)
    e.load_genes()
    e.load_cells(sorted(D.A.nodes()))
    e.assign_expression()
    e.binarize()
    e.difference_matrix()
    
    wbe = cam_lus.wbe(e,D,cells=M.left)
    wbe_data = wbe.get_data()
    tmp = params.fout.replace('.','_wbe.')
    aux.write.from_list(tmp,wbe_data) 

    
    sbe = cam_lus.sbe(e,S,G)       
    sbe_data = sbe.get_data()
    tmp = params.fout.replace('.','_sbe.')
    aux.write.from_list(tmp,sbe_data)
    
    ie = cam_lus.ie(e,D,cells=M.left,iters=params.ie_iter)
    tmp = params.fout.replace('.','_ie.')
    np.savetxt(tmp,ie,delimiter=',',fmt='%1.4f')
