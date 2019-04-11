"""
cam_cluster.py

Looks for combinations of frequently occuring CAM gene combinations

"""

import sys
sys.path.append('./volumetric_analysis')
import argparse

from cam.expression import Matrix

if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('db',
                        action = 'store',
                        help = 'Database name')

    parser.add_argument('exp',
                        action = 'store',
                        help = 'Expression matrix')
    
    parser.add_argument('-t','--exp_thresh',
                        dest='exp_thresh',
                        action = 'store',
                        required=False,
                        default=1,
                        help="Expression threshold"

    params = parser.parse_args()
    
    M = Matrix(params.exp)
    M.binarize(thresh=params.thresh)
    
    
    

