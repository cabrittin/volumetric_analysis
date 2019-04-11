"""
mk_wb_exp_matrix.py

Outputs the expression matrix for WormBase data

"""

import sys
sys.path.append('./volumetric_analysis')
import argparse

import db
import aux
from cam.expression import Expression
from mat_loader import MatLoader

FOUT = 'wb_exp_matrix.csv'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('cam',
                        action = 'store',
                        help = "Path to file with cam genes")

    parser.add_argument('dout',
                        action = 'store',
                        help = "Path to output directory")

    params = parser.parse_args()

    M = MatLoader()
    con = db.connect.default("N2U")
    cur = con.cursor()
    nodes = sorted(db.mine.get_adjacency_cells(cur))
    e = Expression(cur,params.cam,nodes)
    e.assign_expression_patterns(mode='post')
    
    data = e.expression_list()

    data = [d + [100] for d in data]

    fout = params.dout + FOUT

    aux.write.from_list(fout,data)
