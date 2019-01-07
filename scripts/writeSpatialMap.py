"""
writeSpatialMap.py

Writes data for the spatial map analysis.

created: Christopher Brittin
date: 01 November 2018 

"""

import argparse

import db
import aux


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('db',
                        action="store",
                        help="Database name")
    parser.add_argument('fout',
                        action="store",
                        help="Path to output file")
    
    params = parser.parse_args()

    fsplit = params.fout.split('.')
    
    con = db.connect.default(params.db)
    cur = con.cursor()
    pre,post = db.mine.syn_cylinder(cur)
    gap = db.mine.gap_cylinder(cur)

    preout = fsplit[0] + '_pre.' + fsplit[1]
    postout = fsplit[0] + '_post.' + fsplit[1]
    gapout = fsplit[0] + '_gap.' + fsplit[1]

    aux.write.from_list(preout,pre)
    aux.write.from_list(postout,post)
    aux.write.from_list(gapout,gap)
    
