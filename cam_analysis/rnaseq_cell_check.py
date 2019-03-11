"""
rnaseq_cell_check.py

Check which cells are not included in the rnaseq data

created: Christopher Brittin
data: 09 March 2019
"""

import sys
sys.path.append('./volumetric_analysis')
import argparse

import db
import aux
from connectome.load import from_db

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('db',
                        action="store",
                        help="Database name")
    
    parser.add_argument('fin',
                        action="store",
                        help="Path to group file ")

    params = parser.parse_args()

    con = db.connect.default(params.db)
    cur = con.cursor()

    sql = "select distinct(cellsname) from cells"
    cur.execute(sql)
    rna = set([c[0] for c in cur.fetchall()])
    con.close()

    _group = aux.read.into_map(params.fin)
    group = {'map':_group,'key':'class'}

    C = from_db('N2U',adjacency=True,chemical=False,
                electrical=False,group=group,dataType='networkx')

    #groups = set([d[0] for d in aux.read.into_list2(params.fin)])
    groups = set(C.A.nodes())
    print('Cells in rnaseq not in volumetric: %d' %len(rna-groups))
    print(sorted(rna - groups))

    print('Cells in volumetric not in rnaseq: %d' %len(groups-rna))
    print(sorted(groups - rna))
    

