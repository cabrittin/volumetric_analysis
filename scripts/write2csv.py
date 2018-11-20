"""
write2csv.py

Writes adjacency, chemical and electrical networks to edgelist file.

created: Christopher Brittin
date: 01 November 2018 

Synopsis:
  python paper_figures figNum [fout]


Parameters:
  db (str): Database name
  fout (str) : Path to output file 

"""

import argparse
import networkx as nx

from connectome.load import from_db

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('db',
                        action="store",
                        help="Database name")
    parser.add_argument('fout',
                        action="store",
                        help="Path to output file")
    
    params = parser.parse_args()

    C = from_db(params.db,adjacency=True,chemical=True,
                electrical=True, dataType='networkx')

    adj_out = params.fout.replace('.','_adjacency.')
    chem_out = params.fout.replace('.','_chemical.')
    elec_out = params.fout.replace('.','_electrical.')

    nx.write_weighted_edgelist(C.A,adj_out,delimiter=',')
    nx.write_weighted_edgelist(C.C,chem_out,delimiter=',')
    nx.write_weighted_edgelist(C.E,elec_out,delimiter=',')
