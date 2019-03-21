"""
error_analysis.py

Assess adjacency errors

@author Christopher Brittin
@date 2019-03-20
"""

import sys
sys.path.append('./volumetric_analysis')
import argparse
import networkx as nx

import db
from connectome.load import from_db
from trakem2.graph import Graph

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('db',
                        action="store",
                        help="Database")
    
    parser.add_argument('layer',
                        action="store",
                        help="Layer")
    
    parser.add_argument('sif',
                        action="store",
                        help="Path to sif file")

    parser.add_argument('labels',
                        action="store",
                        help="Path to labels file")

    parser.add_argument('fout',
                        action="store",
                        help="Path to output error file")

    params = parser.parse_args()

    con = db.connect.default(params.db)
    cur = con.cursor()
    adj = db.mine.get_adjacency_from_layer(cur,params.layer)
    G = nx.Graph()
    for (a,b,w) in adj: G.add_edge(a,b,weight='w')

    C = from_db(params.db,adjacency=True,chemical=False,
            electrical=False,dataType='networkx')


    T = Graph(params.sif,params.labels)

    man_not_auto = []
    auto_not_man = []

    for (a,b) in T.edges():
        if not G.has_edge(a,b):
            tmp = (a,b,0)
            if C.A.has_edge(a,b):
                tmp = (a,b,1)
            man_not_auto.append(tmp)

    for (a,b) in G.edges():
        if not T.has_edge(a,b):
            auto_not_man.append((a,b,G[a][b]['weight']))

    print(int(len(man_not_auto)) / T.number_of_edges())
    print(int(len(auto_not_man)) / G.number_of_edges())


    for l in man_not_auto:
        print(l)
