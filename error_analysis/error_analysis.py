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
import aux

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
    for (a,b,w) in adj: G.add_edge(a,b,weight=w)

    C = from_db(params.db,adjacency=True,chemical=False,
            electrical=False,dataType='networkx')


    T = Graph(params.sif,params.labels)
    
    print(T.number_of_nodes())
    print(G.number_of_nodes())
    print(set(T.nodes()) - set(G.nodes()))
    man_not_auto = []
    auto_not_man = []

    for (a,b) in T.edges():
        [a,b] = sorted([a,b])
        if not G.has_edge(a,b):
            tmp = (a,b,0)
            if C.A.has_edge(a,b):
                tmp = (a,b,1)
            man_not_auto.append(tmp)

    for (a,b) in G.edges():
        [a,b] = sorted([a,b])
        if not T.has_edge(a,b):
            auto_not_man.append((a,b,G[a][b]['weight']))
    

    man_not_auto.sort(key=lambda x: x[0])
    auto_not_man.sort(key=lambda x: x[0])

    print("Number of auto scored adjacencies: %d" %G.number_of_edges())
    print("Fraction manual not scored in layer: ",int(len(man_not_auto)) / T.number_of_edges())
    print("Fraction auto not scored by manual: ",int(len(auto_not_man)) / G.number_of_edges())

    manual_results = []
    print("Score the reason why the adjacency was not scored.")
    print("Key: (0) should have been scored (1) manual not correct; (2) poor segmentation")

    idx = 0
    N = len(man_not_auto)
    for l in man_not_auto:
        idx += 1
        print("%d/%d: Cells: (%s,%s); Eventually scored: %d" %(idx,N,l[0],l[1],l[2]))
        code = input("Input reason for noauto score: ")
        manual_results.append([params.layer] + list(l) + [int(code)])
    
    auto_results = []
    print("Score the reason why the adjacency was not scored.")
    print("Key: (0) should have been scored (1) auto not correct; (2) poor segmentation")

    idx = 0
    N = len(auto_not_man)
    for l in auto_not_man:
        idx += 1
        print("%d/%d: Cells: (%s,%s)" %(idx,N,l[0],l[1]))
        code = input("Input reason for noauto score: ")
        auto_results.append([params.layer] + list(l) + [int(code)])
        
    
    aux.write.from_list(params.fout + params.layer + '_manual_error.csv',manual_results)
    aux.write.from_list(params.fout + params.layer + "_auto_error.csv",auto_results)

