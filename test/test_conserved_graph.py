import sys
sys.path.append('./volumetric_analysis')
sys.path.append('.')
import networkx as nx

from connectome.load import from_db
from mat_loader import MatLoader

def conserved_graph(G,H,deg,nodes):
    for n in nodes:
        try:
            neigh = set([])
            for h in H:
                neigh = neigh.union(set(list(h.neighbors(n))))
            for m in neigh:
                _deg = 0
                for h in H:
                    if h.has_edge(n,m):
                        _deg += 1
                if _deg >= deg:
                    G.add_edge(n,m)
        except:
            print(n)
            pass  
    
    return G

            
if __name__=="__main__":
    M = MatLoader()
    M.load_left()
    M.load_right()
    M.load_lrmap()
    N = from_db('N2U',adjacency=True,chemical=True,
                electrical=True,dataType='networkx')
    N.split_left_right(M.left,M.right)
    N.map_right_graphs(M.lrmap)

    J = from_db('JSH',adjacency=True,chemical=True,
                electrical=True,dataType='networkx')
    J.split_left_right(M.left,M.right)
    J.map_right_graphs(M.lrmap)
    
    A = nx.Graph()
    C = nx.DiGraph()

    A = conserved_graph(A,[N.Al,N.Ar,J.Al,J.Ar],4,M.left) 
    print(J.Al.number_of_edges(),A.number_of_edges())
    
    C = conserved_graph(C,[N.Cl,N.Cr,J.Cl,J.Cr],4,M.left)


