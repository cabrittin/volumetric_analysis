import sys
sys.path.append('./volumetric_analysis')
sys.path.append('.')
import networkx as nx
from tqdm import tqdm

from connectome.load import from_db
from mat_loader import MatLoader

def consensus_graph(G,H,deg,nodes):
    for n in nodes:
        neigh = set([])
        for h in H:
            if not h.has_node(n): continue
            neigh = neigh.union(set(list(h.neighbors(n))))
        for m in neigh:
            _deg = 0
            for h in H:
                if h.has_edge(n,m):
                    _deg += 1
            if _deg >= deg:
                G.add_edge(n,m)
    return G

DOUT = './mat/consensus_graphs/consensus_%s_deg%d.graphml'

if __name__=="__main__":
    M = MatLoader()
    M.load_left()
    M.load_right()
    M.load_lrmap()
    #M.left.remove('ASEL')
    #M.right.remove('ASER')
    N = from_db('N2U',adjacency=True,chemical=True,
                electrical=True,dataType='networkx')
    N.split_left_right(M.left,M.right)
    N.map_right_graphs(M.lrmap)

    J = from_db('JSH',adjacency=True,chemical=True,
                electrical=True,dataType='networkx')
    J.split_left_right(M.left,M.right)
    J.map_right_graphs(M.lrmap)
    
    for deg in tqdm([1,2,3,4],desc="Degree"):
        A = nx.Graph()
        C = nx.DiGraph()
        E = nx.Graph()
        A = consensus_graph(A,[N.Al,N.Ar,J.Al,J.Ar],deg,M.left) 
        C = consensus_graph(C,[N.Cl,N.Cr,J.Cl,J.Cr],deg,M.left)
        E = consensus_graph(E,[N.El,N.Er,J.El,J.Er],deg,M.left)
        
        aout = DOUT %('adj',deg) 
        cout = DOUT %('chem',deg)
        gout = DOUT %('gap',deg)

        nx.write_graphml(A,aout)
        nx.write_graphml(C,cout)
        nx.write_graphml(E,gout)

        


