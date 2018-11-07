"""
connectome.py

Connectome data structures

Author: Christopher Brittin

"""

import networkx as nx
import aux

SCREEN = ['old']

class Connectome:
    def __init__(self,db,neurons):
        self.db = db
        self.size = len(neurons)
        self.neurons = neurons
        self.C = None
        self.E = None
        self.A = None
        
    def update_neurons(self,_neurons):
        self.neurons = _neurons
        self.size = len(_neurons)

    def remove_cells(self,vertices):
        if self.C:self.remove_vertices(self.C,vertices)
        if self.E:self.remove_vertices(self.E,vertices)
        if self.A:self.remove_vertices(self.A,vertices)

    def remove_vertices(self,G,remove):
        remove = set(remove) & set(G.nodes())
        for n in remove:
            G.remove_node(n)

    def group_cells(self,groups,**kwargs):
        if self.C:
            self.C = self.group_vertices(self.C,groups)
        if self.E:
            self.E = self.group_vertices(self.E,groups)
        if self.A:
            self.A = self.group_vertices(self.A,groups)

    def group_vertices(self,G,GROUPS):
        if nx.is_directed(G):
            H = nx.DiGraph()
        else:
            H = nx.Graph()

        for e in G.edges():
            attr = G[e[0]][e[1]]
            if e[0] in GROUPS:
                n1 = GROUPS[e[0]]
            else:
                n1 = e[0]
            if e[1] in GROUPS:
                n2 = GROUPS[e[1]]
            else:
                n2 = e[1]
            if (n1,n2) in H.edges():
                #print n1,n2,H[n1][n2]
                H[n1][n2]['weight'] += attr['weight']
                H[n1][n2]['count'] += attr['count']
            else:
                H.add_edge(n1,n2,weight = attr['weight'],
                           count=attr['count'])
        return H

    def load_chemical(self,synapses,add_poly=False):
        self.C = nx.DiGraph()
        self.load_edges(self.C,self.neurons,synapses,add_poly=add_poly)

    def load_electrical(self,synapses):
        self.E = nx.Graph()
        self.load_edges(self.E,self.neurons,synapses)
                
    def load_edges(self,G,vertices,edges,add_poly=False):
        for e in edges:
            pre = aux.format.rm_brack(e[0])
            if pre not in vertices: continue
            #i_pre = self.neurons[pre]
            _post = list(set(map(aux.format.rm_brack,e[1].split(','))))
            if len(_post) == 1:
                poly = 'S'
            else:
                poly = 'Sp'
            if self.db == 'N2U' and e[4] in ['VC','DC']:
                w = 2*int(e[2]) - 1
            else:
                w = int(e[2])
                
            for post in _post:
                if post not in vertices: continue 
                #i_post = self.neurons[post]
                if not G.has_edge(pre,post):
                    if add_poly:
                        G.add_edge(pre,post,weight=0.0,count=0,S=0,Sp=0)
                    else:
                        G.add_edge(pre,post,weight=0.0,count=0)
                G[pre][post]['weight'] += w
                G[pre][post]['count'] += 1
                if add_poly:
                    G[pre][post][poly] += 1

    def load_adjacency(self,adjacency,directed=False):
        if directed:
            self.A = nx.DiGraph()
        else:
            self.A = nx.Graph()
        for (i,j,weight,imgNum) in adjacency:
            weight = int(weight)
            count = 1
            if self.db == 'N2U' and 'VC' in imgNum:
                weight *=2
                count = 2
            if not self.A.has_edge(i,j):
                self.A.add_edge(i,j,weight=0,count=0)
            self.A[i][j]['weight'] += weight
            self.A[i][j]['count'] += count
            
    def reduce_to_adjacency(self):
        self.C = self._reduce_to_adjacency(self.A,self.C)
        self.E = self._reduce_to_adjacency(self.A,self.E)

    @staticmethod
    def _reduce_to_adjacency(A,H):
        """
        Eliminates nodes and edges in H not found in A
        """
        G = H.copy()
        GnotAnodes = list(set(G.nodes()) - set(A.nodes()))
        G.remove_nodes_from(GnotAnodes)
        edges = [e for e in G.edges()]
        for (a,b) in edges:
            if not A.has_edge(a,b):
                G.remove_edge(a,b)
        return G
            
                
