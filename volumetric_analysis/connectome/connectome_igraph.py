"""
connectome.py

Connectome data structures

Author: Christopher Brittin

"""

import igraph
import csv
import numpy as np

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
        self.D = None

    def update_cells(self,_neurons):
        self.neurons = _neurons
        self.size = len(_neurons)

    def remove_self_loops(self):
        if self.C: self.C.simplify(multiple=True,loops=True)
        if self.E: self.E.simplify(multiple=True,loops=True)
        
    def remove_cells(self,vertices):
        if self.C:self.C.remove_vertices(vertices)
        if self.E:self.E.remove_vertices(vertices)
        if self.A:self.A.remove_vertices(vertices)

    def group_cells(self,groups,key='group'):
        if self.C:
            self.C.assign_membership(groups,key=key)
            self.C.group_vertices(key)
        if self.E:
            self.E.assign_membership(groups,key=key)
            self.E.group_vertices(key)
        if self.A:
            self.A.assign_membership(groups,key=key)
            self.A.group_vertices(key)
                
    def load_chemical(self,synapses,add_poly=False):
        self.C = Network(directed=True)
        self.C.add_vertices(self.neurons)
        self.load_edges(self.C,self.neurons,synapses,add_poly=add_poly)

    def load_electrical(self,synapses,add_poly=False):
        self.E = Network()
        self.E.add_vertices(self.neurons)
        self.load_edges(self.E,self.neurons,synapses,add_poly=add_poly)
                
    def load_edges(self,G,vertices,edges,add_poly=False):
        eid = -1
        for e in edges:
            pre = aux.format.rm_brack(e[0])
            if pre not in vertices: continue
            #i_pre = self.neurons[pre]
            _post = list(set(map(aux.format.rm_brack,e[1].split(','))))
            if add_poly:
                if len(_post) == 1:
                    poly = 'S'
                else:
                    poly = 'Sp'
            if self.db == 'N2U' and e[4] in ['VC','DC']:
                w = int(e[2])#2*int(e[2]) - 1
            else:
                w = int(e[2])
                
            for post in _post:
                if post not in vertices: continue 
                #i_post = self.neurons[post]
                if not G.are_connected(pre,post):
                    eid += 1
                    G.add_edges([(pre,post)])
                    G.es[eid]['weight'] = 0
                    G.es[eid]['count'] = 0
                    if add_poly:
                        G.es[eid]['S'] = 0
                        G.es[eid]['Sp'] = 0
                _eid = G.get_eid(pre,post)
                G.es[_eid]['weight'] += w
                G.es[_eid]['count'] += 1
                if add_poly:
                    G.es[_eid][poly] += 1
        #G.vs.select(_degree=0).delete()

    def load_adjacency(self,adjacency,directed=False):
        self.A = Network(directed=directed)
        self.A.add_vertices(self.neurons)
        eid = -1
        for (i,j,weight,imgNum) in adjacency:
            weight = weight
            count = 1
            if self.db == 'N2U' and 'VC' in imgNum:
                weight *=1#2
                count = 1
            if not self.A.are_connected(i,j):
                eid += 1
                self.A.add_edges([(i,j)])
                self.A.es[eid]['weight'] = 0
                self.A.es[eid]['count'] = 0
            _eid = self.A.get_eid(i,j)
            self.A.es[_eid]['weight'] += weight
            self.A.es[_eid]['count'] += count

    #def combine_chem_and_elec_(self):
    #    self.D = Network()
    #    self.D.union([self.C,self.E])
            
    def combine_chem_and_elec(self):
        self.D = self.C.copy()
        for e in self.E.es:
            i = self.E.vs[e.source]['name']
            j = self.E.vs[e.target]['name']
            w = 0.5*e['weight']
            c = 0.5*e['count']
            if self.D.are_connected(i,j):
                eid = self.D.get_eid(i,j)
                self.D.es[eid]['weight'] += w
                self.D.es[eid]['count'] += c
            else:
                self.D.add_edge(i,j,weight=w,count=c)

    def reduce_to_adjacency(self):
        if self.C: self.C.reduce_to(self.A)
        if self.E: self.E.reduce_to(self.A)
        
class Network(igraph.Graph):
    def __init__(self,directed=False):
        igraph.Graph.__init__(self,directed=directed)
        self.labels = None

    def get_edge(self,source,target):
        try:
            i = self.vs.select(name=source)[0].index
            j = self.vs.select(name=target)[0].index
            return self.es.select(_source=i,_target=j)[0]
        except IndexError:
            return None

    def get_edge_attr(self,source,target,attr):
        i = self.vs.select(name=source)[0].index
        j = self.vs.select(name=target)[0].index
        w1 = self.es.select(_source=i,_target=j)[attr][0]
        w2 = self.es.select(_source=j,_target=i)[attr][0]
        
    def remove_vertices(self,remove):
        for n in remove:
            for v in self.vs(name=n):
                v.delete()
                
    def assign_membership(self,membership,key='member'):
        if isinstance(membership,list):
            self.assign_membership_list(membership,key=key)
        elif isinstance(membership,dict):
            self.assign_membership_dict(membership,key=key)

    def assign_membership_list(self,membership,key='member'):
        attr = [0]*self.vcount()
        for i in range(len(membership)):
            for j in membership[i]:
                attr[j] = i
        self.vs[key] = attr        

    def assign_membership_dict(self,membership,key='member'):
        attr = [0]*self.vcount()
        for v in self.vs:
            if v['name'] in membership:
                attr[v.index] = membership[v['name']]
            else:
                attr[v.index] = v[key]
        self.vs[key] = attr
        
    def group_vertices(self,key,combine_attrs='first'):
        if isinstance(self.vs[key][0],str):
            vals = list(set(self.vs[key]))
            imap = dict([(vals[i],i) for i in range(len(vals))])
            membership = [imap[v[key]] for v in self.vs]
            self.contract_vertices(membership,combine_attrs=combine_attrs)
        else:
            self.contract_vertices(self.vs[key],combine_attrs=combine_attrs)
        self.simplify(loops=False,combine_edges=sum)
                
    def get_numpy_array(self,directed=False,vertex_order=None,edge_attr=None):
        if vertex_order:
            n = len(vertex_order)
            vhash = dict([(vertex_order[i],i) for i in range(n)])
            A = np.zeros((n,n))
            for e in self.es:
                w = 1
                if edge_attr: w = e[edge_attr]
                A[vhash[e.source]][vhash[e.target]] = w
                if not directed:
                    A[vhash[e.target]][vhash[e.source]] = w
        else:
            n = self.vcount()
            A = np.zeros((n,n))
            for e in self.es:
                w = 1
                if edge_attr: w = e[edge_attr]
                A[e.source][e.target] = w
                if not directed:
                   A[e.target][e.source] = w 
                
        return A

    def symmetrize_matrix(self):
        self.to_undirected(mode="collapse",combine_edges=sum)
        
    def get_neighbors(self,vertex,mode='OUT',vattr='name'):
        _neigh = self.neighbors(vertex, mode=mode)
        return self.vs[_neigh][vattr]
        
    def map_vertex_names(self,vmap):
        G = self.copy()
        for v in G.vs:
            v['name'] = vmap[v['name']]
        return G

    def compute_strength(self,weight='weight'):
        v = range(self.vcount())
        self.vs['out-strength'] = self.strength(v,mode='out',
                                                loops=False,weights=weight)
        self.vs['in-strength'] = self.strength(v,mode='in',
                                               loops=False,weights=weight)
        for e in self.es:
            w = e['weight']
            sout = float(self.vs[e.source]['out-strength'])
            sin = float(self.vs[e.target]['in-strength'])
            e['out-strength'] = w / sout
            e['in-strength'] = w / sin


        
    def reduce_to(self,G):
        eid = []
        mode  = 0
        if not self.is_directed() and G.is_directed():
            mode = 1
        if mode == 1:
            for e in self.es:
                itoj = G.are_connected(e.source,e.target)
                jtoi = G.are_connected(e.target,e.source)
                if not itoj and not jtoi:
                    eid.append(e.index)
        else:
            for e in self.es:
                if not G.are_connected(e.source,e.target):
                    eid.append(e.index)
        self.delete_edges(eid)


    def threshold_edge_greater_than(self,eattr,val):
        es = [e for e in self.es if e[eattr] <= val]
        self.delete_edges(es)
        
    def threshold_edge_less_than(self,eattr,val):
        es = [e for e in self.es if e[eattr] >= val]
        self.delete_edges(es)
