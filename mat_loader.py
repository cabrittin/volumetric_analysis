"""
Class to load mat files

@author Christopher Brittin
@data 09 April 2019
"""

import sys
sys.path.append('./volumetric_analysis')
import networkx as nx

import aux
from connectome import consensus


class Connectome:
    def __init__(self,A,C,E):
        self.A = A
        self.C = C
        self.E = E
        self.neurons = A.nodes()

class MatLoader:
    def __init__(self,mat_files='mat/mat_files.txt'):
        self.load_mat(mat_files=mat_files)
        self.cam = 'mat/cam_isoforms.txt' 
    
    def load_mat(self,mat_files='mat/mat_files.txt'):
        self.mat_files = mat_files
        self.mat = aux.read.into_dict(self.mat_files)
    
    def load_left(self):
        self.left = aux.read.into_list(self.mat['left'])

    def load_right(self):
        self.right = aux.read.into_list(self.mat['right'])

    def load_lrmap(self):
        self.lrmap = aux.read.into_lr_dict(self.mat['lrmap'])
    
    def load_lrmap2(self):
        self.lrmap = aux.read.into_lr_dict(self.mat['lrmap2'])
    
    def load_cam_genes(self):
        self.genes = aux.read.into_list(self.mat['genes'])

    def load_isoforms(self):
        self.isoforms = aux.read.into_dict(self.mat['isoforms'])

    def load_consensus_graphs(self,deg):
        A = nx.read_graphml(self.mat['consensus']%('adj',deg))
        C = nx.read_graphml(self.mat['consensus']%('chem',deg))
        E = nx.read_graphml(self.mat['consensus']%('gap',deg))
        return Connectome(A,C,E)
    
    def load_consensus_master_graphs(self,deg):
        A = nx.read_graphml(self.mat['consensus_master']%('adj',deg))
        C = nx.read_graphml(self.mat['consensus_master']%('chem',deg))
        E = nx.read_graphml(self.mat['consensus_master']%('gap',deg))
        return Connectome(A,C,E)

    def load_consensus_chemical_synapse(self,deg):
        fin = self.mat['consensus_chemical']%deg
        return consensus.convert_xml_to_synapse(fin)

    def load_consensus_chemical_post_synapse(self,deg):
        fin = self.mat['consensus_post']%deg
        return consensus.convert_xml_to_synapse(fin)

    def load_consensus_gap_junctions(self,deg):
        fin = self.mat['consensus_gap']%deg
        return consensus.convert_xml_to_synapse(fin)

    def load_gene_sig_pre(self,deg):
        fin = self.mat['gene_sig_pre']%deg
        return aux.read.into_dict2(fin)

    def load_gene_sig_post(self,deg):
        fin = self.mat['gene_sig_post']%deg
        return aux.read.into_dict2(fin)

    def load_gene_sig_gap(self,deg):
        fin = self.mat['gene_sig_gap']%deg
        return aux.read.into_dict2(fin)
    
    def load_gene_sig_graph(self,deg):
        pre = nx.read_graphml(self.mat['gene_sig_graph']%('pre',deg))
        post = nx.read_graphml(self.mat['gene_sig_graph']%('post',deg))
        return pre,post

    def load_nerve_ring_classes(self):
        data = aux.read.into_dict(self.mat['nrclass'])
        _data = {}
        for d in data: _data[d] = data[d].strip()
        return _data

    def load_reduced_nodes(self):
        return aux.read.into_list(self.mat['reduced_nodes'])

    def load_sa(self,db):
        fin = self.mat['lengths']%(db.lower())
        l = aux.read.into_dict2(fin)
        sa = {k:v[0] for (k,v) in l.items()}
        return sa
    
    def load_cam_class(self,metric,camtype):
        fin = self.mat['cam_class']%(metric,camtype)
        return {k:int(v) for (k,v) in aux.read.into_dict(fin).items()}

if __name__=="__main__":
    M = MatLoader()

    print(M.mat)
    print('Left nodes:')
    M.load_left()
    for n in M.left: print('\t'+n)
    print('Right nodes:')
    M.load_right()
    for n in M.right: print('\t'+n)
    M.load_cam_genes()
    for n in M.genes: print('\t'+n)
    
