import sys
sys.path.append(r'../lib/')
sys.path.append(r'..')
import numpy as np
from tabulate import tabulate
import random

#Brittin modules
import db
import aux

class Gene:
    def __init__(self,idx,gene,pre,post,isoforms):
        self.idx = idx
        self.gene = gene
        self.pre = int(pre)
        self.post = int(post)
        self.isoforms = int(isoforms)
        self.expression = []
        self.isoforms = int(isoforms)
        
    def get_idx(self):
        if isinstance(self.idx,list):
            return random.choice(self.idx)
        else:
            return self.idx

class Expression:
    def __init__(self,cur,cam_file,nodes,splice=False):
        self.cur = cur
        self.splice = splice
        self.cam_file = cam_file
        self.cam = aux.read.into_list2(cam_file)
        nodes = sorted(nodes)
        self.nodes = dict([(nodes[i],i) for i in range(len(nodes))])
        self.M = len(nodes)
        self.exp = dict([('pre',np.zeros([self.M,2])),
                         ('post',np.zeros([self.M,2]))
                         ])
        self.load_genes()
        self.combined_exp = None
             
    def load_genes(self):
        self.genes = {}
        self.gene_idx = {}
        idx = 0
        for (name,pre,post,isoform) in self.cam:
            if self.splice:
                _idx = []
                for j in range(int(isoform)):
                    _idx.append(idx)
                    self.gene_idx[idx] = name
                    idx += 1
                self.genes[name] = Gene(_idx,name,pre,post,isoform)
            else:
                self.genes[name] = Gene(idx,name,pre,post,isoform)        
                self.gene_idx[idx] = name
                idx += 1
        self.N = idx        

    def assign_expression_patterns(self,mode='post'):
        self.E = np.zeros([self.M,self.N],dtype=int)
        for n,idx in self.nodes.items():
            _genes = db.mine.get_cell_genes(self.cur,n)
            for g in _genes:
                if (mode == 'post') and (g in self.genes) and (self.genes[g].post):
                    self.genes[g].expression.append(idx)
                    jdx = self.genes[g].get_idx()
                    self.E[idx,jdx] = 1
                elif (mode == 'pre') and (g in self.genes) and (self.genes[g].pre):
                    self.genes[g].expression.append(idx)
                    jdx = self.genes[g].get_idx()
                    self.E[idx,jdx] = 1
                   
        self.Diff = np.zeros([self.M,self.M])
        for n,idx in self.nodes.items():
            for m,jdx in self.nodes.items():
                self.Diff[idx,jdx] = np.sum(abs(self.E[idx,:] - self.E[jdx,:]))
                
        self.create_links()
        
    def create_links(self):
        self.exp_link = {}
        for n,idx in self.nodes.items():
            self.exp_link[n] = self.assign_tag(self.E[idx,:])
            
    @staticmethod
    def assign_tag(bits):
        return ''.join(map(str,bits))


    def compute_difference(self,n1,n2):
        return self.Diff[self.nodes[n1],self.nodes[n2]]


    def set_combined_expression(self,neurons):
        tmp = np.zeros([1,self.N])
        for n in neurons:
            tmp += self.E[self.nodes[n]]
        tmp[tmp>1] = 1
        self.combined = tmp

    def compute_combined_difference(self,n):
        idx = self.nodes[n]
        return np.sum(abs(self.combined - self.E[idx,:]))


    def gene_per_neuron_count(self):
        data = np.zeros(self.N)
        for idx in xrange(self.M):
            s = int(np.sum(self.E[idx,:]))
            data[s] += 1
        jdx = np.nonzero(data)[0][-1]+1
        return data[:jdx]

    def neuron_per_gene_count(self):
        return np.sum(self.E,axis=0)

    def isoform_per_gene_count(self):
        data = np.zeros(100)
        for g in self.genes:
            s = self.genes[g].isoforms
            data[s] += 1
        jdx = np.nonzero(data)[0][-1]+1
        return data[:jdx]

    def alt_splice_bilateral_dist(self,nclass):
        self.alt_gene_count = {}
        data = []
        for nc in nclass:
            if len(nc) < 3: continue
            if nc[1] not in self.nodes: continue
            if nc[2] not in self.nodes: continue
            alt_genes = []
            for n in nc[1:]:
                idx = self.nodes[n]
                jdx = np.nonzero(self.E[idx,:])[0]
                for j in jdx:
                    g = self.gene_idx[j]
                    if self.genes[g].isoforms > 1:
                        if g not in self.alt_gene_count: self.alt_gene_count[g] = 0
                        self.alt_gene_count[g] += 1
                        alt_genes.append(g)
            data.append([nc[0],sorted(list(set(alt_genes)))])

        return data
    
def get_lr_discrepancies(C,lrd):
    neurons,synapses = 0,0
    discreps = []
    for n in C.A.nodes():
        try:
            for m in C.C.neighbors(n):
                _synapses = 0
                mh = lrd[m]
                if m != mh and not C.C.has_edge(n,mh) and C.A.has_edge(n,mh):
                    discreps.append([C.db,n,m,mh,C.A.has_edge(n,mh),C.C.has_edge(n,mh),
                                     C.C.has_edge(n,m)])
                    _synapses += 1
                    synapses += 1
            if _synapses: neurons += 1
        except:
            pass

    frac_nodes = float(neurons)/C.A.number_of_nodes()
    frac_edges = float(synapses)/C.C.number_of_edges()

    print(tabulate(discreps,
                   headers=['db','Pre (n)','Post (m)',"Homolog post (m')",
                                     "n|--|m'","n-->m'","n-->m"],
                   tablefmt='orgtbl'))
                   
    print('\n-----------------------\n')
    print(tabulate([[C.db,neurons,frac_nodes,synapses,frac_edges]],
                   headers=['db','# nodes','% nodes','# edges','% edges'],
                   tablefmt='orgtbl'))  

    return [frac_nodes,frac_edges]
