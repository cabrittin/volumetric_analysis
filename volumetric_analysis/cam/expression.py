"""
expression.py

Classes for analyzing the CAM expression data

"""

import sys
import numpy as np
from tabulate import tabulate
import random

#Brittin modules
import db
import aux

class Matrix:
    """
    Class to hold expression matrix
    
    ...

    Attributes
    ----------
    fname : str
        path to expression matrix csv cam_file
        Row format [gene,cell,exression level]
    
    
    """
    
    def __init__(self,_fname,cells=None):
        data = aux.read.into_list2(_fname)
        _genes = []
        _cells = []
        for d in data:
            if d[0] not in _genes: _genes.append(d[0])
            if d[1] not in _cells: _cells.append(d[1])

        _genes = sorted(_genes)
        _cells = sorted(_cells)
        if cells: _cells = cells
        self.m = len(_cells)
        self.n = len(_genes)
        self.cells = dict([(_cells[i],i) for i in range(self.m)])
        self.cell_idx = {y:x for x,y in self.cells.items()}
        self.genes = dict([(_genes[i],i) for i in range(self.n)])
        self.genes_idx = {y:x for x,y in self.genes.items()}
        
        self.M = np.zeros((self.m,self.n))

        for [g,c,val] in data:
            if c not in _cells: continue
            self.M[self.cells[c],self.genes[g]] = float(val)

    def binarize(self,thresh=1):
        """
        Creates a binary expression matrix

        Paramters:
        thresh : float
            Values below thresh set to zero. Values above thresh 
            set to 1.
        """
        self.E = np.copy(self.M)
        self.E[self.E<thresh] = 0
        self.E[self.E>0] = 1
    
    def difference_matrix(self):
        """
        Computes the diffence between expression patterns between cells
        """
        self.Diff = np.zeros((self.m,self.m))
        for n,idx in self.cells.items():
            for m,jdx in self.cells.items():
                self.Diff[idx,jdx] = np.sum(abs(self.E[idx,:] - self.E[jdx,:]))
   
    def compute_difference(self,n1,n2):
        """
        Computes gene expression difference between cells n1 and n2

        Parameters
        ----------
        n1 : str
          Cell 1
        n2 : str
          Cell 2
        
        Returns 
        -------
        The modified self.Diff
        """
        return self.Diff[self.cells[n1],self.cells[n2]]


    def cell_expression(self,cell):
        idx = self.cells[cell]
        row = np.where(self.M[idx,:] > 0)[0]
        return [self.genes_idx[i] for i in row]



class Gene:
    """
    Class to represent each gene
    
    ...

    Attributes
    ----------
    idx : int or list 
      Index of gene. Primarily used to distinguish between different isoforms.
      If isoforms, then idx is a list of isoform indicies
    gene : str
      Gene name
    pre  : int
      1 if expressed in presynaptic neurons. Otherwise, 0.
    post : int
      1 if expressed in postsynaptic neurons. Otherwise, 0.
    isoforms : int
      Number of isoforms
    expression : list 
      List of cells in which gene is expressed
      
    Methods
    -------
    get_idx(): returns index. If an isoform, returns a random index.

    """  
    def __init__(self,idx,gene,pre,post,isoforms):
        self.idx = idx
        self.gene = gene
        self.pre = int(pre)
        self.post = int(post)
        self.isoforms = int(isoforms)
        self.expression = []
        
    def get_idx(self):
        """
        Returns index. If an isoform, returns a random index.
        """
        if isinstance(self.idx,list):
            return random.choice(self.idx)
        else:
            return self.idx

class Expression:
    """
    Class to represent the expression data.

    ...
    
    Attributes
    ----------
    cur : MySQLdb.cursor()
      MySQL cursor
    splice : bool
      If true, use alternative splicing.
    cam_file : str
      Path to file with the CAM expression data.
    cam : 2D list
      Data in cam_file
    nodes : dict
      Dictionary that maps cell names to index.
    nodes_idx : dict
      Array of nodes ordered by index
    M : int
      Number of cells
    N : int
      Number of genes
    genes : dict
      Dictionary of gene name to Gene object
    genes_idx : dict
      dictionary of gene index to gene name
    E : numpy array
      Expression matrix
    exp_link : dict
      Dictionary that maps cell to binary gene expression pattern
    Diff : [MxM] array
      Difference between gene expression patterns
    combined : 1D array
      The combined expression of specified cells
    
    Methods
    -------
    load_genes()
       Loads genes from the self.cam_file.
    assign_expression_patterns(mode='post')
       Generates the expression matrix. Mode specifies whether
       to look at postsynapti ('post') for presynaptic ('pre')
       expression. 
    create_links():
       Links the cell to the binary gene expression pattern.
    assign_tag(bits)
       Converts binary expression into a string.
    compute_difference(cell1,cell2)
       Computes gene expression difference between cell1 and cell2
    set_combined_expression(neurons)
       Creates the combined expression patterns of cells in list neuron.
    compute_combined_difference(n)
       Computes the difference between expression of cell n and the 
       combined expression (self.combined)
    gene_per_neuron_count()
       Counts the number genes expressed in each neuron. 
    neuron_per_gene_count()
       Returns the number of neurons expressing a given gene
    isoform_per_gene_count()
       Returns the number of isorforms per gene
    expression_list()
       Returns list of nonzero inputs with [gene,cell]

    """    
    def __init__(self,cur,cam_file,nodes,splice=False):
        """
        Parameters
        ----------
        cur : MySQLdb.cursor()
          MySQL cursor
        cam_file : str
          Path to cam expresssion file
        nodes : list
          list of cell names
        splice : bool (optional)
         If true, use isoform expression. (default = 1)

        """        
        self.cur = cur
        self.splice = splice
        self.cam_file = cam_file
        self.cam = aux.read.into_list2(cam_file)
        nodes = sorted(nodes)
        self.nodes = dict([(nodes[i],i) for i in range(len(nodes))])
        self.nodes_idx = sorted(nodes)
        self.M = len(nodes)
        self.exp = dict([('pre',np.zeros([self.M,2])),
                         ('post',np.zeros([self.M,2]))
                         ])
        self.load_genes()
        self.combined_exp = None
             
    def load_genes(self):
        """
        Loads genes from the self.cam_file.
        """
        
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
        """
        Generates the expression matrix. Mode specifies whether
        to look at postsynapti ('post') for presynaptic ('pre')
        expression.  

        Parameters
        ----------
        mode : str (default='post')
           Specifies either 'pre' or 'post' synaptic expression.
        """
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
        """
        Links the cell to the binary gene expression pattern.
        """
        self.exp_link = {}
        for n,idx in self.nodes.items():
            self.exp_link[n] = self.assign_tag(self.E[idx,:])
            
    @staticmethod
    def assign_tag(bits):
        """
        Converts binary expression into a string
        
        Parameters
        ----------
        bits : numpy row
          Row of 1/0 for gene expression.
        
        """
        return ''.join(map(str,bits))


    def compute_difference(self,n1,n2):
        """
        Computes gene expression difference between cells n1 and n2

        Parameters
        ----------
        n1 : str
          Cell 1
        n2 : str
          Cell 2
        
        Returns 
        -------
        The modified self.Diff
        """
        return self.Diff[self.nodes[n1],self.nodes[n2]]


    def set_combined_expression(self,neurons):
        """
        Creates the combined expression patterns of cells in list neuron.
        
        Parameters
        ----------
        neurons : list
          List of cells to create the combined expression
        
        """
        tmp = np.zeros([1,self.N])
        for n in neurons:
            tmp += self.E[self.nodes[n]]
        tmp[tmp>1] = 1
        self.combined = tmp

    def compute_combined_difference(self,n):
        """
        Computes the difference between expression of cell n and the 
        combined expression (self.combined) 

        Parameters
        ----------
        n : str
          Cell name
        
        Returns
        ---------
        Sum of binary expression pattern
        
        """
        idx = self.nodes[n]
        return np.sum(abs(self.combined - self.E[idx,:]))


    def gene_per_neuron_count(self):
        """
        Return the counts the number genes expressed in each neuron.
        
        """
        data = np.zeros(self.N)
        for idx in range(self.M):
            s = int(np.sum(self.E[idx,:]))
            data[s] += 1
        jdx = np.nonzero(data)[0][-1]+1
        return data[:jdx]

    def neuron_per_gene_count(self):
        """
        Returns the number of neurons expressing a given gene
        """
        return np.sum(self.E,axis=0)

    def isoform_per_gene_count(self):
        """
        Returns the number of isorforms per gene
        """
        data = np.zeros(100)
        for g in self.genes:
            s = self.genes[g].isoforms
            data[s] += 1
        jdx = np.nonzero(data)[0][-1]+1
        return data[:jdx]

    def alt_splice_bilateral_dist(self,nclass):
        """
        Returns the distribution of alternative spliced genes
        in bilaterally symmetric neurons
        
        Parameters:
        nclass : dict
          Specified the bilateral symmetry
        """
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
   
   
    def expression_list(self):
        """
        Returns list of nonzero elements of expression matrix
        """
        nz = np.nonzero(self.E)
        data = []
        for i in range(len(nz[0])):
            data.append([self.gene_idx[nz[1][i]],self.nodes_idx[nz[0][i]]])

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
