"""
cam_lus.py

Submodule for testing different CAM expression models
via the LUS.

Author: Christopher Brittin
Created: 08 February 2018

"""
from tqdm import tqdm
import numpy as np

import db

class Model(object):
    """
    Class for doing the LUS analysis on the different models
    
    Attributes
    ----------
    lus : dict
     LUS scores for cells
    modelName : str
     Name of model being analyzed
    
    Methods
    -------
    add_lust(_neuron,_lus,_mode)
      Adds an LUS score for _neuron
    get_lus(_mode=None,thresh=0)
      Returns the LUS scores abouve thresh
    get_data()
      Returns model results. 
      Row format: cell_name,pre_lus,gap_lus
    
    """
    def __init__(self,_modelName):
        """
        Parameters
        ---------
        _modelName : str
         Name of model

        """
        self.lus = {}
        self.modelName = _modelName

    def add_lus(self,_neuron,_lus,_mode):
        """
        Adds teh LUS score for _neuron
        
        Parameters
        ----------
        _neuron : str
          Name of cell
        _lus : float
          LUS score for the cell
        _mode : str
         Either 'pre', 'post' or 'gap' synaptic mode
        """
        if _neuron not in self.lus:
            self.lus[_neuron] = {'gap':-1,'pre':-1,'post':-1}
        self.lus[_neuron][_mode] = _lus

    def get_lus(self,_mode=None,thresh=0):
        """
        Returns the LUS scores above specified threshold
        
        Parameters
        ----------
        _mode : str
          Either 'pre', 'post' or 'gap' synaptic mode. (default = None)
        thresh : float
          Threshold for LUS score. (default = 0)
        """
        return [self.lus[n][_mode] for n in self.lus
                if self.lus[n][_mode] >= thresh]

    def get_data(self):
        """
        Returns model results. 
        Row format: cell_name,pre_lus,gap_lus       
        """
        data = [[n,self.lus[n]['pre'],self.lus[n]['gap']] for n in self.lus]
        return data
        

def wbe(exp,C,cells = None):
    """
    Whole-cell binary expression model (WBE)

    Parameters
    ----------
    exp : Expression object
    C : Connectome object (Networkx datatype)
    cells: list
     List of cells. By default all cells will be used.
    """

    if not cells: cell = C.A.nodes()
    WBE = Model('WBE')
    for n in cells:
        _num,_den = 0.0,0.0
        if not C.C.has_node(n): continue
        neigh = set(C.A.neighbors(n))

        if C.C.has_node(n):
            syn = set(C.C.neighbors(n))
            lus = _wbe(exp,syn,neigh)
            #print(n,lus,sorted(syn),sorted(set(neigh)))
            WBE.add_lus(n,lus,'pre')

        if C.E.has_node(n):
            syn = set(C.E.neighbors(n))
            lus = _wbe(exp,syn,neigh)
            #print(n,lus,sorted(syn),sorted(set(neigh)))
            WBE.add_lus(n,lus,'gap')
        
    return WBE

def _wbe(exp,syn,neigh):
    """
    Computes LUS for WBE
    
    Parameters
    ----------
    exp : Expression object
    syn : list
      Cells that synapse
    neigh : list
      Cells that are adjacent but not synaptic
    """
    nonsyn = neigh - syn
    lus = -1
    _num,_den = 0.0,0.0
    for s in syn:
        _den += 1
        for ns in nonsyn:
            diff = exp.compute_difference(s,ns)
            if diff < 1:
                _num += 1
                break
    if _den > 0: lus = 1 - (_num/_den)
    return lus
    

def ie(exp,C,cells=None,iters=1000,mode='post'):
    """
    Isoform expression model (IE)

    Parameters
    ----------
    exp : Expression object
    C : Connectome object (Networkx datatype)
    iters : int
      Number of iterations to simulate alt. splicing. (default=1000)
    mode : str
      'pre' or 'post' synaptic expression. 
    """
    exp.load_genes(splice=True)
    lus = np.zeros((iters,2))
    for i in tqdm(range(iters),desc='IE'):
        exp.assign_expression()
        exp.binarize()
        exp.difference_matrix()
        tmp = wbe(exp,C,cells=cells)
        lus[i,0] = np.mean(tmp.get_lus('pre'))
        lus[i,1] = np.mean(tmp.get_lus('gap'))
    return lus

def sbe(exp,syn_chem,syn_gap):
    """
    Subcellular binary expression model

    Parameters
    ----------
    exp : Expression object
    syn_chem : dict
      Chemical synapse data
    syn_gap : dict
      Gap junction data
    """
    lus_pre = _sbe(exp,syn_chem)
    lus_gap = _sbe(exp,syn_gap)
    
    SBE = Model('SBE')
    for n in lus_pre: SBE.add_lus(n,lus_pre[n],'pre')
    for n in lus_gap: SBE.add_lus(n,lus_gap[n],'gap')
    return SBE

def _sbe(exp,synapses):
    """ 
    Subcellular binary expression model for chemical synapses

    Paramters:
    ----------
    exp : Expression object
    synapses: dict
        Synapse data
    """
    lus = {}
    for cell in synapses:
        lus[cell] = [0.,0.]
        for cont in synapses[cell]:
            adj = set(synapses[cell][cont]['neighbors'])
            post = set(synapses[cell][cont]['partners'])
            nonsyn = adj - post
            for s in post:
                lus[cell][0] += 1
                for ns in nonsyn:
                    if ns not in exp.cells.keys(): continue
                    diff = exp.compute_difference(s,ns)
                    if diff < 1:
                        lus[cell][1] += 1
                        break
    LUS = dict([(n,1 - lus[n][1]/lus[n][0]) for n in lus if lus[n][0] > 0])
    return LUS 

