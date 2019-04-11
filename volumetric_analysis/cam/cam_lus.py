"""
cam_lus.py

Submodule for testing different CAM expression models
via the LUS.

Author: Christopher Brittin
Created: 08 February 2018

"""

import progressbar
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
            WBE.add_lus(n,lus,'pre')

        if C.E.has_node(n):
            syn = set(C.E.neighbors(n))
            lus = _wbe(exp,syn,neigh)
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
    

def ie(exp,C,iters=1000,mode='post'):
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
    bar = progressbar.ProgressBar(maxval=iters-1,
                                  widgets=[progressbar.Bar('.','[',']'),
                                           'IE: Alt. splice',
                                           progressbar.Percentage()])
    bar.start()
    lus = np.zeros((iters,2))
    for i in range(iters):
        exp.assign_expression_patterns(mode=mode)
        tmp = wbe(exp,C)
        lus[i,0] = np.mean(tmp.get_lus('pre'))
        lus[i,1] = np.mean(tmp.get_lus('gap'))
        bar.update(i)
    return lus

def sbe(exp,end=500,rmchem=[]):
    """
    Subcellular binary expression model

    Parameters
    ----------
    exp : Expression object
    end : int
       Most posterior section number
    """
    lus_pre = sbe_pre(exp,end=end,rmchem=rmchem)
    lus_gap = sbe_gap(exp,end=end)

    SBE = Model('SBE')
    for n in lus_pre: SBE.add_lus(n,lus_pre[n],'pre')
    for n in lus_gap: SBE.add_lus(n,lus_gap[n],'gap')
    return SBE

def sbe_pre(exp,end=500,rmchem=None):
    """
    Subcellular binary expression model for chemical synapses

    Parameters
    ----------
    exp : Expression object
    end : int
       Most posterior section number
    """   
    syn = get_synapses(exp.cur,exp.nodes,end=end,rmsyn=rmchem)
    idx = 0
    lus = {}
    objs = syn.keys()
    N = len(objs)
    bar = progressbar.ProgressBar(maxval=N-1,
                                  widgets=[progressbar.Bar('.','[',']'),
                                           'SBE: Synapse object',
                                           progressbar.Percentage()])
    bar.start()
    #objs = ['71193']
    for o in objs:
        if not syn[o].pre_obj: continue
        adj = set(db.mine.get_object_adjacency(exp.cur,syn[o].pre_obj))
        post = set(syn[o].post)
        nonsyn = adj - post
        if syn[o].pre not in lus: lus[syn[o].pre] = [0.,0.]
        for s in post:
            lus[syn[o].pre][0] += 1
            for ns in nonsyn:
                if ns not in exp.cells.keys(): continue
                diff = exp.compute_difference(s,ns)
                if diff < 1:
                    lus[syn[o].pre][1] += 1
                    break

        bar.update(idx)
        idx += 1
        
    LUS = dict([(n,1 - lus[n][1]/lus[n][0]) for n in lus if lus[n][0] > 0])
    #print(1-np.mean(LUS),np.std(LUS)/np.sqrt(len(LUS)))
    return LUS

def sbe_gap(exp,end=500):
    """
    Subcellular binary expression model for gap junctions

    Parameters
    ----------
    exp : Expression object
    end : int
       Most posterior section number
    """ 
    syn = get_gapjunctions(exp.cur,exp.nodes,end=end)
    idx = 0
    objs = syn.keys()
    N = len(objs)
    lus = {}
    bar = progressbar.ProgressBar(maxval=N-1,
                                  widgets=[progressbar.Bar('.','[',']'),
                                           'SBE: Synapse gap object',
                                           progressbar.Percentage()])
    bar.start()    
    for o in objs:
        cells = syn[o].get_cells()
        for c in cells:
            try:
                adj = set(db.mine.get_object_adjacency(exp.cur,syn[o].cellObjs[c]))
                s = syn[o].partner[c]
                post = set([s])
                nonsyn = adj - post
                if c not in lus: lus[c] = [0.,0.]
                lus[c][0] += 1
                for ns in nonsyn:
                    diff = exp.compute_difference(s,ns)
                    if diff < 1:
                        lus[c][1] += 1
                        break
                bar.update(idx)
            except:
                print('Mysql: %s' %c)
        idx += 1
    LUS = dict([(n,1 - lus[n][1]/lus[n][0]) for n in lus if lus[n][0] > 0])
    return LUS

class SynObj:
    """
    Class for chemical synapse objects. Used to support SBE analysis
    
    Attributes
    ----------
    object : int
      Elegance object number
    pre : str
      Presynaptci cell name
    pre_obje : int
      Presynaptic object number
    post : dict 
      Dictionary that maps postsynaptic cells to object number
    relabel : 
       Used to for cell name compatability issues

    Methods
    -------
    add_pre(cell,obj)
      Adds the presynaptic cell
    add_post(cell,obj)
      Adds the postsynaptic cell
    """
    def __init__(self,obj):
        """
        Attributes
        ----------
        obj : int
          Elegance object number
        """
        self.object = obj
        self.pre = None
        self.pre_obj = None
        self.post = {}
        self.relabel = {'PVR.':'PVR',
                        'PHAR.':'PHAR'}

    def add_pre(self,cell,obj):
        """
        Adds the presynaptic cell
        
        Parameters
        ----------
        cell : str
           cell name
        obj : int
           Elegance object number
        """
        if cell in self.relabel: cell = self.relabel[cell]
        self.pre = cell
        self.pre_obj = obj

    def add_post(self,cell,obj):
        """
        Adds the postsynaptic cell
        
        Attributes
        ----------
        cell : str
           cell name
        obj : int
           Elegance object number     
        
        """
        if cell in self.relabel: cell = self.relabel[cell]
        if cell: self.post[cell] = obj


    
def get_synapses(cur,nodes,end=500,rmsyn=[]):
    """
    Gets synapse data from the Elegance database
    """
    sql = ("select mid,pre,post1,post2,post3,post4,"
           "preobj,postobj1,postobj2,postobj3,postobj4 "
           "from synapsecombined "
           "join object "
           "on synapsecombined.mid = object.OBJ_Name "
           "join image "
           "on object.IMG_Number = image.IMG_Number "
           "where synapsecombined.type = 'chemical' "
           "and image.IMG_SectionNumber <= %d "
           %end)
    cur.execute(sql)
    syn = {}
    for a in cur.fetchall():
        [obj,pre,p1,p2,p3,p4,preo,p1o,p2o,p3o,p4o] = a
        s = SynObj(obj)
        if pre not in nodes: continue
        s.add_pre(pre,preo)
        if p1 in nodes and (pre,p1) not in rmsyn: s.add_post(p1,p1o)
        if p2 in nodes and (pre,p2) not in rmsyn: s.add_post(p2,p2o)
        if p3 in nodes and (pre,p3) not in rmsyn: s.add_post(p3,p3o)
        if p4 in nodes and (pre,p4) not in rmsyn: s.add_post(p4,p4o)
        if s.pre and s.post: syn[obj] = s

    
    return syn


class GapObj:
    """
    Class for gap junction objects. Used to support SBE analysis.
    
    Attributes
    ----------
    object : int
      Elegance object number
    cellObjs : dict
      Dictionary that maps cell name to object number
    partners : dict
      Dictionary that maps gap junction partners.
    
    Methods
    -------
    add_cell_obj(cell,obj)
      Adds cell and object number
    add_partner(cell1,cell2)
      Maps cell1 to cell2
    get_cells()
      Returns list of cells
    
    """
    def __init__(self,obj):
        """
        Attributes
        ----------
        obj : int
          Elegance object number       
        """
        self.object = obj
        self.cellObjs = {}
        self.relabel = {'PVR.' : 'PVR',
                        'PHAR.': 'PHAR'}
        self.partner = {}
        
    def add_cell_obj(self,cell,obj):
        """
        Adds cell and object number
        
        Attributes
        ----------
        cell : str
          Cell name
        obj : int
          Elegance object number
        
        """
        if cell in self.relabel: cell = self.relabel[cell]
        self.cellObjs[cell] = obj
            
    def add_partners(self,cell1,cell2):
        """
        Maps cell1 to cell2
        
        Attributes
        ----------
        cell1 : str
          Cell name
        cell2 : str
          Cell name
        """
        if cell1 in self.relabel: cell1 = self.relabel[cell1]
        if cell2 in self.relabel: cell2 = self.relabel[cell2]
        self.partner[cell1] = cell2
        self.partner[cell2] = cell1

    def get_cells(self):
        """
        Returns list of cell names
        """
        return self.cellObjs.keys()

def get_gapjunctions(cur,nodes,end=500):
    """
    Gets gap junctions from the Elegance database
    """
    sql = ("select mid,pre,post1,preobj,postobj1 "
           "from synapsecombined "
           "join object "
           "on synapsecombined.mid = object.OBJ_Name "
           "join image "
           "on object.IMG_Number = image.IMG_Number "
           "where synapsecombined.type = 'electrical' "
           "and image.IMG_SectionNumber <= %d "
           %end)
    cur.execute(sql)
    syn = {}
    for a in cur.fetchall():
        [obj,pre,p1,preo,p1o] = a
        g = GapObj(obj)
        if pre in nodes and p1 in nodes:
            g.add_partners(pre,p1)
            g.add_cell_obj(pre,preo)
            g.add_cell_obj(p1,p1o)
            syn[obj] = g
    return syn
        
