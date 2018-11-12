"""
cam_lus.py

Submodule for testing different CAM expression models
via the LUS.

Author: Christopher Brittin
Created: 08 February 2018

"""

import progressbar
import numpy as np

import DB

class Model(object):
    def __init__(self,_modelName):
        self.lus = {}
        self.modelName = _modelName

    def add_lus(self,_neuron,_lus,_mode):
        if _neuron not in self.lus:
            self.lus[_neuron] = {'gap':-1,'pre':-1,'post':-1}
        self.lus[_neuron][_mode] = _lus

    def get_lus(self,_mode=None,thresh=0):
        return [self.lus[n][_mode] for n in self.lus
                if self.lus[n][_mode] >= thresh]

    def get_data(self):
        data = [[n,self.lus[n]['pre'],self.lus[n]['gap']] for n in self.lus]
        return data
        

def wbe(exp,C):
    WBE = Model('WBE')
    for n in C.A.nodes():
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
    
    bar = progressbar.ProgressBar(maxval=iters-1,
                                  widgets=[progressbar.Bar('.','[',']'),
                                           'IE: Alt. splice',
                                           progressbar.Percentage()])
    bar.start()
    lus = np.zeros((iters,2))
    for i in xrange(iters):
        exp.assign_expression_patterns(mode=mode)
        tmp = wbe(exp,C)
        lus[i,0] = np.mean(tmp.get_lus('pre'))
        lus[i,1] = np.mean(tmp.get_lus('gap'))
        bar.update(i)
    return lus

def sbe(exp,end=500):
    lus_pre = sbe_pre(exp,end=500)
    lus_gap = sbe_gap(exp,end=500)

    SBE = Model('SBE')
    for n in lus_pre: SBE.add_lus(n,lus_pre[n],'pre')
    for n in lus_gap: SBE.add_lus(n,lus_gap[n],'gap')
    return SBE

def sbe_pre(exp,end=500):
    syn = get_synapses(exp.cur,exp.nodes,end=end)
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
        adj = set(DB.mine.get_object_adjacency(exp.cur,syn[o].pre_obj))
        post = set(syn[o].post)
        nonsyn = adj - post
        if syn[o].pre not in lus: lus[syn[o].pre] = [0.,0.]
        for s in post:
            lus[syn[o].pre][0] += 1
            for ns in nonsyn:
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
                adj = set(DB.mine.get_object_adjacency(exp.cur,syn[o].cellObjs[c]))
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
    def __init__(self,obj):
        self.object = obj
        self.pre = None
        self.pre_obj = None
        self.post = {}
        self.relabel = {'PVR.':'PVR',
                        'PHAR.':'PHAR'}

    def add_pre(self,cell,obj):
        if cell in self.relabel: cell = self.relabel[cell]
        self.pre = cell
        self.pre_obj = obj

    def add_post(self,cell,obj):
        if cell in self.relabel: cell = self.relabel[cell]
        if cell: self.post[cell] = obj

    def _post(self):
        for p in self.post:
            yield p,self.post[p]

    def get_post_neurons(self):
        return self.post.keys()

    def get_pre_obj(self):
        return self.pre_obj
    
    def get_pre(self):
        return self.pre
    
    def get_post(self):
        return self.post
    
def get_synapses(cur,nodes,end=500):
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
        if pre in nodes: s.add_pre(pre,preo)
        if p1 in nodes: s.add_post(p1,p1o)
        if p2 in nodes: s.add_post(p2,p2o)
        if p3 in nodes: s.add_post(p3,p3o)
        if p4 in nodes: s.add_post(p4,p4o)
        if s.pre and s.post: syn[obj] = s

    
    return syn


class GapObj:
    def __init__(self,obj):
        self.object = obj
        self.cellObjs = {}
        self.relabel = {'PVR.' : 'PVR',
                        'PHAR.': 'PHAR'}
        self.partner = {}
        
    def add_cell_obj(self,cell,obj):
        if cell in self.relabel: cell = self.relabel[cell]
        self.cellObjs[cell] = obj
            
    def add_partners(self,cell1,cell2):
        if cell1 in self.relabel: cell1 = self.relabel[cell1]
        if cell2 in self.relabel: cell2 = self.relabel[cell2]
        self.partner[cell1] = cell2
        self.partner[cell2] = cell1

    def get_cells(self):
        return self.cellObjs.keys()

def get_gapjunctions(cur,nodes,end=500):
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
        
