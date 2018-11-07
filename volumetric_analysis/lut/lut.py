"""
lut.py
Author: Christopher Brittin

Functions for looking up data from files.

"""
import sys
sys.path.append(r'..')
import os
from os.path import expanduser
import aux.read

def db():
    return expanduser("~") + '/.dbinfo.txt'

def skeleton(neuron,db):
    pre = get_lut_pre()    
    skel = pre + '/lut/' + db + '/skeletons.txt'
    skel = aux.read.into_list(skel)
    return matchFile(neuron,skel)
    
def wavefront(neuron,db):
    pre = get_lut_pre()
    obj = pre + '/lut/' + db + '/objFiles.txt'
    obj = aux.read.into_list(obj)
    return matchFile(neuron,obj)

def links(neuron,db):
    pre = get_lut_pre()
    links = pre + '/lut/' + db + '/links.txt'
    links = aux.read.into_list(links)
    return matchFile(neuron,links)

def endpoints(neuron,db):
    pre = get_lut_pre()
    epts = pre + '/lut/' + db + '/endpoints.txt'
    epts = aux.read.into_list2(epts)
    for e in epts:
        if neuron == e[0]:
            ee = [int(p) for p in e[1:] if p]
            return ee

def neuron_class(db):
    #Added neuron class
    pre = get_lut_pre()
    nclass = pre + '/lut/' + db + '/nclass.txt'
    nclass = aux.read.into_dict(nclass)
    return nclass

def loader(db):
    pre = get_lut_pre()
    loader = pre + '/lut/' + db + '/loader.txt'
    return  aux.read.into_dict(loader)
    

def get_lut_pre():
    if os.path.exists('lut'):
        pre = '.'
    else:
        pre = '..'
    return pre

def matchFile(mtch,lst):
    mfile = None
    for l in lst:
        if _matchFile(mtch,l):
            mfile = l
            break
    return mfile

def _matchFile(neuron,fname):
    fname = fname.split('.')
    fname = fname[0].split('/')
    if neuron == fname[-1]:
        return True
    else:
        return False
    
