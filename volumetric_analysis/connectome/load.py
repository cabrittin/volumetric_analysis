"""
conectome.load.py

Module for routine loading of connectome data

Functions:
---------
from_db(db,chemical=False,electrical=False,adjacency=False,add_poly=False,
            touch_density=False,remove=None,group=None,dataType='igraph'):
   Loads connectome data from database. Returns Connectome object.

scrube_neurons(_neurons)
   Removes neurons in the SCREEN (global variable) list.

created: Christopher Brittin
date: 01 November 2018

"""

#Brittin modules
import db 
from connectome.connectome_igraph import Connectome
from connectome.connectome_networkx import Connectome as nxConnectome


SCREEN = ['old','duplicate','Frag','error','unk']

def from_db(_db,chemical=False,electrical=False,adjacency=False,add_poly=False,
            touch_density=False,remove=None,group=None,dataType='igraph'):
    """
    Loads connectome data from database. Returns Connectome object.
    
    Parameters:
    -----------
    _db : str
      database name
    chemical : bool (default False)
      If true, load chemical graph 
    electrical : bool (default False)
      If true, load gap junction graph.
    adjacency : bool (default False)
      If true, load adjacency graph.
    add_poly : bool (default False)
      If true, record number of polyad and monad syanpses.
    touch_density : bool (default (False)
      If true, compute the touch density between adjacent cell in the 
      adjacency graph
    remove : list (default None)
      Remove vertices/cells from the connectome object in the remove list
    group : dict (default None)
      Dictionary (key,val) = (cellName,groupName) to groups cells
    dataType : str (default 'igraph')
      Either use the 'igraph' or 'networkx' graph data structures. 


    Returns:
    --------
    Connectome object
    
    """
    
    if _db == 'N2U':
        end = 325
    else:
        end = 500

    con = db.connect.default(_db)
    cur = con.cursor()
    
    if adjacency:
        neurons = sorted(db.mine.get_adjacency_cells(cur))
        adjacency = db.mine.get_adjacency_data(cur)
        if dataType == 'igraph':
            C = Connectome(_db,neurons)
        elif dataType == 'networkx':
            C = nxConnectome(_db,neurons)
        C.load_adjacency(adjacency,directed=False)
    elif touch_density:
        neurons = sorted(_db.mine.get_adjacency_cells(cur))
        adjacency = db.mine.get_touch_density(cur)
        if dataType == 'igraph':
            C = Connectome(_db,neurons)
        C.load_adjacency(adjacency,directed=True)
    else:
        neurons = sorted(scrub_neurons(db.mine.get_neurons(cur)))
        if dataType == 'igraph':
            C = Connectome(_db,neurons)
        elif dataType == 'networkx':
            C = nxConnectome(_db,neurons)
        #C = Connectome(db,neurons)

    if chemical:
        synapses = db.mine.get_synapse_data(cur,'chemical',end=end)
        C.load_chemical(synapses,add_poly=add_poly)

    if electrical:
        synapses = db.mine.get_synapse_data(cur,'electrical',end=end)
        C.load_electrical(synapses)

    if remove: C.remove_cells(remove)
    if group: C.group_cells(group['map'],key=group['key'])
        
    return C

def scrub_neurons(_neurons):
    """
    Removes neurons in the SCREEN (global variable) list.

    Parameters:
    ----------
    _neurons : list
     List of cell names to screen

    Returns:
    --------
    List of screened cells
    
    """
    neurons = []
    for n in _neurons:
        remove = False 
        for s in SCREEN:
            if s in n:
                remove = True
                break
        if not remove: neurons.append(n)
    return neurons


