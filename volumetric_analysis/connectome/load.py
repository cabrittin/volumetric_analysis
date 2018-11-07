"""
Functions for loading connectome

"""

#Brittin modules
import DB
from Connectome.connectome_igraph import Connectome
from Connectome.connectome_networkx import Connectome as nxConnectome
import aux
import LUT


SCREEN = ['old','duplicate','Frag','error','unk']

def from_db(db,chemical=False,electrical=False,adjacency=False,add_poly=False,
            touch_density=False,remove=None,group=None,dataType='igraph'):
    if db == 'N2U':
        end = 325
    else:
        end = 500

    con = DB.connect.default(db)
    cur = con.cursor()
    
    if adjacency:
        neurons = sorted(DB.mine.get_adjacency_cells(cur))
        adjacency = DB.mine.get_adjacency_data(cur)
        if dataType == 'igraph':
            C = Connectome(db,neurons)
        elif dataType == 'networkx':
            C = nxConnectome(db,neurons)
        C.load_adjacency(adjacency,directed=False)
    elif touch_density:
        neurons = sorted(DB.mine.get_adjacency_cells(cur))
        adjacency = DB.mine.get_touch_density(cur)
        if dataType == 'igraph':
            C = Connectome(db,neurons)
        C.load_adjacency(adjacency,directed=True)
    else:
        neurons = sorted(scrub_neurons(DB.mine.get_neurons(cur)))
        if dataType == 'igraph':
            C = Connectome(db,neurons)
        elif dataType == 'networkx':
            C = nxConnectome(db,neurons)
        #C = Connectome(db,neurons)

    if chemical:
        synapses = DB.mine.get_synapse_data(cur,'chemical',end=end)
        C.load_chemical(synapses,add_poly=add_poly)

    if electrical:
        synapses = DB.mine.get_synapse_data(cur,'electrical',end=end)
        C.load_electrical(synapses)

    if remove: C.remove_cells(remove)
    if group: C.group_cells(group['map'],key=group['key'])
        
    return C

def scrub_neurons(_neurons):
    neurons = []
    for n in _neurons:
        remove = False 
        for s in SCREEN:
            if s in n:
                remove = True
                break
        if not remove: neurons.append(n)
    return neurons


