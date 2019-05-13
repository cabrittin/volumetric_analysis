"""
load_radial_phi_z.py

Load spatial map data (radial,phi,z) to database.   

Parameters
----------
segments : str
  xml file with segment measures
db : str
  database name

created: Christopher Brittin
date: 01 November 2018
"""
import sys
sys.path.append('./volumetric_analysis')

import argparse
from lxml import etree
import numpy as np
from tqdm import tqdm
import multiprocessing as mp

import db
from trakem2.parse import ParseTrakEM2
from aux.format import chunk_list

def upload_radial(_db,P,layers):
    con = db.connect.default(_db)
    cur = con.cursor()
    cur.connection.autocommit(True)
    obj_skipped = []
    for (_l,newl) in tqdm(layers, desc='Layers:'):
        B = P.get_boundaries_in_layer(_l,
                                      area_lists=['Pharynx',
                                                  'Phi_Marker'])
        if not B: continue
        B['Pharynx'][0].set_centroid()
        B['Phi_Marker'][0].set_centroid()
        data = []
        objs = db.mine.get_objects_in_layer(cur,newl)
        for o in objs:
            loc = db.mine.get_object_xyz(cur,o)
            if loc: data.append((o,loc[0],loc[1]))
            else: obj_skipped.append(o)
        #data = db.mine.get_synapse_from_layer(cur,newl)
        #print(_l,len(data))
        if data:
            rad = compute_radial_distance(B['Pharynx'][0].path,data)
            phi = compute_angle(B['Pharynx'][0].cent,
                                B['Phi_Marker'][0].cent,data)
            data = [[d,rad[d],phi[d]] for d in rad]
            db.insert.radial_pharynx(cur,data)
    con.close()

    return obj_skipped
 
def compute_radial_distance(pharynx,data):
    radial = {}
    n = len(pharynx)
    A = np.zeros([n,2])
    for i in range(n):
        A[i,:] = pharynx[i]
    for d in data:
        D = np.ones([n,2])
        D[:,0] = D[:,0]*d[1]
        D[:,1] = D[:,1]*d[2]
        dist = A - D
        dist2 = int(min(np.sqrt(dist[:,0]**2 + dist[:,1]**2)))
        radial[d[0]] = dist2
    return radial

def compute_angle(pharynx,marker,data):
    angle = {}
    pharynx = np.array(pharynx)
    marker = np.array(marker)
    for d in data:
        p = np.array([d[1],d[2]])
        v1 = p - pharynx
        v2 = marker - pharynx
        phi = angle_between(v1,v2)
        direction = angle_direction(marker,p)
        angle[d[0]] = direction*phi
    return angle

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def angle_direction(marker,p):
    if marker[0] - p[0] < 0:
        return -1
    else:
        return 1

def adjust_layer_name(db,_l):
    #Needed for synapses
    newl = _l
    if params.db == 'N2U':
        if 'VC' in _l:
            newl = _l.replace('_VC_','VC')
        else:
            newl = _l.replace('_','NR')
    elif params.db == 'JSH':
        newl = _l.replace('JSH','JSHJSH')
    return newl

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('trakem2',
                        action = 'store',
                        help = "TrakEM2 file")

    parser.add_argument('db',
                        action = 'store',
                        help = 'Database name')
    
    parser.add_argument('-n','--nproc',
                        action = 'store',
                        dest = 'nproc',
                        required = False,
                        default = 1,
                        type = int,
                        help ="Number of cores to run multiprocessing")
                    
    
    params = parser.parse_args()

    print('TrakEM2 file: %s' %params.trakem2)
    P = ParseTrakEM2(params.trakem2)
    P.get_layers()    
    P.get_area_lists()

    layers = sorted(P.layers.keys())
    nlayers = [(l,l) for l in layers] 
    
    #clayers = chunk_list(nlayers,params.nproc)
    
    obj_skipped = upload_radial(params.db,P,nlayers)
    for o in obj_skipped:
        print(o)
    #pool = mp.Pool(processes = params.nproc)
    #results = [pool.apply_async(upload_radial,args=(params.db,P,_layers,)) for _layers in clayers]
    ##adj = [o for p in results for o in p.get()]
    #for p in results: p.get()
    



