"""
load_adjacency.py

Load adjacency data produced by parsetrakem2 module to elegance database.  

Parameters
----------
adj : str
  xml file with adjacency data
segments : str
  xml file with segment measures

created: Christopher Brittin
date: 01 November 2018
"""

import argparse
from lxml import etree
import numpy as np

import db

def get_closest_object(refxy,testxy):
    mindist = 1e9
    minobj = None
    for t in testxy:
        d = distance(refxy,(t[1],t[2]))
        if d < mindist:
            mindist = d
            minobj = t[0]
    return minobj
    
def distance(c1,c2):
    return np.sqrt((c1[0]-c2[0])**2 + (c1[1]-c2[1])**2)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('adj',
                        action="store",
                        help="XML file with adjacency data")

    parser.add_argument('seg',
                        action = 'store',
                        help = "XML file with segmentation measures")

    parser.add_argument('db',
                        action = 'store',
                        help = 'Database name')
                        
    
    params = parser.parse_args()

    con = db.connect.default(params.db)
    cur = con.cursor()
    
    tree = etree.parse(params.seg)
    root = tree.getroot()
    layers = sorted([l.get('name') for l in root.findall('layer')])
    segments = {}
    for _l in layers:
        print(_l)
        if params.db == 'N2U':
            if 'VC' in _l:
                newl = _l.replace('_VC_','VC')
            else:
                newl = _l.replace('_','NR')
        elif params.db == 'JSH':
            newl = _l.replace('JSH','JSHJSH')
        segments[_l] = {}
        l = root.find("layer[@name='%s']" %_l)
        segs = l.findall('segment')
        for s in segs:
            name = s.find('name').text
            index = s.find('index').text
            centx = int(s.find('centx').text)
            centy = int(s.find('centy').text)
            if name not in segments[_l]:
                segments[_l][name] = {}
            #segments[_l][name][index] = [centx,centy]
            xy = db.mine.get_object_xy_in_layer(cur,name,newl)
            minobj = get_closest_object((centx,centy),xy)
            segments[_l][name][index] = minobj

    tree = etree.parse(params.adj)
    root = tree.getroot()
    
    layers = sorted([l.get('name') for l in root.findall('layer')])
    data = []
    for _l in layers:
        l = root.find("layer[@name='%s']" %_l)
        areas = l.findall('area')
        for a in areas:
            c1 = a.find('cell1').text
            c2 = a.find('cell2').text
            i1 = a.find('index1').text
            i2 = a.find('index2').text
            adj = a.find('adjacency').text
            o1 = segments[_l][c1][i1]
            o2 = segments[_l][c2][i2]
            data.append([c1,c2,i1,i2,_l,adj,o1,o2])

    db.insert.adjacency(cur,data)
    con.commit()
    con.close()

            
            
