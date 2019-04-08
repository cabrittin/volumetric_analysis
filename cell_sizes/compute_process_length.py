"""
compute_process_length.py

Computes the process length of cell

@author Christopher Brittin
@date 07 April 2019

"""

import sys
sys.path.append('./volumetric_analysis')
import argparse
import networkx as nx
from itertools import combinations
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import db
from connectome.load import from_db
import aux

TCELL = 'AIBL'

class Skeleton(nx.Graph):
    def __init__(self,contin):
        self.contin = contin
        nx.Graph.__init__(self)
        self.scale = np.array([5,5,18])

    def set_endpoints(self):
        epts = [v for v in self.nodes if self.degree(v) == 1]
        comb = combinations(epts,2)
        maxdist = 0
        maxepts = [0,0]
        for (a,b) in comb:
            try:
                length = nx.shortest_path_length(self,a,b,weight='weight')
                if length > maxdist:
                    maxdist = length
                    maxepts = [a,b]
            except:
                continue

        self.endpts = maxepts
        self.length = maxdist

    def set_distances(self):
        tmp = {}
        for (a,b) in self.edges():
            dx = np.multiply((self.node[a]['loc'] - self.node[b]['loc']),self.scale)
            dist = np.sqrt(np.sum(dx**2))
            self._adj[a][b]['weight'] = dist    
         
    def plot_skeleton(self,ax):
        for (i,j) in self.edges():
            x = [self.node[i]['loc'][0],self.node[j]['loc'][0]]
            y = [self.node[i]['loc'][1],self.node[j]['loc'][1]]
            z = [self.node[i]['loc'][2]*self.scale[2],self.node[j]['loc'][2]*self.scale[2]]
            ax.plot(z,x,y,'b-')

    def plot_max_endpoints(self,ax):
        for i in self.endpts:
            x = self.node[i]['loc'][0]
            y = self.node[i]['loc'][1]
            z = self.node[i]['loc'][2]*self.scale[2]
            ax.plot([z],[x],[y],'ro')
            ax.text(z,x,y,str(i),(1,1,0))
        
    def plot_endpoints(self,ax):
        epts = [v for v in self.nodes if self.degree(v) == 1 and v not in self.endpts]
        for i in epts:
            x = self.node[i]['loc'][0]
            y = self.node[i]['loc'][1]
            z = self.node[i]['loc'][2]*self.scale[2]
            ax.plot([z],[x],[y],'go')
            ax.text(z,x,y,str(i),(1,1,0))
            

if __name__=="__main__":
    
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('db',
                        action="store",
                        help="Database name")
    parser.add_argument('fout',
                        action="store",
                        help="Output file")
    
    params = parser.parse_args()
  
    C = from_db(params.db,adjacency=True)
    
    END = 425
    if params.db == 'N2U': END = 325
    con = db.connect.default(params.db)
    cur = con.cursor()
    cell_data = []
    for cell in sorted(C.neurons):
        contins = db.mine.get_contins(cur,cell)
        for c in contins:
            skel = Skeleton(c)
            edges = db.mine.get_contin_edges(cur,c,end=END)
            skel.add_edges_from(edges)
            for o in skel.nodes():
                xyz = np.array(list(map(int,db.mine.get_object_xyz(cur,o))))
                skel.node[o]['loc'] = xyz
            skel.set_distances()
            skel.set_endpoints()
            if skel.length == 0: continue
            print('Cell: ' + cell + ' Contin: ' + str(c))
            l = int(skel.length)
            print(cell,c,l,skel.endpts[0],skel.endpts[1],0)
            print(cell,c,l,skel.endpts[1],skel.endpts[0],0)
            print(cell,c,l,skel.endpts[0],skel.endpts[1],1)
            print(cell,c,l,skel.endpts[1],skel.endpts[0],1)
            fig = plt.figure(figsize=(10,5))
            ax = fig.add_subplot(111,projection='3d')
            skel.plot_skeleton(ax)
            skel.plot_max_endpoints(ax)
            skel.plot_endpoints(ax)
            ax.set_xlabel('Z')
            ax.set_ylabel('X')
            ax.set_zlabel('Y')
            ax.set_xlim([0,8000])
            ax.set_ylim([0,8000])
            ax.set_zlim([0,8000])
            plt.show()


