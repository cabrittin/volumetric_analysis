"""
cam_contact.py

Look for CAM expression that is over represented in larger contact.

@author Christopher Brittin 
@date 28 April 2019
"""

import sys
sys.path.append('.')
sys.path.append('./volumetric_analysis')
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import networkx as nx
from itertools import combinations
from collections import defaultdict
from scipy.spatial.distance import hamming

from mat_loader import MatLoader
from cam.expression import Matrix
import cam.cam_predict as predict
from connectome.load import from_db

def contact_profile(e,A,cell,alpha=0.33):
    neigh = sorted(A[cell].items(), key=lambda edge: edge[1]['weight'])
    k = len(neigh)
    low = [e.cells[n[0]] for n in neigh[:int(k*alpha)]]
    high = [e.cells[n[0]] for n in neigh[-int(k*alpha):]]
    return low,high
    
def cam_corr(M,cell,test):
    cor = []
    for i in test:
        #r = pearsonr(M[cell,:],M[i,:])[0]
        r = hamming(M[cell,:],M[i,:])
        if np.isnan(r): 
            continue
        else:
            cor.append(r)
    return cor

def overlap(u,v):
    w = u + v
    #m = np.sum(u) + np.sum(v)
    m = float(len(np.where(w > 0)[0]))
    n = len(np.where(w > 1)[0])
    if m > 0: return 1 - n / m
    return np.nan

DB = ['JSH','N2U']
REMOVE = ['VB01', 'VD01']
cmap = {'SMN':0.0,'HMNa':0.1, 'HMNp':0.2,'I1':0.4, 'I2':0.5, 
        'Sa':0.7, 'Sp1':0.9, 'Sp2':0.9, 'XMN':1.0}

idx_gene = range(85,98)
#idx_gene = range(54,85)
#idx_gene = range(54)
#idx_gene = range(98,106)
#idx_gene = range(106)

METHOD = 1

if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('matrix',
                        action = 'store',
                        help = 'Path to matrix file')
    
    parser.add_argument('db',
                        action = 'store',
                        help = 'Database')

    #parser.add_argument('cell',action='store',help='Cell name')

    params = parser.parse_args()
    #cell = params.cell

    ML = MatLoader()
    ML.load_lrmap()
    nodes = sorted(ML.load_reduced_nodes())


    C = from_db(params.db,adjacency=True,dataType='networkx',remove=REMOVE)
    
    e = Matrix(ML.cam,params.matrix)
    e.load_genes()
    e.load_cells(nodes)
    e.assign_expression()
    e.binarize()
    #np.random.shuffle(e.M)
    
    M = e.E[:,idx_gene]
  
    #np.random.shuffle(M)

    if METHOD == 1:
        data = defaultdict(list)
        comb = combinations(nodes,2)
        for (u,v) in comb:
            l = nx.shortest_path_length(C.A,source=u,target=v)
            i = e.cells[u]
            j = e.cells[v]
            #r = pearsonr(M[i,:],M[j,:])[0]
            #r = hamming(M[i,:],M[j,:])
            r = overlap(M[i,:],M[j,:])
            if np.isnan(r): continue
            data[l].append(r)

        for (l,d) in sorted(data.items()):
            mu = np.mean(d)
            std = np.std(d)
            print(l,mu,std,len(d))
    
    if METHOD == 2:
        data = np.zeros((C.A.number_of_edges(),2))
        idx = 0
        dlow,dhigh = [],[]
        for v in C.A.nodes():
            idx = e.cells[v]
            low,high = contact_profile(e,C.A,v,alpha=0.5)
            clow = cam_corr(M,idx,low)
            chigh = cam_corr(M,idx,high)
            if clow: dlow.append(np.mean(clow))
            if chigh: dhigh.append(np.mean(chigh))
        #print(dlow)
        print(np.mean(dlow),np.std(dlow))
        print(np.mean(dhigh),np.std(dhigh))
        #plt.figure()
        #plt.plot([dlow,dhigh])
        #plt.show()

    if METHOD == 3:
        LOW,HIGH = [],[]
        E = M[:,:]
        E[E>0] = 1
        print(E.shape)
        idx = 0
        no_cad = []
        for i in range(M.shape[0]):
            if np.sum(M[i,:]) > 0 : 
                idx += 1
            else: 
                print(nodes[i])
                no_cad.append(nodes[i])
        print(idx)
        nodes = sorted(list(set(nodes) - set(no_cad)))
        k = len(nodes)
        CV = -9 * np.ones((k,k))
        ndict = {nodes[i]:i for i in range(k)}
        comb = combinations(nodes,2)
        for (u,v) in comb:
            i = e.cells[u]
            j = e.cells[v]
            #r = pearsonr(M[i,:],M[j,:])[0]
            #r = hamming(M[i,:],M[j,:])
            r = overlap(M[i,:],M[j,:])
            if np.isnan(r): continue
            CV[ndict[u],ndict[v]] = r
            CV[ndict[v],ndict[u]] = r

        dx = 10
        for i in range(CV.shape[0]):
            rank = np.argsort(CV[i,:])
            low = rank[1:1+dx]
            high = rank[-dx:]
            cell = nodes[i]
            for j in low:
                lcell = nodes[j]
                rcell = ML.lrmap[lcell]
                w1,w2 = 0,0
                if C.A.has_edge(cell,lcell):
                    w1 = C.A[cell][lcell]['weight']
                if C.A.has_edge(cell,rcell):
                    w2 = C.A[cell][rcell]['weight']
                w = max(w1,w2)
                LOW.append(w)
            
            for j in high:
                lcell = nodes[j]
                rcell = ML.lrmap[lcell]
                w1,w2 = 0,0
                if C.A.has_edge(cell,lcell):
                    w1 = C.A[cell][lcell]['weight']
                if C.A.has_edge(cell,rcell):
                    w2 = C.A[cell][rcell]['weight']
                w = max(w1,w2)
                HIGH.append(w)

        print(np.mean(LOW),np.std(LOW))
        print(np.mean(HIGH),np.std(HIGH))
        
        maxw = max(max(LOW),max(HIGH))
        plt.figure()
        plt.hist(HIGH,bins=1000,range=(0,maxw),histtype='step',density=True,cumulative=True,color='b',label='Large distnance')
        plt.hist(LOW,bins=1000,range=(0,maxw),histtype='step',density=True,cumulative=True,color='r',label='Small distance')
        plt.xlabel('Surface area contact')
        plt.ylabel('Cumulative distribution')
        plt.xlim([0,70000])
        plt.ylim([0,1])
        plt.legend()
        plt.show()


        

 



    

