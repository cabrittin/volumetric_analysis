"""
optimized_expression.py

Finds the optimal expression pattern that maximizes coverage of synaptic partners
while minimizing coverage of neighbors. 

@author Christopher Brittin
@date 2019 April 14

"""
import sys
sys.path.append('./volumetric_analysis')
sys.path.append('.')
import argparse
import numpy as np
from random import sample
from tqdm import tqdm
from random import random
import matplotlib.pyplot as plt

from mat_loader import MatLoader
from cam.expression import Matrix

class Optimize(Matrix):
    def __init__(self,fgenes,fexp):
        Matrix.__init__(self,fgenes,fexp)

    def reduced_expression(self,cells):
        idx = [self.cells[c] for c in cells]
        E = self.E[idx,:]
        jdx = np.where(E.any(axis=0))[0]
        E = E[:,jdx]
        genes = [self.gene_idx[j] for j in jdx]
        return genes,E.T

def build_synapse_matices(cells,synapses):
    n = len(cells)
    m = len(synapses)
    cells_idx = {cells[i]:i for i in range(n)}
    
    P = np.zeros((n,m))
    A = np.zeros((n,m))

    jdx = 0
    for cont in synapses:
        partners = set(synapses[cont]['partners'])
        neighbors = set(synapses[cont]['neighbors'])
        nonsyn = neighbors - partners
        for c in partners: P[cells_idx[c],jdx] = 1
        for c in nonsyn: A[cells_idx[c],jdx] = 1
        jdx += 1
    
    return P,A


def anneal(sol):
    old_cost = cost(sol)
    T = 1.0
    T_min = 0.00001
    alpha = 0.9
    while T > T_min:
        i = 1
        while i <= 100:
            new_sol = neighbor(sol)
            new_cost = cost(new_sol)
            ap = acceptance_probability(old_cost, new_cost, T)
            if ap > random():
                sol = new_sol
                old_cost = new_cost
            i += 1
        T = T*alpha
    return sol, cost

def cost(E,S,A,gamma = [0.5,0.4,0.1]):
    syn = np.sum(np.dot(E,S))
    adj = np.sum(np.dot(E,A))
    sparse = np.sum(E != 0)
    J = gamma[0] * syn - gamma[1] * adj - gamma[2] * sparse
    return J
    
def acceptance_probability(new_cost,old_cost,T):
    p =  np.exp((new_cost - old_cost) / T)
    return min(p,1)

def perturb(coord,E):
    rdx = np.random.randint(0,len(coord))
    (a,b) = coord[rdx]
    E[a,b] = int(not E[a,b])
    return E

def run_optimization(E0,S,A,iters=100):
    (n,m) = E0.shape
    E = np.zeros(E0.shape)
    coord = np.nonzero(E0)
    coord = list(zip(coord[0],coord[1]))
    k = int(len(coord) / 2)
    for (a,b) in sample(coord,k): E[a,b] = 1
    if not iters: iters = 10*len(coord)
    
    old_cost = cost(E,S,A)
    T = 1.0
    T_min = 0.00001
    alpha = 0.9
    cost_rec = [old_cost]
    while T > T_min:
        i = 1
        while i <= iters:
            #Perturb E
            rdx = np.random.randint(0,len(coord))
            (a,b) = coord[rdx]
            E[a,b] = int(not E[a,b])
            #rdx = np.random.randint(0,n)
            #old_row = E[rdx,:]
            #nonzero = np.nonzero(old_row)
            #E[rdx,nonzero] = 0

            #Compute new cost 
            new_cost = cost(E,S,A)
            #Decide to accept perturbation
            ap = acceptance_probability(new_cost,old_cost,T)
            if ap > random():
                old_cost = new_cost
                cost_rec.append(new_cost)
            else:
                E[a,b] = int(not(E[a,b]))
                #E[rdx,nonzero] = 1
            i += 1
        T *= alpha
    return E,cost_rec

cam = 'mat/cam_isoforms.txt'

if __name__=='__main__':

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('matrix',
                        action = 'store',
                        help = 'Path to matrix file')

    parser.add_argument('deg',
                        action = 'store',
                        type= int,
                        help = 'Conserved degree')

    
    params = parser.parse_args()

    M = MatLoader()
    M.load_left()
    C = M.load_consensus_graphs(params.deg)
    S = M.load_consensus_chemical_synapse(params.deg)


    e = Optimize(cam,params.matrix)
    e.load_genes()
    e.load_cells(sorted(C.A.nodes()))
    e.assign_expression()
    e.binarize()
    
    cell = 'AIYL'
    neigh = sorted(C.A.neighbors(cell))
    genes,E = e.reduced_expression(neigh) 
    
    P,A = build_synapse_matices(neigh,S[cell])

    Eopt,cost_rec = run_optimization(E,P,A)


    print(np.sum(E),np.sum(Eopt))
    print(np.sum(np.dot(E,P)),np.sum(np.dot(E,A)))
    print(np.sum(np.dot(Eopt,P)),np.sum(np.dot(Eopt,A)))
    idx = set(np.nonzero(Eopt)[0])
    _genes = [genes[i] for i in idx]
    print(len(genes),len(_genes))
    print(sorted(_genes))
    plt.figure()
    plt.plot(cost_rec,'b-')
    plt.show()
