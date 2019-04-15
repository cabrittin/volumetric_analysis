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

def run_optimization(E0,S,A,iters=None):
    E = np.zeros(E0.shape)
    coord = np.nonzero(E0)
    coord = list(zip(coord[0],cord[1]))
    k = int(len(coord) / 2)
    for (a,b) in sample(coord,k): E[a,b] = 1
    if not iters: iters = 10*len(coord)
    
    for i in tqdm(range(iters),desc="Annealing"):
        syn = np.sum(np.dot(E,S))
        adj = np.sum(np.dot(E,A))
        sparse = np.sum(E != 0)
        rdx = np.random.randint(0,k)
        (a,b) = coord[rdx]
        E[a,b] = int(not E[a,b])
         

    #print(E.shape,S.shape,A.shape,syn.shape,adj.shape)
    #print(np.linalg.norm(syn,ord=2),np.linalg.norm(adj,ord=2))
    #print(np.sum(syn),np.sum(adj))
    #print(np.sum(E != 0))

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

    run_optimization(E,P,A)

