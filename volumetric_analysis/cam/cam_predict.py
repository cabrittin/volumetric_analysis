"""
cam_predict.py

Functions for testing CAM predictions

@author Christopher Brittin
@date 19 April 2019
"""
import numpy as np
from random import random
from collections import defaultdict


def gene_differential(E,syn,neigh):
    """
    Returns list of gene (indices) that have higher expression in the 
    synaptic partners compared to (nonsynaptic) neighbors

    Parameters:
    -----------
    E : numpy array
     Expression matrix
    syn : list
      list of synaptic partners at each synapse
    neigh : list
      list of (nonsynaptic) neighbors at each synapse

    Note: syn[i] and neigh[i] correpspond to the ith synapse

    """
    
    (n,m) = E.shape
    k = len(syn)
    syn_count = np.zeros(m)
    neigh_count = np.zeros(m)
    for i in range(k):
        sdx = syn[i]
        ndx = neigh[i]
        ssum = np.sum(E[sdx,:],axis=0)
        ssum[ssum > 0] = 1
        syn_count += ssum
        nsum = np.sum(E[ndx,:],axis=0)
        nsum[nsum > 0] = 1
        neigh_count += nsum
        diff = ssum - nsum
        tmp = E[sdx,:] - nsum 
        #print('gene diff',np.where(diff > 0))
        #print(sdx,ndx,np.where(diff > 1))

    diff = syn_count - neigh_count
    #print('Sum gene diff',np.where(diff > 0))
    #print(k)
    #reldff = 0.5*(syn_count - neigh_count) / (syn_count + neigh_count)
    #print('Rel diff', np.where(reldff))
    return np.where(diff > 0)[0].tolist()

def gene_profile(E,syn,neigh):
    """
    Returns dictionary of cam profiles for postsynaptic partners
    Dict format = {cell:cam_feature_vector}
    
    Parameters:
    -----------
    E : numpy array
     Expression matrix
    syn : list
      list of synaptic partners at each synapse
    neigh : list
      list of (nonsynaptic) neighbors at each synapse

    Note: syn[i] and neigh[i] correpspond to the ith synapse

    """

    (n,m) = E.shape
    k = len(syn)
    profile = defaultdict(lambda:np.zeros(m))
    syn_count = defaultdict(int)
    for i in range(k):
        ssum = np.sum(E[syn[i],:],axis=0)
        ssum[ssum > 0] = 1
        nsum = np.sum(E[neigh[i],:],axis=0)
        nsum[nsum > 0] = 1
        diff = ssum - nsum
        diff[diff < 1] = 0
        for j in syn[i]: 
            profile[j] += diff
            syn_count[j] += 1
    
    for j in profile: profile[j] /= syn_count[j]
    
    return profile



def get_synapse_data(S,e,cpartners=set([]),screen=None,remove_partners = False): 
    """
    Formats the synapse data. Returns list synapse (syn) and neighbors (neigh)
    where cell names have been converted to cell indicies

    Parameters:
    -----------
    S : dictionary
     Dictionary synapse data
    e : expression matrix object
    cpartners : set (optional)
      List of synaptic partners to remove from all neighbors
    screen : string (optional)
      Screen synapse based on data set. Will screen based on image name. Suggest using
      'N2U' for adult synapses and 'JSH' for L4 synapses.
    """
    
    syn,neigh,cneigh = [],[],[]
    for cont in S:
        if screen and screen not in S[cont]['sections'][0]: continue
        partners = set(S[cont]['partners'])
        neighbors = set(S[cont]['neighbors'])
        nonsyn = neighbors - partners
        if remove_partners:
            nonsyn = nonsyn - cpartners
        _cneigh = neighbors & cpartners
        syn.append([e.cells[n] for n in partners if n in e.cells])
        neigh.append([e.cells[n] for n in nonsyn if n in e.cells])
        cneigh.append([e.cells[n] for n in _cneigh if n in e.cells]) 
    return syn,neigh,cneigh


def score_overlap(sig,test):
    """
    Scores the overlap between the gene signature and the test signature

    Parameters:
    -----------
    sig: set
     Gene signature
    test: set
     Test signature
    """

    num = len(sig & test)
    den = float(len(sig))
    return num / den

def get_overlap(sig,E,syn,neigh):
    """
    Returns the overlap between the computed the gene signature and the
    gene expression of synaptic and (nonsynaptic) partners.

    Paramters:
    ----------
    sig : set
     Gene signature
    E : numpy array
     Expression matrix
    syn : list
     List of synaptic partners at each synapse
    neigh : list
     List of neighbors at each synapse

    Return:
    -------
    ssig : list
     List of overlap scores for each synaptic partner at each synapse
    neigh : list
     List of overlap score for each neighbor at each synapse
    idsyn : float
     Fraction of synapses where the highest overlap score is a synaptic partner

    """
    
    k = len(syn)
    den = float(len(sig))
    sig = set(sig)
    ssig,nsig = [],[]
    idsyn = 0
    for i in range(k):
        synscore = [0]
        neighscore = [0]
        for j in syn[i]:
            _ssig = set(np.where(E[j,:] > 0)[0].tolist())
            score = score_overlap(sig,_ssig)
            synscore.append(score)
            ssig.append(score)

        for j in neigh[i]:
            _nsig = set(np.where(E[j,:]> 0)[0].tolist())
            score = score_overlap(sig,_nsig)
            neighscore.append(score)
            nsig.append(score)
        if max(synscore) > max(neighscore): idsyn += 1
    
     
    return ssig,nsig,idsyn/float(k)

def get_overlap_spatial_loc(sig,E,syn,neigh,cneigh):
    """
    Returns the overlap between the computed the gene signature and the
    gene expression of synaptic and (nonsynaptic) partners.

    Paramters:
    ----------
    sig : set
     Gene signature
    E : numpy array
     Expression matrix
    syn : list
     List of synaptic partners at each synapse
    neigh : list
     List of neighbors at each synapse
    cneigh : lsit
     List of neighbors that are synaptic partners elsewhere

    Return:
    -------
    ssig : list
     List of overlap scores for each synaptic partner at each synapse
    neigh : list
     List of overlap score for each neighbor at each synapse
    idsyn : float
     Fraction of synapses where the highest overlap score is a synaptic partner

    """
    
    k = len(syn)
    den = float(len(sig))
    sig = set(sig)
    ssig,nsig = [],[]
    idsyn = 0
    for i in range(k):
        synscore = [0]
        neighscore = [0]
        for j in syn[i]:
            _ssig = set(np.where(E[j,:] > 0)[0].tolist())
            score = score_overlap(sig,_ssig)
            synscore.append(score)
            ssig.append(score)

        for j in neigh[i]:
            _nsig = set(np.where(E[j,:]> 0)[0].tolist())
            if j in cneigh[i]:
                _nsig = _nsig - sig
            score = score_overlap(sig,_nsig)
            neighscore.append(score)
            nsig.append(score)
        if max(synscore) > max(neighscore): idsyn += 1
    
     
    return ssig,nsig,idsyn/float(k)


def get_random_overlap(syn,neigh):
    """
    Returns the probability of randomly selecting the correct synaptic partners  

    Paramters:
    ----------
    syn : list
     List of synaptic partners at each synapse
    neigh : list
     List of neighbors at each synapse

    Return:
    -------
    idsyn : float
     Fraction of synapses where the highest overlap score is a synaptic partner

    """

    k = len(syn)
    idsyn = 0
    for i in range(k):
        num = float(len(syn[i]))
        den = num + len(neigh[i]) 
        r = random()
        thresh = num /den
        #print(r,thresh)
        if r <= thresh: 
            #print(r,thresh,syn[i],neigh[i])
            idsyn += 1
    #print(idsyn,k) 
    return idsyn/float(k)



