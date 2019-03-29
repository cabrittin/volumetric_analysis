"""
compare_neigh_overlap.py

Plots distributions of Jaccard distances for overlapping ipsilateral 
neighborhoods (blue) and homologous contralateral neighborhoods (red) 
in the adult and L4.

crated: Christopher Brittin
data: 01 November 2018

"""
import sys
sys.path.append(r'./volumetric_analysis')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ttest_1samp
from scipy.stats import ttest_ind
from scipy.stats import mannwhitneyu

#Brittin modules
from connectome.load import from_db
from networks.stats import *
from figures.stats import *
import aux

lr_dict = './mat/lr_dict.txt'
left_nodes = './mat/left_nodes.txt'
right_nodes = './mat/right_nodes.txt'

def write_out(fsim,fover,G1,G2,nodes1,nodes2):
    sim = get_neighborhood_similarity(A,reflected,nodes1)
    overlap = get_neighborhood_overlap_similarity(A,left+right)

def get_data(A,reflected,left,right):
    sim = get_neighborhood_similarity(A,reflected,left) 
    sim = [s for s in sim if s > -1]
    overlap = get_neighborhood_overlap_similarity(A,left+right)
    overlap = [o for o in overlap if o > -1]
    return np.array(sim),np.array(overlap)

def arcsine(data):
    return 2*np.arcsin(np.sqrt(data))

def run(fout=None,source_data=None):
    N2U = 'N2U'
    JSH = 'JSH'
    left = aux.read.into_list(left_nodes)
    right = aux.read.into_list(right_nodes)
    _lrd = aux.read.into_dict(lr_dict)
    lrd = {}
    for key,val in _lrd.items():
        lrd[key] = val
        lrd[val] = key
    _remove = ['VC01','VD01','VB01','VB02','HSNL','HSNR','PVNL','PVNR']

    n2u = from_db(N2U,adjacency=True,remove=_remove)
    a = n2u.A.copy()
    reflected = n2u.A.map_vertex_names(lrd)
    vertices = left + right
    ns,no = get_data(n2u.A,reflected,left,right)
    if source_data:
        fsplit = source_data.split('.')
        fsim = fsplit[0] + '_adult_homologous.' + fsplit[1]
        fovr = fsplit[0] + '_adult_overlap.' + fsplit[1]
        sim = get_neighborhood_similarity(n2u.A,reflected,left)
        dsim = []
        for i in range(len(left)):
            dsim.append([left[i],right[i],sim[i]])
        both = left+right
        overlap = get_neighborhood_overlap_similarity(n2u.A,left+right)
        dovr = []
        for i in range(len(both)):
            dovr.append([both[i],overlap[i]])
        aux.write.from_list(fsim,dsim)
        aux.write.from_list(fovr,dovr)
            
    jsh = from_db(JSH,adjacency=True,remove=_remove)
    reflected = jsh.A.map_vertex_names(lrd)
    vertices = left + right
    js,jo = get_data(jsh.A,reflected,left,right)
    #js,jo = arcsine(js),arcsine(jo)
    if source_data:
        fsplit = source_data.split('.')
        fsim = fsplit[0] + '_l4_homologous.' + fsplit[1]
        fovr = fsplit[0] + '_l4_overlap.' + fsplit[1]
        sim = get_neighborhood_similarity(jsh.A,reflected,left)
        dsim = []
        for i in range(len(left)):
            dsim.append([left[i],right[i],sim[i]])
        both = left+right
        overlap = get_neighborhood_overlap_similarity(jsh.A,left+right)
        dovr = []
        for i in range(len(both)):
            dovr.append([both[i],overlap[i]])
        aux.write.from_list(fsim,dsim)
        aux.write.from_list(fovr,dovr)
        
    vertices = sorted((set(n2u.neurons)&set(jsh.neurons))-set(_remove))
    sim = get_neighborhood_similarity(n2u.A,jsh.A,vertices)
    noverlap = get_neighborhood_overlap_similarity(n2u.A,vertices)
    joverlap = get_neighborhood_overlap_similarity(jsh.A,vertices)
    bs = np.array([s for s in sim if s > -1])
    overlap = noverlap + joverlap
    bo = np.array([o for o in overlap if o > -1])
    if source_data:
        fsplit = source_data.split('.')
        fsim = fsplit[0] + '_adult_l4_homologous.' + fsplit[1]
        fovr = fsplit[0] + '_adult_l4_overlap.' + fsplit[1]        
        dsim = []
        for i in range(len(vertices)):
            dsim.append([vertices[i],vertices[i],sim[i]])
        dovr = []
        for i in range(len(vertices)):
            dovr.append(['adult-'+vertices[i],noverlap[i]])
            dovr.append(['l4-'+vertices[i],joverlap[i]])
        aux.write.from_list(fsim,dsim)
        aux.write.from_list(fovr,dovr)
        
    #tval1,upval1 = ttest_ind(ns,no)
    #tval2,upval2 = ttest_ind(js,jo)
    #tval3,upval3 = ttest_ind(bs,bo)
    tval1,upval1 = mannwhitneyu(ns,no)
    tval2,upval2 = mannwhitneyu(js,jo)
    tval3,upval3 = mannwhitneyu(bs,bo)

    print(upval1,upval2,upval3)
    
    data = [ns,no,js,jo,bs,bo]
    pval = [(0,1,upval1),(2,3,upval2),(4,5,upval3)]
    
    print('Stats:')
    print_wilcoxon(data[0],'Adult L/R contra.')
    print_wilcoxon(data[1],'Adult L/R overlap')
    print_wilcoxon(data[2],'L4 L/R contra.')
    print_wilcoxon(data[3],'L4 L/R overlap')
    print_wilcoxon(data[4],'Adult/L4 contra.')
    print_wilcoxon(data[5],'Adult/L4 overlap')

    fig,ax = plt.subplots(1,1,figsize=(12,10))
    plot_overlap_compare(ax,data,pval,fout=fout)
    plt.show()


if __name__ == '__main__':
    run()

