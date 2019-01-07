"""
tpair_adj_deg.py

Plots distribution of the differences in homologous adjacency weight

crated: Christopher Brittin
data: 01 November 2018

"""
import sys
sys.path.append(r'./volumetric_analysis')
import matplotlib.pyplot as plt
from scipy.stats import ttest_rel
from scipy.stats import ttest_ind

#Brittin modules
from connectome.load import from_db
from networks.stats import *
from figures.stats import *
import aux

lr_pairs = './mat/lr_dict.txt'

def write_out(edges,score1,score2,_fout):
    _out = []
    for i in range(len(edges)):
        e = edges[i]
        _out.append([e,score1[i],score2[i],abs(score1[i]-score2[i])])
    aux.write.from_list(_fout,_out)

def run(fout=None,source_data=None):
    N2U = 'N2U'
    JSH = 'JSH'
    _lrd = aux.read.into_dict(lr_pairs)
    lrd = {}
    for key,val in _lrd.items():
        lrd[key] = val
        lrd[val] = key
        
    _remove = ['VC01','VD01','VB01','VB02','HSNL','HSNR','PVNL','PVNR']
    n2u = from_db(N2U,adjacency=True,remove=_remove)
    rn2u = n2u.A.map_vertex_names(lrd)
    ntdbound = get_td_bounds(n2u.A.es['weight'])
    _edges,lntd,rntd = get_corresponding_edge_attr(n2u.A,rn2u)
    _edges,lntd,rntd = filter_corresponding_tds(_edges,lntd,rntd,[ntdbound,ntdbound])
    if source_data:
        fsplit = source_data.split('.')
        dout = fsplit[0] + '_adult_contralateral.' + fsplit[1]
        write_out(_edges,lntd,rntd,dout)

    jsh = from_db(JSH,adjacency=True,remove=_remove)
    rjsh = jsh.A.map_vertex_names(lrd)
    jtdbound = get_td_bounds(jsh.A.es['weight'])
    _edges,ljtd,rjtd = get_corresponding_edge_attr(jsh.A,rjsh)
    _edges,ljtd,rjtd = filter_corresponding_tds(_edges,ljtd,rjtd,[jtdbound,jtdbound])
    if source_data:
        fsplit = source_data.split('.')
        dout = fsplit[0] + '_l4_contralateral.' + fsplit[1]
        write_out(_edges,ljtd,rjtd,dout)

    
    _edges,bntd,bjtd = get_corresponding_edge_attr(n2u.A,jsh.A)
    _edges,bntd,bjtd = filter_corresponding_tds(_edges,bntd,bjtd,[ntdbound,jtdbound])
    if source_data:
        fsplit = source_data.split('.')
        dout = fsplit[0] + '_adult_l4_homologous.' + fsplit[1]
        write_out(_edges,bntd,bjtd,dout)     
    
    ndelta = np.array(lntd) - np.array(rntd)
    jdelta = np.array(ljtd) - np.array(rjtd)
    bdelta = np.array(bntd) - np.array(bjtd)

    data = [ndelta,jdelta,bdelta]

    print('Stats:')
    print_wilcoxon(ndelta,'Adult L/R')
    print_wilcoxon(jdelta,'L4 L/R')
    print_wilcoxon(bdelta,'Adult/L4')
    
    tval1,pval1 = ttest_ind(ndelta,jdelta)
    tval2,pval2 = ttest_ind(jdelta,bdelta)
    tval3,pval3 = ttest_ind(ndelta,bdelta)

    pval = [(0,1,pval1),(1,2,pval2)]

    fig,ax = plt.subplots(1,1,figsize=(12,10))
    tpair_adj_weight(ax,data,pval,fout=fout)
    plt.show()

if __name__ == '__main__':
    run()

