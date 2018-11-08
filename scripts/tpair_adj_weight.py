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
import LUT
import aux

lr_pairs = './mat/lr_dict.txt'

def run(fout=None):
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
    lntd,rntd = get_corresponding_edge_attr(n2u.A,rn2u)
    lntd,rntd = filter_corresponding_tds(lntd,rntd,[ntdbound,ntdbound])

    jsh = from_db(JSH,adjacency=True,remove=_remove)
    rjsh = jsh.A.map_vertex_names(lrd)
    jtdbound = get_td_bounds(jsh.A.es['weight'])
    ljtd,rjtd = get_corresponding_edge_attr(jsh.A,rjsh)
    ljtd,rjtd = filter_corresponding_tds(ljtd,rjtd,[jtdbound,jtdbound])

    bntd,bjtd = get_corresponding_edge_attr(n2u.A,jsh.A)
    bntd,bjtd = filter_corresponding_tds(bntd,bjtd,[ntdbound,jtdbound])

    ndelta = np.array(lntd) - np.array(rntd)
    jdelta = np.array(ljtd) - np.array(rjtd)
    bdelta = np.array(bntd) - np.array(bjtd)

    data = [ndelta,jdelta,bdelta]

    print(ttest_rel(lntd,rntd))
    print(ttest_rel(ljtd,rjtd))
    print(ttest_rel(bntd,bjtd))

    tval1,pval1 = ttest_ind(ndelta,jdelta)
    tval2,pval2 = ttest_ind(jdelta,bdelta)
    tval3,pval3 = ttest_ind(ndelta,bdelta)

    pval = [(0,1,pval1),(1,2,pval2)]

    fig,ax = plt.subplots(1,1,figsize=(12,10))
    tpair_adj_weight(ax,data,pval,fout=fout)
    plt.show()

if __name__ == '__main__':
    run()

