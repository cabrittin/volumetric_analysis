"""
tpair_confrac.py

Plot distribution of the differences between the arcsine of 
connectivity fractions of homologous neurons.

created: Christopher Brittin
date: 01 November 2018

"""
import sys
sys.path.append(r'./volumetric_analysis')
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import ttest_rel
from scipy.stats import ttest_ind

#Brittin modules
from connectome.load import from_db
from networks.stats import *
from figures.stats import *
import aux

mpl.rcParams['xtick.labelsize'] = 32
mpl.rcParams['ytick.labelsize'] = 32

lr_pairs = './mat/lr_neurons.txt'

def run(fout=None):
    N2U = 'N2U'
    JSH = 'JSH'
    lr = aux.read.into_list2(lr_pairs)
    lr = list(zip(*lr))
    _remove = ['VC01','VD01','VB01','VB02']

    n2u = from_db(N2U,adjacency=True,chemical=True,
                  electrical=True,remove=_remove)
    lncf = get_cf(n2u,vertices = lr[0],_arcsine=True)
    rncf = get_cf(n2u,vertices = lr[1],_arcsine=True)
    
    jsh = from_db(JSH,adjacency=True,chemical=True,
                  electrical=True,remove=_remove)
    ljcf = get_cf(jsh,vertices = lr[0],_arcsine=True)
    rjcf = get_cf(jsh,vertices = lr[1],_arcsine=True)

    cells = sorted((set(n2u.neurons)&set(jsh.neurons))-set(_remove))
    bncf = get_cf(n2u,vertices = cells,_arcsine=True)
    bjcf = get_cf(jsh,vertices = cells,_arcsine=True)
    
    data = [np.array(lncf['gap']) - np.array(rncf['gap']),
            np.array(ljcf['gap']) - np.array(rjcf['gap']),
            np.array(bncf['gap']) - np.array(bjcf['gap']),
            np.array(lncf['pre']) - np.array(rncf['pre']),
            np.array(ljcf['pre']) - np.array(rjcf['pre']),
            np.array(bncf['pre']) - np.array(bjcf['pre']),
            np.array(lncf['post']) - np.array(rncf['post']),
            np.array(ljcf['post']) - np.array(rjcf['post']),
            np.array(bncf['post']) - np.array(bjcf['post'])]

    print(ttest_rel(lncf['pre'],rncf['pre']))
    print(ttest_rel(ljcf['pre'],rjcf['pre']))
    print(ttest_rel(bncf['pre'],bjcf['pre']))

    print(ttest_rel(lncf['post'],rncf['post']))
    print(ttest_rel(ljcf['post'],rjcf['post']))
    print(ttest_rel(bncf['post'],bjcf['post']))

    tval0,pval0 = ttest_ind(data[0],data[2])
    tval1,pval1 = ttest_ind(data[0],data[1])
    tval1,pval2 = ttest_ind(data[1],data[2])
    tval1,pval3 = ttest_ind(data[3],data[4])
    tval1,pval4 = ttest_ind(data[4],data[5])
    tval1,pval5 = ttest_ind(data[6],data[7])
    tval1,pval6 = ttest_ind(data[7],data[8])

    pval = [
        (0,1,pval1),(1,2,pval2),
        (3,4,pval3),(4,5,pval4),
        (6,7,pval5),(7,8,pval6)]

    fig,ax = plt.subplots(1,1,figsize=(12,10))
    #tpair_confrac(ax,data,pval,fout=fout)
    tpair_syn(ax,data,pval,
              fout=fout,
              ylabel='arcsine(CF) difference',
              title = 'Homologous connectivity fractions (CF)',
              xticklabels = ['$C^{\mathrm{gap}}$',
                             '$C^{\mathrm{pre}}$',
                             '$C^{\mathrm{post}}$'],
              ylim=[-1,1])
    plt.tight_layout()
    plt.show()

if __name__=="__main__":
    run()
