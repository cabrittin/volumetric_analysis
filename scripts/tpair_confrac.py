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
from scipy.stats import wilcoxon
from scipy.stats import mannwhitneyu

#Brittin modules
from connectome.load import from_db
from networks.stats import *
from figures.stats import *
import aux

mpl.rcParams['xtick.labelsize'] = 32
mpl.rcParams['ytick.labelsize'] = 32

lr_pairs = './mat/lr_neurons.txt'

ARCSINE = False

def run(fout=None):
    N2U = 'N2U'
    JSH = 'JSH'
    lr = aux.read.into_list2(lr_pairs)
    lr = list(zip(*lr))
    _remove = ['VC01','VD01','VB01','VB02']

    n2u = from_db(N2U,adjacency=True,chemical=True,
                  electrical=True,remove=_remove)
    lncf = get_cf(n2u,vertices = lr[0],_arcsine=ARCSINE)
    rncf = get_cf(n2u,vertices = lr[1],_arcsine=ARCSINE)
    ancf = get_cf(n2u,vertices = set(n2u.neurons) - set(_remove),_arcsine=ARCSINE)
    nstd = {'gap': np.std(ancf['gap']),
            'pre': np.std(ancf['pre']),
            'post': np.std(ancf['post'])}

    jsh = from_db(JSH,adjacency=True,chemical=True,
                  electrical=True,remove=_remove)
    ljcf = get_cf(jsh,vertices = lr[0],_arcsine=ARCSINE)
    rjcf = get_cf(jsh,vertices = lr[1],_arcsine=ARCSINE)
    ajcf = get_cf(jsh,vertices = set(jsh.neurons) - set(_remove),_arcsine=ARCSINE)
    jstd = {'gap' : np.std(ajcf['gap']),
            'pre' : np.std(ajcf['pre']),
            'post': np.std(ajcf['post'])}

    cells = sorted((set(n2u.neurons)&set(jsh.neurons))-set(_remove))
    bncf = get_cf(n2u,vertices = cells,_arcsine=ARCSINE) 
    bjcf = get_cf(jsh,vertices = cells,_arcsine=ARCSINE)
    bstd = {'gap' : np.std(np.concatenate((ancf['gap'],ajcf['gap']))),
            'pre' : np.std(np.concatenate((ancf['pre'],ajcf['pre']))),
            'post': np.std(np.concatenate((ancf['post'],ajcf['post'])))}



    data = [(lncf['gap'] - rncf['gap'])/nstd['gap'],
            (ljcf['gap'] - rjcf['gap'])/jstd['gap'],
            (bncf['gap'] - bjcf['gap'])/bstd['gap'],
            (lncf['pre'] - rncf['pre'])/nstd['pre'],
            (ljcf['pre'] - rjcf['pre'])/jstd['pre'],
            (bncf['pre'] - bjcf['pre'])/bstd['pre'],
            (lncf['post'] - rncf['post'])/nstd['post'],
            (ljcf['post'] - rjcf['post'])/jstd['post'],
            (bncf['post'] - bjcf['post'])/bstd['post']]

    print('Stats:')
    print_wilcoxon(data[0],'Adult L/R gap')
    print_wilcoxon(data[1],'L4 L/R gap')
    print_wilcoxon(data[2],'Adult/L4 gap')
    print_wilcoxon(data[3],'Adult L/R pre')
    print_wilcoxon(data[4],'L4 L/R pre')
    print_wilcoxon(data[5],'Adult/L4 pre')
    print_wilcoxon(data[6],'Adult L/R post')
    print_wilcoxon(data[7],'L4 L/R post')
    print_wilcoxon(data[8],'Adult/L4 post')


    tval0,pval0 = mannwhitneyu(data[0],data[2])
    tval1,pval1 = mannwhitneyu(data[0],data[1])
    tval1,pval2 = mannwhitneyu(data[1],data[2])
    tval1,pval3 = mannwhitneyu(data[3],data[4])
    tval1,pval4 = mannwhitneyu(data[4],data[5])
    tval1,pval5 = mannwhitneyu(data[6],data[7])
    tval1,pval6 = mannwhitneyu(data[7],data[8])

    pval = [
        (0,1,pval1),(1,2,pval2),
        (3,4,pval3),(4,5,pval4),
        (6,7,pval5),(7,8,pval6)]

    fig,ax = plt.subplots(1,1,figsize=(12,10))
    #tpair_confrac(ax,data,pval,fout=fout)
    tpair_syn(ax,data,pval,
              fout=fout,
              ylabel='Normalized CF difference',
              title = 'Homologous connectivity fractions (CF)',
              xticklabels = ['$C^{\mathrm{gap}}$',
                             '$C^{\mathrm{pre}}$',
                             '$C^{\mathrm{post}}$'],
              ylim=[-3,3])

    plt.tight_layout()
    plt.show()

if __name__=="__main__":
    run()
