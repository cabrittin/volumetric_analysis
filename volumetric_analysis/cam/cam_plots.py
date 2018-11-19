
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib as mpl
import numpy as np

mpl.rcParams['xtick.labelsize'] = 24
mpl.rcParams['ytick.labelsize'] = 24

PRESYN_COL = '#FA5252'
POSTSYN_COL = '#F9E530'
GAPJUNC_COL = '#47A705'

def plot_expression_patterns(ax,Exp,fout=None):
    cmap = ListedColormap(['w', 'k'])
    ax.matshow(Exp.E,cmap=cmap)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xlabel('CAM gene',fontsize=32)
    ax.set_ylabel('Neuron',fontsize=32)
    ax.set_title('Expression patterns',fontsize=32)
    plt.tight_layout()
    if fout: plt.savefig(fout)
    
def plot_genes_in_neurons(ax,Exp,fout=None):
    gpn = Exp.gene_per_neuron_count()
    plot_bar(ax,gpn,'# CAM genes expressed in given neuron','# neurons')
    plt.tight_layout()
    if fout: plt.savefig(fout)

def plot_isoforms_per_gene(ax,Exp,fout=None):
    ipg = Exp.isoform_per_gene_count()
    plot_bar(ax,ipg,'# isoforms', '# genes')
    ax.set_ylim([0,25])
    ax.set_xlim([0.5,10])    
    plt.tight_layout()
    if fout: plt.savefig(fout)
    
def plot_neurons_expressing_genes(ax,Exp,fout=None):
    npg = Exp.neuron_per_gene_count()
    plot_hist(ax,npg,'# neurons expressing a given CAM gene', 'Distribution')
    plt.tight_layout()
    if fout: plt.savefig(fout)

def plot_alt_spliced_genes(ax,Exp,nclass,fout=None):
    data = Exp.alt_splice_bilateral_dist(nclass)
    for [n,genes] in data:
        if genes:
            genes = ','.join(genes)
            
    count = [(g,c)
             for (g,c) in sorted(Exp.alt_gene_count.items())]
    [x,y] = zip(*count)
    
    plot_bar_chart(ax,x,y,
                   ylabel = '# bilateral neurons expressing gene',
                   title = ('Alt. spliced genes expressed in '
                            'bilaterally symmetric neurons'),
                   fs=21)
    for label in ax.get_xticklabels():
        label.set_fontstyle('italic')
    plt.tight_layout()
    if fout: plt.savefig(fout)

def plot_bar_chart(ax,x,y,width=0.45,ylabel=None,title=None,
                   color='k',facecolor='#CAC9C9',linewidth=3,
                   fs=32):
    N = len(x)
    ind = np.arange(N)
    rects = ax.bar(ind,y,width,color=color,facecolor=facecolor,
                   linewidth=linewidth)
    ax.set_xticks(ind + width/2)
    ax.set_xticklabels(x)
    labels = ax.get_xticklabels()
    plt.setp(labels, rotation=30)
    ax.set_xlim([-1,N])
    if ylabel: ax.set_ylabel(ylabel,fontsize=fs)
    if title: ax.set_title(title,fontsize=fs)
                   
    
def plot_bar(ax,data,xlabel,ylabel):
    data_x = np.arange(len(data))
    ax.bar(data_x,data,align='center',color='black')
    ax.set_xlabel(xlabel,fontsize=24)
    ax.set_ylabel(ylabel,fontsize=24)
    ax.set_xticks(data_x)
    ax.set_xlim([-1,len(data)])

def plot_hist(ax,data,xlabel,ylabel):
    nbins = np.max(data)
    ax.hist(data,bins=nbins,range=(0,nbins),
            histtype='step',linewidth=4,color='k',
            cumulative=1,
            normed = 1)
    ax.set_xlabel(xlabel,fontsize=24)
    ax.set_ylabel(ylabel,fontsize=24)
    ax.set_ylim([0,1])
    ax.set_xlim([0,nbins])


def plot_lr_discrepancies(ax,n2u,jsh,fout=None):
    N = 2
    ind = np.arange(N)
    width = 0.25

    rects1 = ax.bar(ind,n2u,width,color='#5A5A5A')
    rects2 = ax.bar(ind+width,jsh,width,color='#A3A3A3')

    ax.set_ylabel('Fraction of network',fontsize=28)
    ax.set_title('Inconsistent left/right connectivity',fontsize=28)
    #ax.set_xticks(ind + width / 2)
    ax.set_xticks(ind + width)
    ax.set_xticklabels(('neurons','synaptic partners'))
    ax.set_ylim([0,0.5])
    ax.set_xlim([-0.5,2])
    ax.legend((rects1[0],rects2[0]),('Adult','L4'),fontsize=24)

    autolabel(ax,rects1)
    autolabel(ax,rects2)
    plt.tight_layout()
    if fout: plt.savefig(fout)

def plot_cam_lus(ax,models,yerr,fout=None):
    error_kw = {'capsize': 10, 'capthick': 2,
                'ecolor': 'black','elinewidth':2}
    N = len(models[0])
    
    ind = np.arange(N)
    width = 0.25

    rects1 = ax.bar(ind,models[0],width,color=PRESYN_COL,
                   yerr=yerr[0],error_kw=error_kw)
    rects2 = ax.bar(ind+width,models[1],width,color=GAPJUNC_COL,
                    yerr=yerr[1],error_kw=error_kw)
    ax.set_ylabel('LUS',fontsize=28)
    ax.set_xticks(ind + width*0.5)
    ax.set_xticklabels(('WBE','SBE','IE'))
    ax.set_ylim([0,1.])
    ax.set_xlim([-0.5,N])
    ax.legend((rects1[0], rects2[0]), ('Presyn.', 'Gap J.'),fontsize=24)
    #autolabel(ax,rects1)
    #autolabel(ax,rects2)
    ax.set_title('Combinatorial CAM models',fontsize=28)
    plt.tight_layout()
    if fout: plt.savefig(fout)
    
def autolabel(ax,rects):
    """
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                '%1.2f' % float(height),
                ha='center', va='bottom',fontsize=18)


