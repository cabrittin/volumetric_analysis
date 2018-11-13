import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.mlab as mlab
from matplotlib_venn import venn2
import numpy as np
from scipy.stats import skewnorm
from scipy.optimize import curve_fit
from scipy.stats import norm
from sklearn.neighbors import KernelDensity
from sklearn.mixture import GaussianMixture
import scipy.stats as scpstats
from matplotlib.ticker import MaxNLocator
from scipy.stats import pearsonr


mpl.rcParams['xtick.labelsize'] = 24
mpl.rcParams['ytick.labelsize'] = 24

PRESYN_COL = '#FA5252'
POSTSYN_COL = '#F9E530'
GAPJUNC_COL = '#47A705'
ALL_COL = 'k'

ADULT_COL = '#FFC300'
L4_COL = '#3380FF'
AL_COL = '#14FA29'

SEN_COL = '#910000'
INT_COL = '#1700FF'
MOT_COL = '#1D6500'

def plot_adj_degree(ax,deg,fout=None,nbins=100,hrange=(0,200),
                    xlim=[0,120],ylim=None,density=True,fit_mode=None,
                    xlabel=None,ylabel=None,xfs=32,yfs=32):
    """
    Plot adjacency degree.

    Args:
       ax: matplolib axis handle
       deg: list of degrees

    """
    y,bins,_ = ax.hist(deg,bins=nbins,range=hrange,
                       histtype='step',density=density,
                       cumulative=False,color='k',
                       linewidth=3,label='empirical data')

    y = np.array(deg).reshape(-1,1)
    x = np.linspace(0,120,121).reshape(-1,1)
    if fit_mode=='KDE':
        kde = KernelDensity(kernel='gaussian',bandwidth=5.).fit(y)
        log_dens = kde.score_samples(x)
        ax.plot(x,np.exp(log_dens),'r-',label='fit',linewidth=6)

    if fit_mode=='GMM':
        g = GaussianMixture(n_components=3,covariance_type='full')
        g.fit(y)
        log_dens = g.score_samples(x)
        ax.plot(x,np.exp(log_dens),'r-',label='GMM fit (n=3)',linewidth=6)
        
    ax.legend(loc='upper left',fontsize=min(24,yfs))
    if xlim: ax.set_xlim(xlim)
    if ylim: ax.set_ylim(ylim)
    if xlabel: ax.set_xlabel(xlabel,fontsize=xfs)
    if ylabel: ax.set_ylabel(ylabel,fontsize=yfs)
    plt.tight_layout()
    if fout: plt.savefig(fout)


def plot_syn_degree(ax,deg,fout=None,nbins=100,hrange=(0,200),
                    xlim=[0,40],ylim=None,density=True,fit_mode=None,
                    xlabel=None,ylabel=None,xfs=38,yfs=38):
    """
    Plot adjacency degree.

    Args:
       ax: matplolib axis handle
       deg: list of degrees

    """
    y,bins,_ = ax.hist(deg,bins=nbins,range=hrange,
                       histtype='step',density=density,
                       cumulative=False,color='k',
                       linewidth=3,label='empirical data')

    y = np.array(deg).reshape(-1,1)
    x = np.linspace(0,120,121).reshape(-1,1)
    if fit_mode=='KDE':
        kde = KernelDensity(kernel='gaussian',bandwidth=5.).fit(y)
        log_dens = kde.score_samples(x)
        ax.plot(x,np.exp(log_dens),'r-',label='KDE fit',linewidth=6)

    if fit_mode=='GMM':
        g = GaussianMixture(n_components=3,covariance_type='full')
        g.fit(y)
        log_dens = g.score_samples(x)
        ax.plot(x,np.exp(log_dens),'r-',label='GMM fit (n=3)',linewidth=6)
        
    ax.legend(loc='upper right',fontsize=min(24,yfs))
    if xlim: ax.set_xlim(xlim)
    if ylim: ax.set_ylim(ylim)
    if xlabel: ax.set_xlabel(xlabel,fontsize=xfs)
    if ylabel: ax.set_ylabel(ylabel,fontsize=yfs)
    plt.tight_layout()
    if fout: plt.savefig(fout)    

def plot_adj_weight(ax,deg,fout=None,nbins=1000,hrange=(0,200),
                    xlim=[0,1000],ylim=None,density=True,fit_mode=None):
    """
    Plot adjacency weight.

    Args:
       ax: matplolib axis handle
       deg: list of degrees

    """
    y,bins,_ = ax.hist(deg,bins=nbins,range=hrange,
                       histtype='step',density=density,
                       cumulative=False,color='k',
                       linewidth=3,label='empirical data')

    y = np.array(deg).reshape(-1,1)
    x = np.linspace(0,120,121).reshape(-1,1)
        
    #ax.legend(loc='upper left',fontsize=32)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel('Adjacency degree, $d$',fontsize=38)
    ax.set_ylabel('Probability',fontsize=38)
    plt.tight_layout()
    if fout: plt.savefig(fout)
    
def plot_dist(ax,deg,fout=None,nbins=100,hrange=(0,1),
              xlim=None,ylim=None,density=False,cumulative=False,
              xlabel=None,ylabel=None,fs=38):

    y,bins,_ = ax.hist(deg,bins=nbins,range=hrange,
                       histtype='step',density=density,
                       cumulative=cumulative,color='k',
                       linewidth=3,label='empirical data')

    if xlim: ax.set_xlim(xlim)
    if ylim: ax.set_ylim(ylim)
    if xlabel: ax.set_xlabel(xlabel,fontsize=fs)
    if ylabel: ax.set_ylabel(ylabel,fontsize=fs)
    if fout: plt.savefig(fout)
    
def plot_lognorm_probplot(ax,td,fout=None,fs=32,
                          ylabel=None):
    tdlog = np.log(td)
    a,b = scpstats.probplot(tdlog,dist='norm',plot=ax)
    ax.get_lines()[0].set_markerfacecolor('w')
    ax.get_lines()[0].set_markeredgecolor('k')
    ax.get_lines()[1].set_color('k')
    ax.get_lines()[1].set_linewidth(4)
    ax.axvspan(-2,2,facecolor='#9B9B9B',alpha=0.5)
    ax.set_title('Normal probability plot',fontsize=fs)
    if ylabel: ax.set_ylabel(ylabel,fontsize=fs)
    ax.set_xlabel('$\mathcal{N}(0,1)$ ordered statistic medians',
                  fontsize=fs)
    ax.text(0.30,0.15,'95% of physical contacts',
            verticalalignment='bottom',
            transform=ax.transAxes,fontsize=int(0.5*fs))
    if fout:plt.savefig(fout)


def plot_cf_dist(ax,cf,nbins=101,density=True,cumulative=False,
                 hrange=(0,1),linewidth=3,fs=38):
    
    
    ax.hist(cf[0],bins=nbins,range=hrange,histtype='step',
            density=density,cumulative=cumulative,color=PRESYN_COL,
            linewidth=3,label='$C^{\mathrm{pre}}$')
    
    ax.hist(cf[1],bins=nbins,range=hrange,histtype='step',
            density=density,cumulative=cumulative,color=POSTSYN_COL,
            linewidth=3,label='$C^{\mathrm{post}}$')

    ax.hist(cf[2],bins=nbins,range=hrange,histtype='step',
            density=density,cumulative=cumulative,color=GAPJUNC_COL,
            linewidth=3,label='$C^{\mathrm{post}}$')  
    
    ax.set_xlabel('Connectivity fraction',fontsize=fs)
    ax.set_ylabel('Fraction of neurons',fontsize=fs)

    ax.legend(fontsize=24)


def plot_cf_fit(ax,cf,nbins=100,linewidth=3,xlim=None,ylim=None,fs=38):
    x = np.linspace(0,1,nbins)
    muall = np.mean(cf['all'])
    mupre = np.mean(cf['pre'])
    mupost = np.mean(cf['post'])
    mugap = np.mean(cf['gap'])
    plot_skew_norm_fit(ax,cf['all'],x,color=ALL_COL,linewidth=linewidth,
                       label='$C^{\mathrm{all}}$, $\mu=$%1.2f'%muall)    
    plot_skew_norm_fit(ax,cf['gap'],x,color=GAPJUNC_COL,linewidth=linewidth,
                       label='$C^{\mathrm{gap}}$, $\mu=$%1.2f'%mugap)
    plot_skew_norm_fit(ax,cf['pre'],x,color=PRESYN_COL,linewidth=linewidth,
                       label='$C^{\mathrm{pre}}$, $\mu=$%1.2f'%mupre)
    plot_skew_norm_fit(ax,cf['post'],x,color=POSTSYN_COL,linewidth=linewidth,
                       label='$C^{\mathrm{post}}$, $\mu=$%1.2f'%mupost)
    if xlim: ax.set_xlim(xlim)
    if ylim: ax.set_ylim(ylim)
    ax.set_xlabel('Connectivity fraction',fontsize=fs)
    ax.set_ylabel('Probability density', fontsize=fs)
    ax.legend(fontsize=32)
    
def plot_skew_norm_fit(ax,data,bins,color='k',
                       linestyle='-',linewidth=3,label=None):
    s,loc,scale = skewnorm.fit(data)
    pdf = skewnorm.pdf(bins,s,scale =scale,loc=loc)
    ax.plot(bins,pdf,color=color,linestyle=linestyle,
            linewidth=linewidth,label=label) 

def plot_norm_fit(ax,data,bins,color='k',linestyle='--',linewidth=3,label=None):
    mu = np.mean(data)
    sig = np.std(data)
    pdf = mlab.normpdf(bins,mu,sig)
    ax.plot(bins, pdf,color=color,linestyle=linestyle,
            linewidth=linewidth,label=label)   

def plot_boxplots(ax,data,labels=None,pval=[],showfliers=False,
                  annotscale=1.25,positions=None,width=0.6,
                  xlim=None,ylim=None,xlabel=None,ylabel=None,
                  title=None,yfs=32,xfs=32,tfs=32,fout=None,
                  colors = None):
    flierprops = dict(markersize=4,marker='d',markerfacecolor='k')
    medianprops = dict(linestyle='-',linewidth=5,color='k')
    whiskerprops = dict(linestyle='-',linewidth=3,color='k')
    boxprops = dict(linewidth=2,color='k',facecolor='#ABABAB')
    capprops = dict(linewidth=3)
    bp = ax.boxplot(data,positions=positions,vert=True,
                    patch_artist=True,
                    labels=labels,
                    medianprops=medianprops,
                    whiskerprops=whiskerprops,
                    boxprops=boxprops,
                    capprops=capprops,
                    showfliers=showfliers,
                    flierprops=flierprops,
                    widths=width)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.tick_params(axis='x',direction='out')
    ax.tick_params(axis='y',length=0)
    ax.grid(axis='y',color='0.9',linestyle='-',linewidth=1)
    if colors:
        for patch,color in zip(bp['boxes'],colors):
            patch.set_facecolor(color)        
    
    for (i,j,p) in pval:
        print(stars(p))
        cap1 = bp['caps'][2*i+1].get_ydata()
        cap2 = bp['caps'][2*j+1].get_ydata()
        y_max = np.max(np.concatenate((cap1,cap2)))
        y_min = np.min(np.concatenate((cap1,cap2)))
        ax.annotate("", xy=(i+1, y_max), xycoords='data',
                    xytext=(j+1, y_max), textcoords='data',
                    arrowprops=dict(arrowstyle="-", ec='k',
                                    connectionstyle="bar,fraction=0.2"))
        ax.text(0.5*(i+j+2), y_max+annotscale, stars(p),
                fontsize=24,
                horizontalalignment='center',
                verticalalignment='center')  

    if ylim: ax.set_ylim(ylim)
    if xlim: ax.set_xlim(xlim)
    if ylabel: ax.set_ylabel(ylabel,fontsize=yfs)
    ax.tick_params(axis='x',labelsize=xfs)
    if title: ax.set_title(title,fontsize=tfs,y=1.04)
    if fout: plt.savefig(fout)
    return bp

def plot_overlap_compare(ax,data,pval,fout=None):
    labels=None
    pos = [1,2,3,4,5,6]
    colors = [SEN_COL,INT_COL,
              SEN_COL,INT_COL,
              SEN_COL,INT_COL]
    
    bp = plot_boxplots(ax,data,labels=labels,positions=pos,
                       pval=pval,annotscale=0.08,
                       ylabel='Jaccard distance',
                       title = 'Overlapping vs. homologous neighborhoods',
                       showfliers=True,width=0.2,colors=colors)

    ax.set_xticklabels(['Adult L\R',
                        'L4 L\R',
                        'Adult/L4'])
    ax.set_xticks([1.5,3.5,5.5])
    ax.axvspan(0,2.5,facecolor='#C3C3C3')
    ax.axvspan(2.5,4.5,facecolor='#D8D7D7')
    ax.axvspan(4.5,7,facecolor='#C3C3C3')
    ax.set_ylim([0,1])        
    _A, = ax.plot([1,1],SEN_COL)
    _L, = ax.plot([1,1],INT_COL)
    leg =ax.legend((_A, _L),('Homologous neighborhoods',
                             'Overlapping neighborhoods'),fontsize=18)
    for legobj in leg.legendHandles:
        legobj.set_linewidth(4.0)
    _A.set_visible(False)
    _L.set_visible(False)
    plt.tight_layout() 
    if fout: plt.savefig(fout)


def dist_adj_subgroups(ax,data,fout=None):
    labels=None
    pos = [1.5,2,2.5,3.5,4,4.5]
    colors = [SEN_COL,INT_COL,MOT_COL,
              SEN_COL,INT_COL,MOT_COL]
    
    bp = plot_boxplots(ax,data,labels=labels,positions=pos,
                       ylabel='Adjacency degree',
                       title='Adjacency degree by functional class',
                       showfliers=True,width=0.2,colors=colors)

    ax.set_xticklabels(['Adult','L4'])
    ax.set_xticks([2, 4])
    ax.axvspan(0,3,facecolor='#C3C3C3')
    ax.axvspan(3,5,facecolor='#D8D7D7')
    ax.set_ylim([0,120])        
    _A, = ax.plot([1,1],SEN_COL)
    _L, = ax.plot([1,1],INT_COL)
    _AL, = ax.plot([1,1],MOT_COL)
    leg =ax.legend((_A, _L,_AL),('Sensory', 'Interneuron','Motor'),
                   loc='upper center',fontsize=24)
    for legobj in leg.legendHandles:
        legobj.set_linewidth(4.0)
    _A.set_visible(False)
    _L.set_visible(False)
    _AL.set_visible(False)
    plt.tight_layout() 
    if fout: plt.savefig(fout)


def dist_adj_subgroups2(ax,data,fout=None):
    labels=None
    pos =[]
    colors = []
    N = int(len(data) / 2)
    dx = 0.15
    for i in range(1,N+1):
        pos.append(i-dx)
        pos.append(i+dx)
        colors.extend([ADULT_COL,L4_COL])

    
    bp = plot_boxplots(ax,data,labels=labels,positions=pos,
                       ylabel='Adjacency degree',
                       title='Adjacency degree by group',
                       showfliers=True,width=0.2,colors=colors)

    
    ax.set_xticklabels(['Sp','Sa','I1','I2','SMN','HMNp','HMNa'])
    ax.set_xticks([1,2,3,4,5,6,7])

    
    ax.axvspan(0,1.5,facecolor='#C3C3C3')
    ax.axvspan(1.5,2.5,facecolor='#D8D7D7')
    ax.axvspan(2.5,3.5,facecolor='#C3C3C3')
    ax.axvspan(3.5,4.5,facecolor='#D8D7D7')
    ax.axvspan(4.5,5.5,facecolor='#C3C3C3')
    ax.axvspan(5.5,6.5,facecolor='#D8D7D7')
    ax.axvspan(6.5,8,facecolor='#C3C3C3')
    
    ax.set_ylim([0,120])


    _A, = ax.plot([1,1],ADULT_COL)
    _L, = ax.plot([1,1],L4_COL)
    leg =ax.legend((_A, _L),('Adult', 'L4'),
                   loc='upper right',fontsize=24)
    for legobj in leg.legendHandles:
        legobj.set_linewidth(4.0)
    _A.set_visible(False)
    _L.set_visible(False)


    plt.tight_layout() 
    if fout: plt.savefig(fout)
    
def dist_syn_subgroups(ax,data,fout=None):
    labels=None

    pos = [1.5,2,2.5,3.5,4,4.5,5.5,6,6.5]
    colors = [SEN_COL,INT_COL,MOT_COL,
              SEN_COL,INT_COL,MOT_COL,
              SEN_COL,INT_COL,MOT_COL]
    
    bp = plot_boxplots(ax,data,labels=labels,positions=pos,
                       ylabel='Connectivity fraction',
                       title='Connectivity fraction by modality',
                       showfliers=True,width=0.2,colors=colors)

    ax.set_xticklabels(['gap j.',
                        'presyn.',
                        'postsyn.'])
    ax.set_xticks([2, 4, 6])
    ax.axvspan(0,3,facecolor='#C3C3C3')
    ax.axvspan(3,5,facecolor='#D8D7D7')
    ax.axvspan(5,8,facecolor='#C3C3C3')
    #ax.axhline(0,color='r',linewidth=3,linestyle='--')
    ax.set_ylim([0,40])        
    _A, = ax.plot([1,1],ADULT_COL)
    _L, = ax.plot([1,1],L4_COL)
    leg =ax.legend((_A, _L),('Adult', 'L4'),
                   loc='upper left',fontsize=24)
    for legobj in leg.legendHandles:
        legobj.set_linewidth(4.0)
    _A.set_visible(False)
    _L.set_visible(False)
    #ax.axhline(0,color='r',linewidth=3,linestyle='--')
    plt.tight_layout() 
    if fout: plt.savefig(fout)    
    
def plot_confrac_subgroups(ax,data,fout=None):
    labels=None
    pos = [1.5,2,2.5,3.5,4,4.5,5.5,6,6.5]
    colors = [SEN_COL,INT_COL,MOT_COL,
              SEN_COL,INT_COL,MOT_COL,
              SEN_COL,INT_COL,MOT_COL]
    
    bp = plot_boxplots(ax,data,labels=labels,positions=pos,
                       ylabel='Connectivity fraction',
                       title='Connectivity fraction by functional class',
                       showfliers=True,width=0.2,colors=colors)

    ax.set_xticklabels(['$C^{\mathrm{gap}}$',
                        '$C^{\mathrm{pre}}$',
                        '$C^{\mathrm{post}}$'])
    ax.set_xticks([2, 4, 6])
    ax.axvspan(0,3,facecolor='#C3C3C3')
    ax.axvspan(3,5,facecolor='#D8D7D7')
    ax.axvspan(5,8,facecolor='#C3C3C3')
    #ax.axhline(0,color='r',linewidth=3,linestyle='--')
    ax.set_ylim([0,0.6])        
    _A, = ax.plot([1,1],SEN_COL)
    _L, = ax.plot([1,1],INT_COL)
    _AL, = ax.plot([1,1],MOT_COL)
    leg =ax.legend((_A, _L,_AL),('Sensory', 'Interneuron','Motor'),
                   loc='upper left',fontsize=24)
    for legobj in leg.legendHandles:
        legobj.set_linewidth(4.0)
    _A.set_visible(False)
    _L.set_visible(False)
    _AL.set_visible(False)
    #ax.axhline(0,color='r',linewidth=3,linestyle='--')
    plt.tight_layout() 
    if fout: plt.savefig(fout)



def tpair_adj_deg(ax,data,pval,fout=None):
    labels = ['Adult L/R','L4 L/R','Adult/L4']
    pos = [1,2,3]
    colors = [ADULT_COL,L4_COL,AL_COL]
    
    bp = plot_boxplots(ax,data,labels=labels,positions=pos,
                       pval=pval,annotscale=10,
                       ylabel='Degree difference',
                       title = 'Homologous adjacency degree',
                       showfliers=True,width=0.2,colors=colors)
    ax.axhline(0,color='r',linewidth=3,linestyle='--')
    ax.set_ylim([-40,40])
    if fout: plt.savefig(fout)

def tpair_adj_weight(ax,data,pval,fout=None):
    labels = ['Adult L/R','L4 L/R','Adult/L4']
    pos = [1,2,3]
    colors = [ADULT_COL,L4_COL,AL_COL]
    
    bp = plot_boxplots(ax,data,labels=labels,positions=pos,
                       pval=pval,annotscale=1.0,
                       ylabel='log(surface area) difference',
                       title = 'Homologous surface area contacts',
                       showfliers=True,width=0.2,colors=colors)
    ax.axhline(0,color='r',linewidth=3,linestyle='--')
    ax.set_ylim([-4,4])
    if fout: plt.savefig(fout)

              
def tpair_syn(ax,data,pval,fout=None,ylabel=None,
              title=None,xticklabels=None,ylim=[-1,1]):
    labels=None
    pos = [1,2,3,4,5,6,7,8,9]
    colors = [ADULT_COL,L4_COL,AL_COL,
              ADULT_COL,L4_COL,AL_COL,
              ADULT_COL,L4_COL,AL_COL]
    
    bp = plot_boxplots(ax,data,labels=labels,positions=pos,
                       pval=pval,annotscale=0.1,
                       ylabel=ylabel,
                       title=title,
                       showfliers=True,width=0.2,colors=colors)

    if xticklabels: ax.set_xticklabels(xticklabels)
    ax.set_xticks([2, 5, 8])
    ax.axvspan(0,3.5,facecolor='#C3C3C3')
    ax.axvspan(3.5,6.5,facecolor='#D8D7D7')
    ax.axvspan(6.5,10,facecolor='#C3C3C3')
    ax.axhline(0,color='r',linewidth=3,linestyle='--')
    ax.set_ylim(ylim)        
    _A, = ax.plot([1,1],ADULT_COL)
    _L, = ax.plot([1,1],L4_COL)
    _AL, = ax.plot([1,1],AL_COL)
    leg =ax.legend((_A, _L,_AL),('Adult L/R', 'L4 L/R','Adult/L4'),fontsize=18)
    for legobj in leg.legendHandles:
        legobj.set_linewidth(4.0)
    _A.set_visible(False)
    _L.set_visible(False)
    _AL.set_visible(False)
    ax.axhline(0,color='r',linewidth=3,linestyle='--')
    plt.tight_layout() 
    if fout: plt.savefig(fout)
    

def tpair_syn_adj_ratio(ax,data,fout=None):
    labels=None
    pos = [1,2,3,4,5,6]
    colors = [ADULT_COL,L4_COL,
              ADULT_COL,L4_COL,
              ADULT_COL,L4_COL,]
    
    bp = plot_boxplots(ax,data,labels=labels,positions=pos,
                       ylabel='synapse-to-adjacency ratio',
                       title = 'Adjacency contact occupied by synapse',
                       showfliers=True,width=0.2,colors=colors)

    ax.set_xticklabels(['gap j.',
                        'presyn.',
                        'postsyn'])
    ax.set_xticks([1.5, 3.5, 5.5])
    ax.axvspan(0,2.5,facecolor='#C3C3C3')
    ax.axvspan(2.5,4.5,facecolor='#D8D7D7')
    ax.axvspan(4.5,7,facecolor='#C3C3C3')
    ax.set_ylim([0,1])        
    _A, = ax.plot([1,1],ADULT_COL)
    _L, = ax.plot([1,1],L4_COL)
    leg =ax.legend((_A, _L),('Adult', 'L4'),fontsize=18,loc='upper left')
    for legobj in leg.legendHandles:
        legobj.set_linewidth(4.0)
    _A.set_visible(False)
    _L.set_visible(False)
    ax.axhline(0.5,color='r',linewidth=3,linestyle='--')
    plt.tight_layout() 
    if fout: plt.savefig(fout)
    

def plot_adj_syn_mi(ax,data,fout=None):
    N = len(data['chem'])
    ax.plot(range(N),data['chem'],color=PRESYN_COL,
            linewidth=3,label='Chemical')
    N = len(data['gap'])
    ax.plot(range(N),data['gap'],color=GAPJUNC_COL,
            linewidth=3,label='Gap Junc.')    
    ax.set_ylabel('Mutual Info. (bits)',fontsize=32)
    ax.set_xlabel('k',fontsize=32)
    ax.set_ylim([0,0.5])
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_title('Contact surface area vs. synaptic volume',fontsize=32)
    plt.legend(loc='upper right',fontsize=24)
    plt.tight_layout()

def plot_adj_syn_weight_mi(ax,krange,data,fout=None):
    N = len(data['chem'])+1
    ax.plot(krange,data['chem'],color=PRESYN_COL,
            linewidth=3,label='Chemical')
   # N = len(data['gap'])+1
   # ax.plot(krange,data['gap'],color=GAPJUNC_COL,
   #         linewidth=3,label='Gap Junc.')    
    ax.set_ylabel('Mutual Info. (bits)',fontsize=32)
    ax.set_xlabel('k',fontsize=32)
    #ax.set_ylim([0,2])
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_title('Contact surface area vs. synaptic volume',fontsize=32)
    plt.legend(loc='upper right',fontsize=24)
    plt.tight_layout()

def plot_venn(ax,subsets,labels,title=None,color=None):
    v = venn2(subsets=subsets,ax=ax,set_labels=labels)
    if color:
        v.get_patch_by_id('100').set_color(color[0])
        v.get_patch_by_id('010').set_color(color[1])
        v.get_patch_by_id('11').set_color(color[2])
    for text in v.set_labels: text.set_fontsize(24)
    for text in v.subset_labels: text.set_fontsize(14)   
    if title: ax.set_title(title,fontsize=24)

def plot_corr(ax,l1,l2,label=None,_x=0.05,_y=0.9,corr_label='r',
              col='r',marker='o'):
    x1,x2 = [],[]
    for n in l1:
        if n != 'mean':
            x1.append(l1[n])
            x2.append(l2[n])

    r = pearsonr(x1,x2)
    
    ax.plot(x1,x2,color=col,marker=marker,linestyle='',label=label,
            markersize=9)
    ax.plot([0,1],[0,1],'k')
    ax.text(_x,_y,'%s = %1.3f' %(corr_label,r[0]**2),
            verticalalignment='bottom',
            transform=ax.transAxes,fontsize=24) 
    
    
def stars(p):
   if p < 0.0001:
       return "****"
   elif (p < 0.001):
       return "***"
   elif (p < 0.01):
       return "**"
   elif (p < 0.05):
       return "*"
   else:
       return "ns"
