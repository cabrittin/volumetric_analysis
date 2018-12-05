import sys
sys.path.append(r'./volumetric_analysis')
import numpy as np
import db
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['xtick.labelsize'] = 24
mpl.rcParams['ytick.labelsize'] = 24

NRTHRESH = 150
INPLANETHRESH = 3000

def plot_hist(ax,data,nbins=25,hrange=(0,25),color='k',label=None):
    ax.hist(data,bins=nbins,range=hrange,histtype='step',
            density=True,cumulative=1,linewidth=3,label=label)

def plot_syn(ax,syn):
    plot_hist(ax,syn[0],color='r',label='NR in plane (n = %d)'%len(syn[0]))
    plot_hist(ax,syn[1],color='b',label='NR 45 deg (n = %d)'%len(syn[1]))
    plot_hist(ax,syn[2],color='g',label='Not NR (n = %d)'%len(syn[2]))
    ax.set_xlim([1,10])
    ax.set_ylim([0,1])   
    ax.legend(loc='lower right',fontsize=24)
    
def get_synapses(cur,stype):
    syn = db.mine.order_synapses_by_section_number(cur,stype=stype)
    in_plane = []
    nr_45 = []
    not_nr = []
    for s in syn:
        w = int(s[2])
        if s[6] < NRTHRESH:
            if s[8] < INPLANETHRESH:
                in_plane.append(w)
            else:
                nr_45.append(w)
        else:
            not_nr.append(w)
        
    #print(len(in_plane),np.mean(in_plane),np.std(in_plane))
    #print(len(nr_45),np.mean(nr_45),np.std(nr_45))
    #print(len(not_nr),np.mean(not_nr),np.std(not_nr))

    return [in_plane,nr_45,not_nr]

_db = 'N2U'
con = db.connect.default(_db)
cur = con.cursor()

chem = get_synapses(cur,'chemical')
elec = get_synapses(cur,'electrical')

fig,ax = plt.subplots(1,2,figsize=(15,10),sharey=True,sharex=True)
plot_syn(ax[0],chem)
plot_syn(ax[1],elec)
ax[0].set_ylabel('Cumulative distribution',fontsize=24)
ax[0].set_xlabel('Synapse weight (# EM sections)',fontsize=24)
ax[1].set_xlabel('Synapse weight (# EM sections)',fontsize=24)
ax[0].set_title('Chemical synapses',fontsize=24)
ax[1].set_title('Gap junctions',fontsize=24)
plt.savefig('results/nr_vs_nonnr_syn.png')
plt.show()
