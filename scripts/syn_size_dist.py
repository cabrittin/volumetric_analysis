import sys
sys.path.append(r'./volumetric_analysis')
import numpy as np
import db
import matplotlib.pyplot as plt
import matplotlib as mpl
import aux
import argparse

mpl.rcParams['xtick.labelsize'] = 24
mpl.rcParams['ytick.labelsize'] = 24

N2U_NRTHRESH = 150
N2U_INPLANETHRESH = 3000

N930_NRTHRESH = 920
N930_INPLANETHRSH = 8000

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
    
def get_synapses(cur,stype,nrthresh,inplanethresh):
    syn = db.mine.order_synapses_by_section_number(cur,stype=stype)
    data = []
    in_plane = []
    nr_45 = []
    not_nr = []
    for s in syn:
        w = int(s[2])
        if s[6] < nrthresh:
            if s[8] < inplanethresh:
                in_plane.append(w)
                loc = 'NR in plane'
            else:
                nr_45.append(w)
                loc = 'NR 45 deg'
        else:
            not_nr.append(w)
            loc = 'Not NR'
        data.append(list(s) + [loc])
    #print(len(in_plane),np.mean(in_plane),np.std(in_plane))
    #print(len(nr_45),np.mean(nr_45),np.std(nr_45))
    #print(len(not_nr),np.mean(not_nr),np.std(not_nr))

    return [in_plane,nr_45,not_nr],data

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Outputs list of synapse positions")
    parser.add_argument('db',
                        action="store",
                        help="Database name"
                        )
   
    params = parser.parse_args()

    if params.db == 'N2U':
        nrthresh = N2U_NRTHRESH
        inplane_thresh = N2U_INPLANETHRESH
    elif params.db == 'n930':
        nrthresh = N930_NRTHRESH
        inplane_thresh = N930_INPLANETHRESH
    
    con = db.connect.default(params.db)
    cur = con.cursor()

    chem,csyn = get_synapses(cur,'chemical',nrthresh,inplane_thresh)
    elec,celec = get_synapses(cur,'electrical',nrthresh,inplane_thresh)

    syn = csyn + celec
    aux.write.from_list('results/%s_syn_list_with_loc.csv'%params.db,syn)

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
