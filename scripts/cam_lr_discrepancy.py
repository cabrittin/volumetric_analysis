"""
cam_lr_discrepancy.py

Plots fraction of cells that synapse onto a left(right) cell while being also being a neighbor to the right(left) cell

"""
import sys
sys.path.append('./volumetric_analysis')
sys.path.append('.')
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 24
mpl.rcParams['ytick.labelsize'] = 24


from mat_loader import MatLoader

FOUT = '/home/cabrittin/Dropbox/PhD/sr_vol/figs2/fig9/cam_lr_discrepancy.png' 

if __name__=='__main__':
    M = MatLoader()
    C = M.load_consensus_graphs(4)
    M.load_left()
    M.load_lrmap()

    count = [0,0,0,0]
    for n in M.left:
        if C.C.has_node(n):
            count[0] += 1
            partners = set(C.C.neighbors(n))
            neighbors = set(C.A.neighbors(n))
            nonsyn = neighbors - partners
            for p in partners:
                pmap = M.lrmap[p]
                if pmap in nonsyn:
                    count[1] += 1
                    break
        
        if C.E.has_node(n):
            count[2] += 1
            partners = set(C.E.neighbors(n))
            nonsyn = neighbors - partners
            for p in partners:
                pmap = M.lrmap[p]
                if pmap in nonsyn:
                    count[3] += 1
                    break


    print(count)
    
    cfrac = float(count[1]) / count[0]
    gfrac = float(count[3]) / count[2] 

    data = [cfrac,gfrac]

    fig,ax = plt.subplots(1,1,figsize=(9,10))
    ind = [1,2]
    width = 0.3

    rects = ax.bar(ind,data,width,color='k')
    ax.set_ylabel('Fraction of cell classes that make synapse',fontsize=28)
    ax.set_xticks(ind)
    ax.set_xticklabels(('Chemical\n($n$=%d)'%count[0],'Gap J.\n($n$=%d)'%count[1]))
    ax.set_ylim([0,1])
    ax.set_title('Fraction of cell classes that synapse onto\n left(right) but not right(left) neighbor\n',fontsize=28)
    plt.tight_layout()
    plt.savefig(FOUT)
    plt.show()

    
 
