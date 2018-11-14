"""
hierarchy.py

Module for plotting network hierarchy

Created: Christopher Brititn
date: 14 November 2018

"""
import matplotlib.pyplot as plt

def plot_varshney(ax,nodes,xyz,edges,nclass,col):
    N = len(nodes)
    nk = dict([(nodes[i],i) for i in range(N)])
    for i in range(N):
        ax.text(xyz[i,3],xyz[i,2],nodes[i],
                color=col[nclass[nodes[i]]],
                fontsize=14)

    for (i,j) in edges:
        zi,xi = int(xyz[nk[i],2]),xyz[nk[i],3]
        zj,xj = int(xyz[nk[j],2]),xyz[nk[j],3]
        col = 'k'
        if zi < zj: 
            col = '#F7819F' #Up the hierarchy
        elif zi == zj: 
            col = '#BDBDBD' #Same hierarchy
        elif zi > zj: 
            col = '#9FF781' #Down the hierarchy
        #print i,j,zi,zj,col
        plt.plot([xi,xj],[zi,zj],color=col,
                 alpha=0.2)

    plt.xlim([-1.1,1.1])
    plt.ylim([0,10])
    plt.ylabel('Hierarchy',fontsize=24)
    plt.show()
