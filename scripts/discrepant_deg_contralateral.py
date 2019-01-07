"""
discrepant_deg_contralateral.py

Fraction of connections that are discrepant. Connections are classified 
as left discrepant (found only on the left side), right discrepant 
(found only on the right side) or left/right reproducible (found on both 
sides of the animal). Numbers represent the number of total connections.

created: Christopher Brittin
date: 01 November 2018

"""
import sys
sys.path.append(r'./volumetric_analysis')
import matplotlib.pyplot as plt
import matplotlib as mpl
import networkx as nx
from matplotlib.patches import Patch

#Brittin modules
from connectome.load import from_db
import connectome.synspecificity as synspec
from networks.stats import *
from figures.stats import *
import aux

mpl.rcParams['xtick.labelsize'] = 32
mpl.rcParams['ytick.labelsize'] = 32

ADULT_COL = '#FFC300'
L4_COL = '#3380FF'
AL_COL = '#14FA29'

left_nodes = './mat/left_nodes.txt'
right_nodes = './mat/right_nodes.txt'
lr_dict = './mat/lr_dict.txt'

def get_discrepant_directed(G,lr,left,right):
    _left,_right,_both = 0,0,0
    _edges = []
    for (u,v) in G.edges():
        if (u not in left) and (u not in right): continue
        ur,vr = lr[u],lr[v]
        if G.has_edge(ur,vr):
            _both += 1
        else:
            _edges.append([u,v,G[u][v]['weight']])
            if u in left:
                _left += 1
            else:
                _right += 1
    return [_left,_right,_both,float(_left + _right + _both)],_edges
                
def get_discrepant_undirected(G,lr,left,right):
    _left,_right,_both = 0,0,0
    _edges = []
    for (u,v) in G.edges():
        if (((u not in left) and (u not in right)) or
            ((v not in left) and (v not in right))):
            continue
        ur,vr = lr[u],lr[v]
        if G.has_edge(ur,vr):
            _both += 1
        else:
            _edges.append([u,v,G[u][v]['weight']])
            if (u in left) and (v in left):
                _left += 1
            elif (u in right) and (v in right):
                _right += 1
            else:
                _left += 0.5
                _right += 0.5
    return [_left,_right,_both,_left + _right + _both],_edges   

def format_data(syn,gap,adj):
    data = np.zeros((3,3))
    data[0,0] = gap[0]/gap[3]
    data[0,1] = syn[0]/syn[3]
    data[0,2] = adj[0]/adj[3]
    data[1,0] = gap[1]/gap[3]
    data[1,1] = syn[1]/syn[3]
    data[1,2] = adj[1]/adj[3]
    data[2,0] = gap[2]/gap[3]
    data[2,1] = syn[2]/syn[3]
    data[2,2] = adj[2]/adj[3]
    return data

def run(fout=None,source_data=None):
    left = aux.read.into_list(left_nodes)
    right = aux.read.into_list(right_nodes)
    _lrd = aux.read.into_dict(lr_dict)
    lrd = {}
    for key,val in _lrd.items():
        lrd[key] = val
        lrd[val] = key

    _remove = ['VC01','VD01','VB01','VB02']

    n2u = from_db('N2U',adjacency=True,chemical=True,electrical=True,
                  dataType='networkx',remove=_remove)
    jsh = from_db('JSH',adjacency=True,chemical=True,electrical=True,
                  dataType='networkx',remove=_remove)


    nsyn,sedges = get_discrepant_directed(n2u.C,lrd,left,right)
    ngap,gedges = get_discrepant_undirected(n2u.E,lrd,left,right)
    nadj,aedges = get_discrepant_undirected(n2u.A,lrd,left,right)
    print(nsyn,ngap,nadj)
    if source_data:
        dsource = 'adult'
        fsplit = source_data.split('.')
        fchem = fsplit[0] + '_' + dsource + '_chem.' + fsplit[1]
        fgap = fsplit[0] + '_' + dsource + '_gap.' + fsplit[1]
        fadj = fsplit[0] + '_' + dsource + '_adj.' + fsplit[1]
        aux.write.from_list(fchem,sedges)
        aux.write.from_list(fgap,gedges)
        aux.write.from_list(fadj,aedges)
        
    jsyn,sedges = get_discrepant_directed(jsh.C,lrd,left,right)
    jgap,gedges = get_discrepant_undirected(jsh.E,lrd,left,right)
    jadj,aedges = get_discrepant_undirected(jsh.A,lrd,left,right)
    print(jsyn,jgap,jadj)
    print(nsyn,ngap,nadj)
    if source_data:
        dsource = 'l4'
        fsplit = source_data.split('.')
        fchem = fsplit[0] + '_' + dsource + '_chem.' + fsplit[1]
        fgap = fsplit[0] + '_' + dsource + '_gap.' + fsplit[1]
        fadj = fsplit[0] + '_' + dsource + '_adj.' + fsplit[1]
        aux.write.from_list(fchem,sedges)
        aux.write.from_list(fgap,gedges)
        aux.write.from_list(fadj,aedges)
    _n2u = format_data(nsyn,ngap,nadj)
    _jsh = format_data(jsyn,jgap,jadj)

    print(_n2u)
    width = 0.25
    dx = 0.1
    ind = np.arange(3)
    col = ['#E69F00','#56B4E9','#F0E442']

    fig,ax = plt.subplots(1,1,figsize=(12,10))
    ax.axvspan(-0.5,0.5,facecolor='#C3C3C3')
    ax.axvspan(0.5,1.5,facecolor='#D8D7D7')
    ax.axvspan(1.5,2.5,facecolor='#C3C3C3')

    ax.bar(ind-width/2-dx,_n2u[0,:],width,color=col[0],hatch='/')
    ax.bar(ind-width/2-dx,_n2u[1,:],width,color=col[1],bottom=_n2u[0,:],hatch='/')
    ax.bar(ind-width/2-dx,_n2u[2,:],width,color=col[2],bottom=_n2u[0,:] + _n2u[1,:])

    ax.bar(ind+width/2+dx,_jsh[0,:],width,color=col[0],hatch='/')
    ax.bar(ind+width/2+dx,_jsh[1,:],width,color=col[1],bottom=_jsh[0,:],hatch='/')
    ax.bar(ind+width/2+dx,_jsh[2,:],width,color=col[2],bottom=_jsh[0,:] + _jsh[1,:])

    ax.set_xticklabels(['Adult\n(n=%d)'%ngap[3],
                        'L4\n(n=%d)'%jgap[3],
                        'Adult\n(n=%d)'%nsyn[3],
                        'L4\n(n=%d)'%jsyn[3],
                        'Adult\n(n=%d)'%nadj[3],
                        'L4\n(n=%d)'%jadj[3]])
    ax.set_xticks([-0.25,0.25,0.75,1.25,1.75,2.25])
    plt.xticks(fontsize=22)
    #plt.xticks(rotation=45,fontsize=18)

    ax.set_yticklabels([0,0.25,0.5,0.75,1.])
    ax.set_yticks([0,0.25,0.5,0.75,1.])

    ax.set_ylim([0,1.3])
    ax.set_xlim([-0.5,2.5])

    ax.set_ylabel('Fraction of connections',fontsize=32)

    ax.text(-0.2,1.1,'Gap J.',fontsize=34)
    ax.text(0.65,1.1,'Chemical',fontsize=34)
    ax.text(1.65,1.1,'Adjacency',fontsize=34)

    legend_elements = [Patch(facecolor=col[0], edgecolor='k',
                             label='Left discrepant',hatch='//'),
                       Patch(facecolor=col[1], edgecolor='k',
                             label='Right discrepant',hatch='//'),
                       Patch(facecolor=col[2], edgecolor='k',
                             label='L/R reproducible')]
    ax.legend(handles=legend_elements, loc='upper center',ncol=3,fontsize=18)
    plt.tight_layout()
    if fout: plt.savefig(fout)
    plt.show()

if __name__ == '__main__':
    run()
    



