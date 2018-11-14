"""
varshney_hierarchy.py




"""

import matplotlib.pyplot as plt

from connectome.load import load_lite
import networks.hierarchy as hierarchy
from figures.hierarchy import plot_varshney
import aux

DB = 'N2U'
GROUP = './mat/n2u_class_group.txt'
REMOVE = './mat/n2u_remove.txt'
THRESH = 0.1
L0 = ['MUBODY','MUHEAD','anal','int','intestine','sph','um','vm']
col = {'S':'r','I':'b','M':'#D7DF01','Mu':'#8904B1'}
NCLASS = './mat/n2u_class_simmu.txt'

def run():
    _group = aux.read.into_map(GROUP)
    group = {'map':_group,'key':'class'}
    remove = aux.read.into_map(REMOVE)
    C = load_lite(DB,chemical=True,electrical=True,
                  group=group,remove=remove,
                  dataType='networkx')
    C.combine_chem_and_elec()
    H = hierarchy.hierarchy(C.D,L0,THRESH)
    print(H)
    nodes,xyz = hierarchy.varshney_modified(H,C.D)
    nclass = aux.read.into_dict(NCLASS)
    fig,ax = plt.subplots(1,1,figsize=(20,10))
    plot_varshney(ax,nodes,xyz,C.D.edges(),nclass,col)
    #plt.show()
    
if __name__ == '__main__':
    run()


    
