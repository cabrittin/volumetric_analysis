"""
spatialMap.py

Submodule for plotting spatial locations of neuron locations and synapses.

Author: Christopher Brittin
Created: 07 February 2018

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from random import shuffle

#Brittin Modules
import aux
import db

mpl.rcParams['xtick.labelsize'] = 24
mpl.rcParams['ytick.labelsize'] = 24
drad = 0.2



def plot_syn_phiz(db,pre,post,gap,display,nclass,ccode,left,right,fout=None):
    display = display.split(',')
    nclass = aux.read.into_dict(nclass)
    ccode = aux.read.into_dict(ccode)
    left = aux.read.into_list(left)
    right = aux.read.into_list(right)
    for c in ccode: ccode[c] = ccode[c].lower()

    #loc = neuron_cylinder(db)
    #shuffle(loc)
    
    #print(len(loc),len(pre),len(post),len(gap))
    
    #loc = scale_data(loc)
    pre = scale_data(pre)
    post = scale_data(post)
    gap = scale_data(gap)
    
    #loc = split_left_right(loc,left,right)
    pre = split_left_right(db,pre,left,right)
    post = split_left_right(db,post,left,right)
    gap = split_left_right(db,gap,left,right)
    if display[0] not in ['All','all','ALL']:
        #loc = parse_data(loc,nclass,display)
        pre = parse_data(pre,nclass,display)
        post = parse_data(post,nclass,display)
        gap = parse_data(gap,nclass,display)
        
    fs = 32     
    width,height = 10,10
    bkgcol = '#EAEAEA'
    #fig = plt.figure(figsize=(25,15))
    #ax=plt.subplot(231)
    fig,ax = plt.subplots(1,1,figsize=(width,height))
    ax.patch.set_facecolor(bkgcol)
    plot_rz(ax,pre,nclass,ccode,True)
    for x in [20,90,150]: ax.axhline(x,c='r',lw=2,linestyle='--')
    ax.text(8300,25,'EM I',ha='right',fontsize=fs)
    ax.text(8300,91,'EM II',ha='right',fontsize=fs)
    ax.text(8300,151,'EM III',ha='right',fontsize=fs)
    #ax.set_xlim([-8000,8000])
    ax.set_xticks([-8000,-4000,0,4000,8000])
    ax.set_xticklabels([-8,-4,0,4,8])
    ax.set_title('Presynaptic neurons',fontsize=fs)
    ax.set_ylabel('<-Anterior/Posterior->, $\mu$m',fontsize=fs)
    ax.set_xlabel('Radial distance, $\mu$m',fontsize=fs)
    plt.tight_layout()
    if fout:
        tmp = fout.replace('.png','_i.png')
        plt.savefig(tmp)
    
    fig,ax=plt.subplots(1,1,figsize=(width,height))
    ax.patch.set_facecolor(bkgcol)
    plot_phiz(ax,pre,nclass,ccode,True)
    #ax.set_title('Neuron location',fontsize=24)
    for x in [20,90,150]: ax.axhline(x,c='r',lw=2,linestyle='--')
    #ax.text(8300,21,'EM I',fontsize=24,ha='right')
    #ax.text(8300,91,'EM II',fontsize=24,ha='right')
    #ax.text(8300,136,'EM III',fontsize=24,ha='right')    
    ax.set_ylabel('<-Anterior/Posterior->, $\mu$m',fontsize=fs)
    ax.set_xlabel('Azimuth angle, $\phi$',fontsize=fs)
    #ax.set_title('Azimuth angle',fontsize=28)
    plt.tight_layout()
    if fout:
        tmp = fout.replace('.png','_ii.png')
        plt.savefig(tmp)
    
    
    fig,ax = plt.subplots(1,1,figsize=(width,height))
    ax.patch.set_facecolor(bkgcol)
    plot_rz(ax,post,nclass,ccode,False)
    #ax.set_ylabel('<-Anterior/Posterior->, $\mu$m',fontsize=28)
    ax.set_xticks([-8000,-4000,0,4000,8000])
    ax.set_xticklabels([-8,-4,0,4,8])
    ax.set_title('Postsynaptic neurons',fontsize=fs)
    ax.set_xlabel('Radial distance, $\mu$m',fontsize=fs)
    plt.tight_layout()
    if fout:
        tmp = fout.replace('.png','_iii.png')
        plt.savefig(tmp)

    fig,ax = plt.subplots(1,1,figsize=(width,height))
    ax.patch.set_facecolor(bkgcol)
    plot_phiz(ax,post,nclass,ccode,False)
    #ax.set_title('Synaptic output',fontsize=24)
    ax.set_xlabel('Azimuth angle, $\phi$',fontsize=fs)
    plt.tight_layout()
    if fout:
        tmp = fout.replace('.png','_iv.png')
        plt.savefig(tmp)
        
    
    fig,ax=plt.subplots(1,1,figsize=(width,height))
    ax.patch.set_facecolor(bkgcol)
    plot_rz(ax,gap,nclass,ccode,False)
    #ax.set_ylabel('<-Anterior/Posterior->, $\mu$m',fontsize=28)
    ax.set_xticks([-8000,-4000,0,4000,8000])
    ax.set_xticklabels([-8,-4,0,4,8])
    ax.set_title('Gap junctions',fontsize=fs)
    ax.set_xlabel('Radial distance, $\mu$m',fontsize=fs)
    plt.tight_layout()
    if fout:
        tmp = fout.replace('.png','_v.png')
        plt.savefig(tmp)       

    fig,ax = plt.subplots(1,1,figsize=(width,height))
    ax.patch.set_facecolor(bkgcol)
    plot_phiz(ax,gap,nclass,ccode,False)
    #ax.set_title('Synaptic input',fontsize=fs)
    ax.set_xlabel('Azimuth angle, $\phi$',fontsize=fs)
    plt.tight_layout()
    if fout:
        tmp = fout.replace('.png','_vi.png')
        plt.savefig(tmp)
    """
    fig,ax=plt.subplots(1,1,figsize=(width,height))
    ax.patch.set_facecolor(bkgcol)
    plot_rz(ax,gap,nclass,ccode,False)
    ax.set_ylabel('<-Anterior/Posterior->, $\mu$m',fontsize=28)
    ax.set_xlim([-8000,8000])
    ax.set_xticks([-8000,-4000,0,4000,8000])
    ax.set_xticklabels([-8,-4,0,4,8])
    #ax.set_title('Synaptic input',fontsize=28)
    ax.set_xlabel('Radial distance, $\mu$m',fontsize=fs)
    plt.tight_layout()
    if fout:
        tmp = fout.replace('.png','_vii.png')
        plt.savefig(tmp)    

    fig,ax = plt.subplots(1,1,figsize=(width,height))
    ax.patch.set_facecolor(bkgcol)
    plot_phiz(ax,gap,nclass,ccode,False)
    #ax.set_title('Synaptic output',fontsize=24)
    ax.set_xlabel('Azimuth angle, $\phi$',fontsize=fs)
    plt.tight_layout()
    if fout:
        tmp = fout.replace('.png','_viii.png')
        plt.savefig(tmp)   

    """
        
    #plt.show()


def parse_data(data,ncat,display):
    data2 = []
    for d in data:
        if ncat[d[0]] in display:
            data2.append(d)
    return data2    

def split_left_right(db,data,left,right):
    data2 = []
    if db == 'N2U':
        maxd = max([d[1] for d in data])
        
    for d in data:
        if d[0] in left:
            if db == 'N2U':
                d[1] = d[1] - maxd
            elif db == 'JSH':
                d[1] = -1*d[1]
            d[2] -= np.pi - drad
            #d[2] = -1
            data2.append(d)
        elif d[0] in right:
            d[2] += np.pi + drad
            #d[2] = 1
            data2.append(d)
    return data2
    
def scale_data(data):
    #Conversts pixes to nanometers
    data2 = []
    scale = 5
    for d in data:
        d[1] = float(d[1])*scale
        data2.append(d)
    return data2


def plot_rz(ax,data,ncat,col,yticks):
    #ax.set_xticks([-8000,-6000,-4000,-2000,0,2000,4000,6000,8000])
    #ax.set_xticklabels([8000,6000,4000,2000,0,2000,4000,6000,8000])

    if yticks:
        ax.set_yticks(range(181)[::20])
        ylabels = range(181)[::20]
        ylabels = [0.09*y for y in ylabels]
        ylabels = ['%2.1f'%y for y in ylabels]
        ax.set_yticklabels(ylabels)
    else:
        ax.set_yticklabels([])        

    for d in data:
        if d[0] in ncat:
            ax.plot(d[1],d[3],marker='s',markersize=6,
                    mec=col[ncat[d[0]]],mfc=col[ncat[d[0]]])
    
    ax.axvline(0,c='k',lw=2)
    ax.text(0.02,0.16,'Left neurons',ha='left',
            transform=ax.transAxes,fontsize=24)
    ax.text(0.99,0.16,'Right neurons',ha='right',
            transform=ax.transAxes,fontsize=24)
    
 
def plot_phiz(ax,data,ncat,col,yticks):
    ax.set_xticks([-2*np.pi-drad,-np.pi-drad,-drad,
                    drad,np.pi+drad,2*np.pi+drad])
    ax.set_xticklabels(['-$\pi$','0','$\pi$',
                        '-$\pi$','0','$\pi$'])
    
    if yticks:
        ax.set_yticks(range(181)[::20])
        ylabels = range(181)[::20]
        ylabels = [0.09*y for y in ylabels]
        ylabels = ['%2.1f'%y for y in ylabels]
        ax.set_yticklabels(ylabels)
    else:
        ax.set_yticklabels([])
    
    for d in data:
        if d[0] in ncat:
            ax.plot(d[2],d[3],marker='s',markersize=6,
                    mec=col[ncat[d[0]]],mfc=col[ncat[d[0]]])
    
    ax.axvline(0,c='k',lw=2)
    ax.text(0.02,0.16,'Left neurons',ha='left',
            transform=ax.transAxes,fontsize=24)
    ax.text(0.99,0.16,'Right neurons',ha='right',
            transform=ax.transAxes,fontsize=24)
    ax.set_xlim([-2.1*np.pi,2.1*np.pi])   
