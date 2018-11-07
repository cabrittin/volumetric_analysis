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
import DB
import LUT

mpl.rcParams['xtick.labelsize'] = 24
mpl.rcParams['ytick.labelsize'] = 24
drad = 0.2



def plot_syn_phiz(db,display,nclass,ccode,fout=None):
    display = display.split(',')
    nclass = aux.read.into_dict(nclass)
    ccode = aux.read.into_dict(ccode)    
    for c in ccode: ccode[c] = ccode[c].lower()

    pre,post = syn_cylinder(db)
    gap = gap_cylinder(db)
    loc = neuron_cylinder(db)
    shuffle(loc)
    
    print(len(loc),len(pre),len(post),len(gap))
    
    loc = scale_data(loc)
    pre = scale_data(pre)
    post = scale_data(post)
    gap = scale_data(gap)
    
    loc = split_left_right(db,loc)
    pre = split_left_right(db,pre)
    post = split_left_right(db,post)
    gap = split_left_right(db,gap)
    if display[0] not in ['All','all','ALL']:
        loc = parse_data(loc,nclass,display)
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


        
    #plt.show()


def parse_data(data,ncat,display):
    data2 = []
    for d in data:
        if ncat[d[0]] in display:
            data2.append(d)
    return data2    

def syn_cylinder(db):
    """
    Gets cylindrical coordinates of synapses in db
    """
    con = DB.connect.default(db)
    cur = con.cursor()
    sql = ("select synapsecombined.pre,"
           "synapsecombined.post,"
           "radialPharynx.distance,"
           "radialPharynx.phi,"
           "image.IMG_SectionNumber "
           "from synapsecombined "
           "join radialPharynx on radialPharynx.OBJ_Name = synapsecombined.mid "
           "join object on object.OBJ_Name = radialPharynx.OBJ_Name "
           "join image on image.IMG_Number=object.IMG_Number "
           "where synapsecombined.type='chemical'"
           )

    cur.execute(sql)
    pre,post = [],[]
    for a in cur.fetchall():
        r = int(a[2])
        phi = float(a[3])
        z = int(a[4])
        pre.append([a[0],r,phi,z])
        for p in a[1].split(','):
            post.append([p,r,phi,z])
    return pre,post

def gap_cylinder(db,dr=0,dphi=0.1):
    """
    Gets cylindrical coordinates of synapses in db
    """
    con = DB.connect.default(db)
    cur = con.cursor()
    sql = ("select synapsecombined.pre,"
           "synapsecombined.post,"
           "radialPharynx.distance,"
           "radialPharynx.phi,"
           "image.IMG_SectionNumber "
           "from synapsecombined "
           "join radialPharynx on radialPharynx.OBJ_Name = synapsecombined.mid "
           "join object on object.OBJ_Name = radialPharynx.OBJ_Name "
           "join image on image.IMG_Number=object.IMG_Number "
           "where synapsecombined.type='electrical'"
           )

    cur.execute(sql)
    gap = []
    for a in cur.fetchall():
        r = int(a[2])
        phi = float(a[3])
        z = int(a[4])
        gap.append([a[0],r-dr,phi-dphi,z])
        #gap.append([a[1],r+dr,phi+dphi,z])
        
    return gap


def neuron_cylinder(db):
    """
    Gets cylindrical coordinates of neuron locations in db
    """
    con = DB.connect.default(db)
    cur = con.cursor()
    sql = ("select contin.CON_AlternateName,"
           "radialPharynx.distance,"
           "radialPharynx.phi,"
           "image.IMG_SectionNumber "
           "from radialPharynx "
           "join object on object.OBJ_Name=radialPharynx.OBJ_Name "
           "join contin on contin.CON_Number=object.CON_Number "
           "join image on image.IMG_Number=object.IMG_Number "
           )
    cur.execute(sql)
    data = [[a[0],int(a[1]),float(a[2]),int(a[3])] for a in cur.fetchall()]
    return data      

def split_left_right(db,data):
    lr = LUT.lut.neuron_class(db)
    left = aux.read.into_list(lr['left_nodes'])
    right = aux.read.into_list(lr['right_nodes'])

    data2 = []
    for d in data:
        if d[0] in left:
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
