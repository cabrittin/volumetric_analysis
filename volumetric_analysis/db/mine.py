"""
mine.py
Author: Christopher Brittin

Functions for mining Elegance database

"""


import aux

def get_db(cur):
    cur.execute("select database()")
    return cur.fetchone()[0]

def get_neurons(cur):
    sql = ("select distinct(CON_AlternateName) "
           "from contin "
           "where type like 'neuron' "
           "and CON_Remarks like '%OK%'")
    cur.execute(sql)
    
    return list(set([aux.format.rm_brack(a[0]) for a in cur.fetchall()]))

def get_img_number(cur,**kwargs):
    args = aux.format.get_args(kwargs,{'start':-1,'end':1e6})
    sql = ("select IMG_Number "
           "from image "
           "where IMG_Series in ('NR','VC') " 
           "and IMG_SectionNumber >= %d "
           "and IMG_SectionNumber <= %d" %(int(args['start']),int(args['end'])))
    cur.execute(sql)
    return [i[0] for i in cur.fetchall()]

def get_obj_area(cur):
    sql = ("select object,imgNum,area "
           "from dimensions ")
    cur.execute(sql)
    return dict([(a[0],{'image':a[1],'area':a[2]}) for a in cur.fetchall()])
           

def get_synapse_contins(cur,stype,**kwargs):
    args = aux.format.get_args(kwargs,{'start':None,'end':None,'syn_list':None})
    if args['start'] or args['end']: 
        images = get_img_number(cur,start=args['start'],end=args['end'])
        images = ''.join(["'","','".join(images),"'"])
        sql = ("select synapsecombined.continNum "
               "from synapsecombined "
               "join object on object.CON_Number = synapsecombined.continNum "
               "where synapsecombined.type like '%s' "
               "and object.IMG_Number in (%s) " %(stype,images))
    else:
        sql = ("select synapsecombined.continNum "
               "from synapsecombined "
               "join object on object.CON_Number = synapsecombined.continNum "
               "where synapsecombined.type like '%s' " %(stype))			
    cur.execute(sql)
    return ','.join([str(a[0]) for a in cur.fetchall()])   

def get_synapse_data(cur,stype,**kwargs):		
    contins = get_synapse_contins(cur,stype,**kwargs)
    sql = ("select pre,post,sections,continNum,series "
           "from synapsecombined "
           "where continNum in (%s) " %(contins))	
    cur.execute(sql)
    return cur.fetchall()
    
def get_adjacency_data(cur):
    sql = ("select pre,post,weight,imgNum "
           "from adjacency2")
    cur.execute(sql)
    return [(a[0],a[1],int(a[2]),a[3]) for a in cur.fetchall()]

def get_adjacency_cells(cur):
    sql = ("select distinct(pre) from adjacency2 "
           "union "
           "select distinct(post) from adjacency2")
    cur.execute(sql)
    return [a[0] for a in cur.fetchall()]


def get_touch_density(cur,key="pixels_norm"):
    sql = ("select pre,post,%s "
           "from touch_density" %key)
    cur.execute(sql)
    return [(a[0],a[1],float(a[2])) for a in cur.fetchall()]

def maxAnterior(cur,cell):
    sql = ("select min(image.IMG_SectionNumber) "
           "from radialPharynx "
           "join object on object.OBJ_Name=radialPharynx.OBJ_Name "
           "join contin on contin.CON_Number=object.CON_Number "
           "join image on image.IMG_Number=object.IMG_Number "
           "where contin.CON_AlternateName like '%%%s%%' "
           %cell
           )
    
    cur.execute(sql)
    return cur.fetchone()[0]

def maxPosterior(cur,cell):
    sql = ("select max(image.IMG_SectionNumber) "
           "from radialPharynx "
           "join object on object.OBJ_Name=radialPharynx.OBJ_Name "
           "join contin on contin.CON_Number=object.CON_Number "
           "join image on image.IMG_Number=object.IMG_Number "
           "where contin.CON_AlternateName like '%%%s%%' "
           %cell
           )
    
    cur.execute(sql)
    return cur.fetchone()[0]

def neuron_cylinder(cur):
    """
    Gets cylindrical coordinates of neuron locations in db
    """
    sql = ("select contin.CON_AlternateName,"
           "radialPharynx.distance,"
           "radialPharynx.phi,"
           "image.IMG_SectionNumber "
           "from radialPharynx "
           "join object on object.OBJ_Name=radialPharynx.OBJ_Name "
           "join contin on contin.CON_Number=object.CON_Number "
           "join image on image.IMG_Number=object.IMG_Number "
           "where contin.type like 'neuron'"
           )
    cur.execute(sql)
    data = [[a[0],int(a[1]),float(a[2]),int(a[3])] for a in cur.fetchall()]
    return data     
    
def gap_cylinder(cur,dr=0,dphi=0.1):
    """
    Gets cylindrical coordinates of synapses in db
    """
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
        
    return gap    

def syn_cylinder(cur):
    """
    Gets cylindrical coordinates of synapses in db
    """
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
