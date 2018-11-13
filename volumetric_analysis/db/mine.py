"""
mine.py

Module for extracting and writing data to the database

Functions:
---------
get_db(cur)
   Returns database name

get_neurons(cur)
   Returns list of all cell names

get_img_number(cur,**kwargs)
   Returns list of image(sections) names from the NR an VC
   series. Use kwargs to define range = (start,end) where 
   start and end are ints.

get_obj_area(cur)
   Returns list of objects and associated areas

get_synapse_contins(cur,stype,**kwargs)
   Returns list of synapse contin numbers

get_synapse_data(cur,stype,**kwargs)
   Returns list of synapses

get_adjacency_data(cur)
   Returns list of adjacency data

get_adjacency_cells(cur)
   Return list of cell in the adjacency data

get_touch_density(cur,key="pixels_norm")
   Returns list of 'touch_density' measures. Possible key
   values are [pixels,pixels_norm,segments,segments_norm]

maxAnterior(cur,cell)
   Returns the max anterior section of cell

maxPosterior(cur,cell)
   Returns the max posterior sections of cell

neuron_cylinder(cur)
   Returns the cylindrical locations of neuron segment centroids

gap_cylinder(cur)
   Return the cylindrical coordinates of the gap junction

syn_cylinder(cur)
   Return the cylindrical coordinates of synapses

get_cell_genes(cur,cell):
    Return the genes expressed in cell

created: Christopher Brittin
date: 01 November 2018
"""


import aux

def get_db(cur):
    """
    Returns database name

    Parameters:
    -----------
    cur : MySQLdb cursors

    """
    cur.execute("select database()")
    return cur.fetchone()[0]

def get_neurons(cur):
    """
    Returns list of all contin names

    Parameters:
    -----------
    cur : MySQLdb cursors

    """
    sql = ("select distinct(CON_AlternateName) "
           "from contin "
           "where type like 'neuron' "
           "and CON_Remarks like '%OK%'")
    cur.execute(sql)
    
    return list(set([aux.format.rm_brack(a[0]) for a in cur.fetchall()]))

def get_img_number(cur,**kwargs):
    """
    Returns list of image(sections) names from the NR an VC series

    Parameters:
    cur : MySQLdb cursors
    start : int
       First section number (default -1)
    end : int
       Last section number (default 1e6)
    """   
    args = aux.format.get_args(kwargs,{'start':-1,'end':1e6})
    sql = ("select IMG_Number "
           "from image "
           "where IMG_Series in ('NR','VC') " 
           "and IMG_SectionNumber >= %d "
           "and IMG_SectionNumber <= %d" %(int(args['start']),int(args['end'])))
    cur.execute(sql)
    return [i[0] for i in cur.fetchall()]

def get_obj_area(cur):
    """
    Returns dictionary of (key,value) = (object_number,[image_number,area])

    Parameters:
    -----------
    cur : MySQLdb cursor

    """
    sql = ("select object,imgNum,area "
           "from dimensions ")
    cur.execute(sql)
    return dict([(a[0],{'image':a[1],'area':a[2]}) for a in cur.fetchall()])
           

def get_synapse_contins(cur,stype,**kwargs):
    """
    Returns list of synapse contin numbers

    Parameters:
    -----------
    cur : MySQLdb cursors
    stype : str
      synapses type; ['chemical','electrical']
    start : int
      beginning section number (default -1)
    end : int
      last section number (default 1e6)

    """
    args = aux.format.get_args(kwargs,{'start':None,'end':None,})
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
    """
    Returns list of synapses
    [pre_cell,post_cell,number_of_sections,continNum,series]
    
    Parameters:
    -----------
    cur : MySQLdb cursor
    stype : str
      synapse type; ['chemical','electrical']
    start : int
      beginning section number (default -1)
    end : int
      last section number (default 1e6)    
    
    """
    contins = get_synapse_contins(cur,stype,**kwargs)
    sql = ("select pre,post,sections,continNum,series "
           "from synapsecombined "
           "where continNum in (%s) " %(contins))	
    cur.execute(sql)
    return cur.fetchall()
    
def get_adjacency_data(cur):
    """
    Return list of adjacency data,
    [cell1,cell2,amount_of_contact,image_number]

    Parameters:
    -----------
    cur : MySQLdb cursor

    """
    sql = ("select pre,post,weight,imgNum "
           "from adjacency2")
    cur.execute(sql)
    return [(a[0],a[1],int(a[2]),a[3]) for a in cur.fetchall()]

def get_adjacency_cells(cur):
    """
    Returns list of cells in the adjacency data

    Parameters:
    -----------
    cur : MySQLdb cursor

    """
    sql = ("select distinct(pre) from adjacency2 "
           "union "
           "select distinct(post) from adjacency2")
    cur.execute(sql)
    return [a[0] for a in cur.fetchall()]


def get_touch_density(cur,key="pixels_norm"):
    """
    Returns list of 'touch_density' measures;
    [cell1,cell2,touch_density]

    Parameters:
    -----------
    cur : MySQLdb cursor
    key : str
     touch density measure [pixels,pixels_norm,segments,segments_norm]
    
    """
    sql = ("select pre,post,%s "
           "from touch_density" %key)
    cur.execute(sql)
    return [(a[0],a[1],float(a[2])) for a in cur.fetchall()]

def maxAnterior(cur,cell):
    """
    Returns max anterior section of cell

    Parameters:
    -----------
    cur : MySQLdb cursor
    cell : str
      Name of cell
    """
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
    """
    Returns the max posterior sections of cell

    Parameters:
    -----------
    cur : MySQLdb cursor
    cell : str
      name of cell
    """
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
    Returns the cylindrical locations of neuron segment centroids
    [cell_name,radius,phi,z]

    Parameters:
    ----------
    cur : MySQLdb
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
    Return the cylindrical coordinates of the gap junction
    [pre_cell,post_cell,radius,phi,z]

    Parameters:
    -----------
    cur : MySQLdb cursor
    dr : int
      shift for radius values (default 0)
    dphi : float
      shift for phi values (default 0.1)

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
    Return the cylindrical coordinates of synapses
    [pre_cell,post_cell,radius,phi,z]

    Parameters:
    -----------
    cur : MySQLdb cursor

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

def get_cell_genes(cur,cell):
    """
    Return the genes expressed in cell

    Parameters
    ----------
    cur : MySQLdb cursor
    cell : str
      Name of cell
    """
    sql = ("select genename from WB.association "
           "join WB.cells on "
           "WB.association.idcell = WB.cells.idcells "
           "where WB.cells.name like '%s'"
           %cell)
    cur.execute(sql)
    return [a[0] for a in cur.fetchall()]
