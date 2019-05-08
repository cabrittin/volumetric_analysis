"""
mk_rnaseq_exp_matrix.py

Outputs the expression matrix for rnaseq data

"""

import sys
sys.path.append('./volumetric_analysis')
import argparse

import db
import aux

FOUT = 'rnaseq_exp_matrix'

def get_time_points(cur):
    sql = "select idtime,timename from time"
    cur.execute(sql)
    return [t for t in cur.fetchall()]

def get_gene_id(cur,gene):
    sql = "select idgenes from genes where genesname = '%s'" %gene
    cur.execute(sql)
    #print(gene,cur.fetchone()[0])
    idgene = cur.fetchone()
    if idgene: 
        return idgene[0]
    else: 
        return None

def get_gene_expression(cur,idgene,timept,cilb95 = 0):
    sql = ("select cells.cellsname,expression.tpmadj from expression "
            "join cells on cells.idcells = expression.idcells "
            "where expression.idtime = %d "
            "and expression.idgenes = '%s' "
            "and expression.cilb95 > %d "
            %(timept,idgene,cilb95) )
    cur.execute(sql)
    return cur.fetchall()

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('db',
                        action="store",
                        help="Database name")
    
    parser.add_argument('group',
                        action="store",
                        help="Path to group file")

    parser.add_argument('genes',
                        action = "store",
                        help="Path to cam gene list")

    parser.add_argument('dout',
                        action = "store",
                        help="Path to output director")

    params = parser.parse_args()


    con = db.connect.default(params.db)
    cur = con.cursor()
    tpts = get_time_points(cur)

    group = aux.read.into_dict2(params.group)
    gmap = aux.read.into_map(params.group)

    print(group)
    genes = aux.read.into_list(params.genes)
    genemap = {}
    rec = []
    for g in genes:
        idgene = get_gene_id(cur,g)
        if not idgene:
            print(g)
        else:
            genemap[g]=idgene
    for (idtime,tpt) in tpts:
        #print(tpt)
        matrix = []
        for g in genemap:
            #print('\t'+g)
            cells = get_gene_expression(cur,genemap[g],idtime)
            for c in cells:
                if c[0] in group:
                    for _c in group[c[0]]:
                        #print('asdfasd',_c)
                        matrix.append([g,_c,c[1]])
                elif c[0] in gmap:
                    #print(rec.append(c[0]))
                    matrix.append([g,c[0],c[1]])
        if not matrix: continue
        fout = params.dout.replace('.csv','_' + tpt + '.csv')
        aux.write.from_list(fout,matrix)







