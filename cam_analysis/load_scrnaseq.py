"""
load_scrnaseq.py

Loads data from GSE126954

See https://www.biorxiv.org/content/10.1101/565549v2

created: Christopher Brittin
date: 09 March 2019
"""
import sys
sys.path.append('./volumetric_analysis')
import argparse

import db

class Data:
    def __init__(self):
        self.dmap = {}
        self.idx = 0
    
    def add(self,x):
        if x not in self.dmap:
            self.dmap[x] = self.idx
            self.idx += 1

    def get_rows(self):
        return sorted([(y,x) for (x,y) in list(self.dmap.items())])

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('db',
                        action="store",
                        help="Database name")
    
    parser.add_argument('fin',
                        action="store",
                        help="Path the table of terminal expression")

    params = parser.parse_args()

    con = db.connect.default(params.db)
    cur = con.cursor()

    #data = aux.read.into_list2(params.fin,delimiter='\t')

    cells = Data()
    time = Data()
    genes = Data()
    genes.wb = {}
    exp = []
    idx = 0
    
    with open(params.fin,'r') as data:
        next(data)
        for _d in data:
            d = _d.split('\t')
            _gene = d[0]
            _wb = d[1]
            if ':' not in d[2]: continue
            [_cell,_time] = d[2].split(':')
            tpm_raw = d[3]
            tpm_adj = d[4]
            tpm_est = d[5]
            cil95 = d[6]
            cil80 = d[7]
            ciu80 = d[8]
            ciu95 = d[9].split('\n')[0]
            genes.add(_gene)
            genes.wb[_gene] = _wb
            idgene = genes.dmap[_gene]
            time.add(_time)
            idtime = time.dmap[_time]
            cells.add(_cell)
            idcell = cells.dmap[_cell]
            
            exp.append((idx,idgene,idcell,idtime,tpm_raw,tpm_adj,tpm_est,cil95,cil80,ciu80,ciu95))
            idx += 1

    cur.execute("SET FOREIGN_KEY_CHECKS=0")
    
    sql = "INSERT INTO cells (idcells,cellsname) values (%s,%s)"
    cur.executemany(sql,cells.get_rows())

    sql = "INSERT INTO time (idtime,timename) values (%s,%s)"
    cur.executemany(sql,time.get_rows())

    gwb = [(i,g,genes.wb[g]) for (i,g) in genes.get_rows()]
    sql = "INSERT INTO genes (idgenes,genesname,wbname) values (%s,%s,%s)"
    cur.executemany(sql,gwb)

    sql = ("INSERT INTO expression "
            "(idexpression,idgenes,idcells,idtime,tpmraw,tpmadj,tpmmed,cilb95,cilb80,ciub80,ciub95) " 
            "values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)")
    cur.executemany(sql,exp)

    con.commit()
