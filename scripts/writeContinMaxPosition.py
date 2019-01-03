"""
writeContinMaxPosition.py

Write the max anterior and posterior positions for cell in the nerve ring, 
broken down by contin.

created: Christopher Brittin
date: 03 January 2019

"""

import argparse
import progressbar

import db
import aux

neuron_class = './mat/nerve_ring_classes.txt'

DENDRITE = ['Sp1','Sp2']

def run(_db,fout,maxpos=200,scale=1):
    con = db.connect.default(_db)
    cur = con.cursor()
    nclass = aux.read.into_dict(neuron_class)

    bar = progressbar.ProgressBar(maxval=len(nclass),
                                  widgets=[progressbar.Bar('.','[',']'),
                                           'Max position',
                                           progressbar.Percentage()])    
    data = []
    idx = 0
    bar.start()
    for n in sorted(nclass):
        contins = db.mine.get_contins(cur,n)
        for c in contins:
            apos = db.mine.maxAnterior(cur,n,contin=c)
            ppos = db.mine.maxPosterior(cur,n,contin=c)
            if not apos: apos = -1
            if not ppos: ppos = -1
            data.append([n,nclass[n],c,apos,ppos])
        idx += 1
        bar.update(idx)

    aux.write.from_list(fout,data)        
    con.close()
    
if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('db',
                        action="store",
                        help="Database name")
    parser.add_argument('fout',
                        action="store",
                        help="Path to output file")

    params = parser.parse_args()
    if params.db == 'N2U':
        maxpos = 182
    elif params.db == 'JSH':
        maxpos = 200
    
    run(params.db,params.fout,maxpos=maxpos)
    print("")
