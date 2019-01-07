"""
writeMaxPosition.py

Writes the max anteriro and posterior positions in the nerve ring. Creates
source data files to be used as input for scripts/max_anterior.py.

created: Christopher Brittin
date: 02 January 2019 

"""

import argparse
import progressbar

#Brittin modules
import db
import aux

MAXPOSITION = 200
SCALE = 80/1000.
SCREEN = ['SABD']

neuron_class = './mat/nerve_ring_classes.txt'

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
        apos = db.mine.maxAnterior(cur,n)
        ppos = db.mine.maxPosterior(cur,n)
        if apos:
            apos = (maxpos - int(apos))*scale
        else:
            apos = -1
        if ppos:
            ppos = (min(maxpos,int(ppos)))*scale
        else:
            ppos = -1
        data.append([n,nclass[n],apos,ppos])
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
