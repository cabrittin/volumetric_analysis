"""
load_scrnaseq.py

Loads data from GSE126954

See https://www.biorxiv.org/content/10.1101/565549v2

created: Christopher Brittin
date: 09 March 2019
"""
import argparse

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('fin',
                        action="store",
                        help="Path the table of terminal expression")

    params = parser.parse_args()
    
    cells = []

    with open(params.fin,'r') as data:
        next(data)
        for _d in data:
            d = _d.split('\t')
            if ':' in d[2]:
                [_cell,_time] = d[2].split(':')
            else:
                _cell = d[2]
            if _cell not in cells: cells.append(_cell)


    for c in sorted(cells):
        print(c)

    print(len(cells))
