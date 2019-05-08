"""
combine_expression_files.py

Combine expression files

"""

import sys
sys.path.append('./volumetric_analysis')
import argparse
from collections import defaultdict

import aux

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('exp_file',
                        action = 'store',
                        help = 'Path to file with list of how to combine files')

    parser.add_argument('duration',
                        action = 'store',
                        help = 'total time covered')

    parser.add_argument('fout',
                        action = 'store',
                        help = 'Path ot output file')

    
    params = parser.parse_args()

    exp = aux.read.into_list2(params.exp_file)
    
    dur = float(params.duration) 
    data = defaultdict(lambda: defaultdict(lambda:0))
    for (f,weight) in exp:
        c = float(weight) / dur
        for [gene,cell,w] in aux.read.into_list2(f): 
            data[gene][cell] += c*float(w) 
    
    with open(params.fout,'w') as fout:
        for g in sorted(data):
            for c in sorted(data[g]):
                tmp = ','.join([g,c,str(data[g][c])+'\n'])
                fout.write(tmp)

