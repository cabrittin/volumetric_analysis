"""
paper_figure.py

Generates plots and tables for:

Brittin, Cook, Hall, Emmons and Cohen. 'Volumetric reconstruction of 
Caenorhabditis elegans nerve ring supports combinatorial CAM expression 
model for synaptic specificity'. (2018) Under review. 

created: Christopher Brittin
date: 01 November 2018 

Synopsis:
  python paper_figures figNum [fout]


Parameters:
  fig (str): Figure number or table number from manuscript.
             e.g. f1a is Figure 1a
                  fs2a is Figure S2a
                  ts1 is Table S1
  -o, --output (str): Output file name (optional).


"""


import sys
sys.path.append(r'./scripts/')
import argparse



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('fig',
                        action="store",
                        help=("Figure number or table number."
                              "Examples: "
                              "To generate Figure 1a pass f1a. "
                              "To generate Figure S2a pass fs2a. "
                              "To generate Table S1 pass ts1." )
                        )
   
    parser.add_argument('-o','--output',
                        dest = 'fout',
                        action="store",
                        required= False,
                        default = None,
                        help="Output file, if you wish to save the output."
                        )

    params = parser.parse_args()

    if params.fig =='f3a':
        import dist_adj_subgrp2
        dist_adj_subgrp2.run(params.fout)

 