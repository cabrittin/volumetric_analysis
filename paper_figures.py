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
    elif params.fig == 'f3b':
        import tpair_adj_deg
        tpair_adj_deg.run(params.fout)
    elif params.fig == 'f3c':
        import tpair_adj_weight
        tpair_adj_weight.run(params.fout)
    elif params.fig == 'f3d':
        import compare_neigh_overlap
        compare_neigh_overlap.run(params.fout)
    elif params.fig == 'f4a':
        import dist_confrac_subgrp
        dist_confrac_subgrp.run(params.fout)
    elif params.fig == 'f4b':
        import tpair_confrac
        tpair_confrac.run(params.fout)
    elif params.fig == 'f5a':
        import discrepant_deg_contralateral
        discrepant_deg_contralateral.run(params.fout)
    elif params.fig in ['f5b','f5c']:
        import discrepant_adj_account
        discrepant_adj_account.run(params.fout)
    elif params.fig == 'f5d':
        import dist_adj_weight_decision
        dist_adj_weight_decision.run(params.fout)
    elif params.fig == 'f5e':
        import logistic_test
        logistic_test.run(params.fout)
    elif params.fig == 'f6a':
        import synaptic_specificity
        synaptic_specificity.run(params.fout)
    elif params.fig == 'f6b':
        import pre_post_specificity
        pre_post_specificity.run(params.fout)
    elif params.fig == 'f6c':
        import synapse_pos_specificity
        synapse_pos_specificity.run(params.fout)
    elif params.fig == 'f6d':
        import tpair_syn_adj_ratio
        tpair_syn_adj_ratio.run(params.fout)
    elif params.fig == 'f7a':
        import cam_expression_matrix
        cam_expression_matrix.run(params.fout)
    elif params.fig == 'f7b':
        print("Figure generated with cytoscape")
    elif params.fig == 'f7c':
        import cam_lus
        cam_lus.run(params.fout)
    elif params.fig == 'fs1':
        print("Images take from our web app at "
              "http://wormwiring.org/apps/neuronVolume.")
    elif params.fig == 'fs2a':
        import dist_adj_deg
        dist_adj_deg.run(params.fout)
    elif params.fig == 'fs2b':
        import dist_adj_weight
        dist_adj_weight.run(params.fout)
    elif params.fig == 'fs2c':
        import dist_syn_deg
        dist_syn_deg.run(params.fout)
    elif params.fig == 'fs4':
        import dist_confrac_muscle_correct
        dist_confrac_muscle_correct.run(params.fout)
    
