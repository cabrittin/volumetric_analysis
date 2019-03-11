"""
test_exp_matrix.py

Test expression matrix object

"""

import sys
sys.path.append('./volumetric_analysis')
import argparse
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from cam.expression import Matrix

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('matrix',
                        action = 'store',
                        help = "Path to file with expression matrix")

    params = parser.parse_args()


    M = Matrix(params.matrix)

    print(M.cells)
    print(M.genes)
    print(M.M.shape)
    M.M[M.M<500] = 0
    M.M[M.M>0] = 1

    cmap = ListedColormap(['w', 'k'])
    plt.matshow(M.M,cmap=cmap)
    plt.show()

