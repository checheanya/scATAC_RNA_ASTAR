import numpy as np
import pandas as pd
from sklearn.feature_selection import r_regression
import multiprocessing as mp
import os
from multiprocessing import Pool


# initiallizing the pairwise corr array
first = r_regression(matrix_for_genes, matrix_for_cells[0])
corr_matrix_genes = np.array(first)

def func_calc(gene_a):
    return r_regression(matrix_for_genes, matrix_for_cells[gene_a])   # X = (n_samples, n_features), Y = n samples

if __name__ == '__main__':
    N = mp.cpu_count()

    matrix_for_genes = np.load('/home/annac/datasets/kidney_10x/matrix_for_genes')
    matrix_for_cells = np.load('/home/annac/datasets/kidney_10x/matrix_for_cells')

    with mp.Pool(processes = N) as p:
        result = p.map(func_calc, range(1, 67150))

# saving the fina;e corr martix 
final_matrix = np.array(result)
np.save('corr_matrix_genes', final_matrix)
