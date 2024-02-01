import inc_dec    # The code to test

import unittest   # The test framework
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


class Test_scanpy_3k_PBMC(unittest.TestCase):
    results_file = 'write/pbmc3k.h5ad'  # the file that will store the analysis results
    def __init__(self, methodName: str = "runTest") -> None:
        self.test_adata = sc.read_10x_mtx(
                                'data/filtered_gene_bc_matrices/hg19/',  # the directory with the `.mtx` file
                                var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
                                cache=True)                              # write a cache file for faster subsequent reading

        self.adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
        self.test_adata.write_h5ad(self.results_file)
    


    




    def test_load_data(self):
        self.assertEqual(self.adata.read_h5ad('test.h5ad'), )

    def test_load_premade_data(self):
        self.assertEqual(, np.all(adata.X == scv.datasets.pancreas().X))
    
    def test_(self):
        self.assertEqual(, np.all(adata.X == scv.datasets.pancreas().X))

if __name__ == '__main__':
    unittest.main()