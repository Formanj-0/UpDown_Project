#import inc_dec    # The code to test

import unittest   # The test framework
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


class Test_genericPP(unittest.TestCase):
    self.adata = sc.datasets.blobs(random_state = 0)
    self.adata.write('test.h5ad')

    def test_load_data(self):
        self.assertEqual(self.adata.read_h5ad('test.h5ad'), )

    def test_load_premade_data(self):
        self.assertEqual(, np.all(adata.X == scv.datasets.pancreas().X))
    
    def test_(self):
        self.assertEqual(, np.all(adata.X == scv.datasets.pancreas().X))

if __name__ == '__main__':
    unittest.main()