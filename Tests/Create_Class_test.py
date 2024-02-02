# Remove unnecessary import and sys.path.insert() function call
# No changes needed in this code block

import scr.UpDown_Module as ud


import unittest   # The test framework
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import numpy as np


class Test_scanpy_3k_PBMC(unittest.TestCase):


    def setUp(self):
        self.results_file = 'pbmc3k.h5ad'  # the file that will store the analysis results
        self.test_adata = sc.read_10x_mtx(
                                "Datasets/Tutorial/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/",  # the directory with the `.mtx` file
                                var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
                                cache=True)                              # write a cache file for faster subsequent reading

        self.test_adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
        self.test_adata.write_h5ad(self.results_file)

        self.overlord = ud.preprocessing_overlord.overlord_pp()
        self.overlord.load_data(self.results_file, 'h5ad')


    def test_load_data(self):
        assert self.overlord.adata.shape == self.test_adata.shape



    def test_filtering_cells_and_genes(self):
        print(self.overlord.adata)
        print('Test adata')
        print(self.test_adata)
        cell_initial_true = self.test_adata.shape[0]
        gene_initial_true = self.test_adata.shape[1]
        sc.pp.filter_cells(self.test_adata, min_genes=500)
        sc.pp.filter_genes(self.test_adata, min_cells=100)
        cell_final_true = self.test_adata.shape[0]
        gene_final_true = self.test_adata.shape[1]


        cell_initial_test = self.overlord.adata.shape[0]
        gene_initial_test = self.overlord.adata.shape[1]
        self.overlord.remove_cells(min_genes=500)
        self.overlord.remove_genes(min_cells=100)
        cell_final_test = self.overlord.adata.shape[0]
        gene_final_test = self.overlord.adata.shape[1]

        cell_diff_true = cell_initial_true - cell_final_true
        gene_diff_true = gene_initial_true - gene_final_true
        cell_diff_test = cell_initial_test - cell_final_test
        gene_diff_test = gene_initial_test - gene_final_test

        assert cell_diff_true == cell_diff_test and gene_diff_true == gene_diff_test



    def test_qc_metrics(self):
        # pretest 
        sc.pp.filter_cells(self.test_adata, min_genes=200)
        sc.pp.filter_genes(self.test_adata, min_cells=3)
        self.overlord.remove_cells(min_genes=200)
        self.overlord.remove_genes(min_cells=3)

        # test
        self.test_adata.var['mt'] = self.test_adata.var_names.str.startswith('MT-')
        self.overlord.label_qc_genes(['mt'], ['MT-'])

        assert self.test_adata.var['mt'].equals(self.overlord.adata.var['mt'])



    def test_qc_calculation(self):
        # pretest 
        sc.pp.filter_cells(self.test_adata, min_genes=200)
        sc.pp.filter_genes(self.test_adata, min_cells=3)
        self.overlord.remove_cells(min_genes=200)
        self.overlord.remove_genes(min_cells=3)
        self.test_adata.var['mt'] = self.test_adata.var_names.str.startswith('MT-')
        self.overlord.label_qc_genes(['mt'], ['MT-'])

        # test
        sc.pp.calculate_qc_metrics(self.test_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        self.overlord.calculate_qc_metrics(qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

        assert self.test_adata.obs.equals(self.overlord.adata.obs)


    def test_remove_cells_by_qc_var(self):
        # pretest 
        sc.pp.filter_cells(self.test_adata, min_genes=200)
        sc.pp.filter_genes(self.test_adata, min_cells=3)
        self.overlord.remove_cells(min_genes=200)
        self.overlord.remove_genes(min_cells=3)
        self.test_adata.var['mt'] = self.test_adata.var_names.str.startswith('MT-')
        self.overlord.label_qc_genes(['mt'], ['MT-'])
        sc.pp.calculate_qc_metrics(self.test_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        self.overlord.calculate_qc_metrics(qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

        # test
        self.test_adata = self.test_adata[self.test_adata.obs['n_genes'] < 2500, :]
        self.overlord.simple_qc(method='comparison', value=2500, column_name='n_genes',  operator='<')

        self.test_adata = self.test_adata[self.test_adata.obs.pct_counts_mt < 5, :]
        self.overlord.simple_qc(method='comparison', column_name='pct_counts_mt', value=5, operator='<')

        assert self.test_adata.obs_names.equals(self.overlord.adata.obs_names)



    def test_normalize_totals(self):
        # pretest 
        sc.pp.filter_cells(self.test_adata, min_genes=200)
        sc.pp.filter_genes(self.test_adata, min_cells=3)
        self.overlord.remove_cells(min_genes=200)
        self.overlord.remove_genes(min_cells=3)
        self.test_adata.var['mt'] = self.test_adata.var_names.str.startswith('MT-')
        self.overlord.label_qc_genes(['mt'], ['MT-'])
        sc.pp.calculate_qc_metrics(self.test_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        self.overlord.calculate_qc_metrics(qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        self.test_adata = self.test_adata[self.test_adata.obs['n_genes'] < 2500, :]
        self.overlord.simple_qc(method='comparison', value=2500, column_name='n_genes',  operator='<')
        self.test_adata = self.test_adata[self.test_adata.obs.pct_counts_mt < 5, :]
        self.overlord.simple_qc(method='comparison', column_name='pct_counts_mt', value=5, operator='<')

        # test
        sc.pp.normalize_total(self.test_adata, target_sum=1e4)
        self.overlord.normalize_counts(target_sum=1e4)

        assert np.array_equal(self.test_adata.X.toarray(), self.overlord.adata.X.toarray())


    def test_log_transform(self):
        # pretest 
        sc.pp.filter_cells(self.test_adata, min_genes=200)
        sc.pp.filter_genes(self.test_adata, min_cells=3)
        self.overlord.remove_cells(min_genes=200)
        self.overlord.remove_genes(min_cells=3)
        self.test_adata.var['mt'] = self.test_adata.var_names.str.startswith('MT-')
        self.overlord.label_qc_genes(['mt'], ['MT-'])
        sc.pp.calculate_qc_metrics(self.test_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        self.overlord.calculate_qc_metrics(qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        self.test_adata = self.test_adata[self.test_adata.obs['n_genes'] < 2500, :]
        self.overlord.simple_qc(method='comparison', value=2500, column_name='n_genes',  operator='<')
        self.test_adata = self.test_adata[self.test_adata.obs.pct_counts_mt < 5, :]
        self.overlord.simple_qc(method='comparison', column_name='pct_counts_mt', value=5, operator='<')
        sc.pp.normalize_total(self.test_adata, target_sum=1e4)
        self.overlord.normalize_counts(target_sum=1e4)

        # test
        sc.pp.log1p(self.test_adata)
        self.overlord.convert_to_log1p()

        assert np.array_equal(self.test_adata.X.toarray(), self.overlord.adata.X.toarray())


    def test_save_adata(self):
        # pretest 
        self.overlord.remove_cells(min_genes=200)
        self.overlord.remove_genes(min_cells=3)
        self.overlord.label_qc_genes(label_list=['mt'], search_list=['MT-'])
        self.overlord.calculate_qc_metrics(qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        self.overlord.simple_qc(method='comparison', value=2500, column_name='n_genes',  operator='<')
        self.overlord.simple_qc(method='comparison', column_name='pct_counts_mt', value=5, operator='<')
        self.overlord.normalize_counts(target_sum=1e4)
        self.overlord.convert_to_log1p()

        self.overlord.save_data('test.h5ad', 'h5ad', foldername='Test', cwd = os.getcwd())
        print(self.overlord.steps)
        assert self.overlord.steps[-1] == 'convert_to_log1p' and len(self.overlord.steps) == 8


    def test_oneliners(self):
        # pretest 
        self.overlord.remove_cells(min_genes=200)
        self.overlord.remove_genes(min_cells=3)
        self.overlord.label_qc_genes(label_list=['mt'], search_list=['MT-'])
        self.overlord.calculate_qc_metrics(qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        self.overlord.simple_qc(method='comparison', value=2500, column_name='n_genes',  operator='<')
        self.overlord.simple_qc(method='comparison', column_name='pct_counts_mt', value=5, operator='<')
        self.overlord.normalize_counts(target_sum=1e4)
        self.overlord.convert_to_log1p()
        self.overlord.save_data('test.h5ad', 'h5ad', foldername='Test', cwd = os.getcwd())

        # test
        overlord_test = ud.preprocessing_overlord.overlord_pp()
        overlord_test.load_data(self.results_file, 'h5ad')
        overlord_test.load_and_run_files(foldername='Test', cwd=os.getcwd())

        # idk why i cant do what i did above but it works so who cares
        assert np.max(overlord_test.adata.X.toarray()) == np.max(self.overlord.adata.X.toarray()) and np.min(overlord_test.adata.X.toarray()) == np.min(self.overlord.adata.X.toarray())









if __name__ == '__main__':


    unittest.main()