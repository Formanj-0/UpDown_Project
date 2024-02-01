import scanpy as sc

def load_data(data_path, data_type):
    if data_type == 'h5ad':
        adata = sc.read_h5ad(data_path)
    elif data_type == 'csv':
        adata = sc.read_csv(data_path)
    elif data_type == 'txt':
        adata = sc.read_text(data_path)
    elif data_type == '10x':
        adata = sc.read_10x_mtx(data_path)
    else:
        print('Data type not recognized')
    return adata

def remove_cells(adata, min_genes=None, max_genes=None, min_counts=None, max_counts=None):
    """
    Removes cells from the self.adata object based on the specified gene and count thresholds.
    Updates the self.adata object.
    """
    # Filter cells based on gene counts
    sc.pp.filter_cells(adata, min_genes=min_genes, max_genes=max_genes)
    # Filter cells based on total counts
    sc.pp.filter_cells(adata, min_counts=min_counts, max_counts=max_counts)
    return adata








