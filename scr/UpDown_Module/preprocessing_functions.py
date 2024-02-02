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
    if min_genes is not None:
        sc.pp.filter_cells(adata, min_genes=min_genes)
    if max_genes is not None:
        sc.pp.filter_cells(adata, max_genes=max_genes)
    # Filter cells based on total counts
    if min_counts is not None:
        sc.pp.filter_cells(adata, min_counts=min_counts)
    if max_counts is not None:
        sc.pp.filter_cells(adata, max_counts=max_counts)
    return adata

def remove_genes(adata, min_cells=None, max_cells=None, min_counts=None, max_counts=None):
    """
    Removes genes from the self.adata object based on the specified gene and count thresholds.
    Updates the self.adata object.
    """
    # Filter genes based on gene counts
    if min_cells is not None:
        sc.pp.filter_genes(adata, min_cells=min_cells)
    if max_cells is not None:
        sc.pp.filter_genes(adata, max_cells=max_cells)
    # Filter genes based on total counts
    if min_counts is not None:
        sc.pp.filter_genes(adata, min_counts=min_counts)
    if max_counts is not None:
        sc.pp.filter_genes(adata, max_counts=max_counts)
    return adata



def remove_cells_by_qc_var(adata, qc_var_name, value, operator):
    """
    Removes cells from the self.adata object based on the specified quality control variable and operator.
    Updates the self.adata object.
    """
    if operator == '<':
        adata = adata[adata.obs[qc_var_name] < value, :]
    elif operator == '>':
        adata = adata[adata.obs[qc_var_name] > value, :]
    elif operator == '<=':
        adata = adata[adata.obs[qc_var_name] <= value, :]
    elif operator == '>=':
        adata = adata[adata.obs[qc_var_name] >= value, :]
    elif operator == '==':
        adata = adata[adata.obs[qc_var_name] == value, :]
    else:
        print('Invalid operator')
    return adata




