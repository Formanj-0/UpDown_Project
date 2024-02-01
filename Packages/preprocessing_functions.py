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








