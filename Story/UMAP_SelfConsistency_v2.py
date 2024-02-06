# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 19:25:13 2024

@author: Jack
"""

import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import pandas as pd
import os
import anndata as ad
from scipy import interpolate

os.chdir(r'/home/formanj/scRNAseq_model/UpDownProject/Story')
#%% Functions 
def dayseperation(adata):
    # dealing with HTOs
    HTOs = [HTOs for HTOs in adata.var_names if 'DAY' in HTOs]
    print(HTOs)
    htodata = adata[:, HTOs]
    sc.pp.normalize_total(htodata, target_sum=1e4) #normalize every cell to 10,000 UMI
    for hto in HTOs:
        plt.plot(htodata[0:20, hto].X.todense())
    plt.show()
    
    adata.obs['time'] = np.zeros(adata.shape[0])
    for i in range(adata.shape[0]):
        groupindex = np.argmax(htodata[i,:].X.toarray())
        adata.obs['time'][i] = htodata.var_names[groupindex]
    
    n_hto = adata.obs['time'].value_counts().to_dict()
    
    print(n_hto)
    
    plt.bar(range(len(n_hto)), list(n_hto.values()), tick_label=list(n_hto.keys()))
    plt.show()
    
    # removes HTOs so that they dont effect future processes, especially after normalization
    non_HTOs = [name for name in adata.var_names if name not in HTOs]
    adata = adata[:, non_HTOs].copy()
    return adata



#%% Load Data
dataloc1 = r"GSM6118768_Galdos_Seq_Run1_filtered_feature_bc_matrix.h5"
dataloc2 = r"GSM6118769_Galdos_Seq_Run2_filtered_feature_bc_matrix.h5"
dataloc3 = r"GSM6118770_Galdos_Seq_Run3_filtered_feature_bc_matrix.h5"

aset1 = sc.read_10x_h5(dataloc1, gex_only=False)
aset2 = sc.read_10x_h5(dataloc2, gex_only=False)
aset3 = sc.read_10x_h5(dataloc3, gex_only=False)


aset1.uns['dataset_name'] = 'set 1'
aset2.uns['dataset_name'] = 'set 2'
aset3.uns['dataset_name'] = 'Set 3'

aset1.obs['batch'] = ['batch 1']*aset1.shape[0]
aset2.obs['batch'] = ['batch 2']*aset2.shape[0]
aset3.obs['batch'] = ['batch 3']*aset3.shape[0]

adata = [aset1, aset2, aset3]

for i_d in range(len(adata)):
    adata[i_d].var_names_make_unique()

#%% Down Sample
for i_d in range(len(adata)):
    sc.pp.filter_cells(adata[i_d], min_genes=200)
    sc.pp.filter_genes(adata[i_d], min_cells=3)

#%% Label days
for i_d in range(len(adata)):
    adata[i_d] = dayseperation(adata[i_d])
    

#%% Quality Control
for i_d in range(len(adata)):
    adata[i_d].var['mt'] = adata[i_d].var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata[i_d], qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    #sc.pl.violin(adata[i_d], ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)
    
    plt.scatter(adata[i_d].obs['n_genes_by_counts'], adata[i_d].obs['total_counts'], c=adata[i_d].obs['pct_counts_mt'])
    plt.xlabel('genes by counts')
    plt.ylabel('total number of reads')
    plt.title(adata[i_d].uns['dataset_name'] + ' before QC')
    plt.colorbar()
    plt.show()
    
    upper_lim = np.quantile(adata[i_d].obs.n_genes_by_counts.values, .98)
    lower_lim = np.quantile(adata[i_d].obs.n_genes_by_counts.values, .02)
    print(f'{lower_lim} to {upper_lim}')
    adata[i_d] = adata[i_d][(adata[i_d].obs.n_genes_by_counts < upper_lim) & (adata[i_d].obs.n_genes_by_counts > lower_lim)]
    adata[i_d] = adata[i_d][adata[i_d].obs.pct_counts_mt < 20]
    
    plt.scatter(adata[i_d].obs['n_genes_by_counts'], adata[i_d].obs['total_counts'], c=adata[i_d].obs['pct_counts_mt'])
    plt.xlabel('genes by counts')
    plt.ylabel('total number of reads')
    plt.title(adata[i_d].uns['dataset_name'] + ' After QC')
    plt.colorbar()
    plt.show()
    
    # remove stimulated group from data
    adata[i_d] = adata[i_d][[i for i in adata[i_d].obs_names if 'WTC' in adata[i_d].obs['time'][i]], :]


#%% Normilization
# Unequal waiting of genes
for i_d in range(len(adata)):
    sc.pp.normalize_total(adata[i_d], target_sum=1e4)
    sc.pp.log1p(adata[i_d])
    sc.pp.highly_variable_genes(adata[i_d], min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=5000)
    adata[i_d].raw = adata[i_d]
    adata[i_d] = adata[i_d][:, adata[i_d].var.highly_variable]
    #sc.pp.regress_out(adata[i_d], ['total_counts', 'pct_counts_mt']) # regresses out effect
    #sc.pp.scale(adata[i_d], max_value=10) # treats all genes the same
    

#%% Integration
# Simple integration
adata = ad.concat(adata) # compress to one anndata object
"""
We need better data integration

"""
# Combat Integration
sc.pp.combat(adata, key='batch')

#%% PCA
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)


#%% Visualize 
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution = 0.25) # Tune this to get biologically relavent groups
sc.pl.umap(adata, color=['leiden'])
sc.pl.umap(adata, color=['time'])

#%% Save Checkpoint
save = adata.copy()
#%% Restart 
adata = save
#%% Lines
inds_starting = np.where(adata.obs['time'] == 'WTC_DAY0')
inds_final = np.where(adata.obs['time'] == 'WTC_DAY30')

adata.obs['Lines'] = ['Main Data']*adata.shape[0]

nlines = 100
npoints = 1000
newdata = np.zeros([nlines*npoints, adata.shape[1]])
listofcelltypes = ['Null']*nlines*npoints
for i in range(nlines):
    i_start = inds_starting[0][int(np.random.rand()*len(inds_starting))]
    i_final = inds_final[0][int(np.random.rand()*len(inds_final))]
    
    startcell = adata[i_start]
    finalcell = adata[i_final]
    
    delta = np.subtract(finalcell.X, startcell.X)/npoints
    
    for j in range(npoints):
        newdata[i*npoints + j] = np.add(startcell.X, delta*j)
        listofcelltypes[i*npoints + j] = 'Line ' + str(i)
    

lineadata = ad.AnnData(newdata)
lineadata.obs['Lines'] = listofcelltypes
lineadata.var_names = adata.var_names

sc.tl.ingest(lineadata, adata, obs='leiden')

adata = ad.concat([adata, lineadata], join='outer')


sc.pl.umap(adata, color=['Lines'])


plt.scatter(adata.obsm['X_umap'][:,0], adata.obsm['X_umap'][:,1], s=0.25)
plt.scatter(lineadata.obsm['X_umap'][:,0], lineadata.obsm['X_umap'][:,1], label=lineadata.obs['Lines'], s=1)
plt.show()

#%% Connecting the Lines
plt.scatter(adata.obsm['X_umap'][:,0], adata.obsm['X_umap'][:,1], s=0.25)
for i in range(nlines):
    data = adata[adata.obs['Lines'] == 'Line ' + str(i)]
    #spl = interpolate.splrep(data.obsm['X_umap'][:,0], data.obsm['X_umap'][:,1])
    plt.plot(data.obsm['X_umap'][:,0], data.obsm['X_umap'][:,1])
plt.show()


#%% Restart 
adata = save

#%% Balls
nballs = 10
npoints = 100
epsilon = np.log10(100+1) # counts
newdata = np.zeros([nballs*npoints, adata.shape[1]])
listofcelltypes = ['Null']*nballs*npoints
for i in range(nballs):
    i_start = int(np.random.rand()*adata.shape[0])
    
    startcell = adata[i_start]
    
    for j in range(npoints):
        delta = ((np.random.rand(adata.shape[1]) - 0.5)*2)*epsilon
        newdata[i*npoints + j] = np.add(startcell.X, delta)
        listofcelltypes[i*npoints + j] = 'Ball ' + str(i)
    

balladata = ad.AnnData(newdata)
balladata.obs['Balls'] = listofcelltypes
balladata.var_names = adata.var_names

sc.tl.ingest(balladata, adata, obs='leiden')

adata = ad.concat([adata, balladata], join='outer')


sc.pl.umap(adata, color=['Balls'])


plt.scatter(adata.obsm['X_umap'][:,0], adata.obsm['X_umap'][:,1], s=0.25)
plt.scatter(balladata.obsm['X_umap'][:,0], balladata.obsm['X_umap'][:,1], label=balladata.obs['Balls'], s=1)
plt.show()













