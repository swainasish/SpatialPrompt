import numpy as np
import pandas as pd
import os
import anndata
import scanpy as sc
dataset_dir = "spatial_prompt/"
dt_dir = "benchmark/"
import time
#%% import data 
libd_data_dir = dataset_dir+"benchmark/spatiallibd/st_data/spatialLIBD/"
libd_st_obj = sc.read_visium(libd_data_dir+"151673/")
libd_st_meta = pd.read_table(libd_data_dir+"151673/"+"metadata.tsv")
libd_st_obj.obs.loc[:,"layer"]=pd.Categorical(libd_st_meta.layer_guess_reordered)
libd_x = libd_st_obj.obs.array_row
libd_y = libd_st_obj.obs.array_col
libd_st_array = libd_st_obj.X.toarray()
#%% 
adata = libd_st_obj.copy()
t1 = time.time()
adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4) 
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
# sc.tl.umap(adata)
sc.tl.leiden(adata)
t2=time.time()
print(t2-t1)