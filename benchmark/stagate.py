#%%
dataset_dir = "spatial_prompt/"
import pandas as pd
import numpy as np
import scanpy as sc
import STAGATE
import time
import tensorflow as tf
tf.compat.v1.disable_eager_execution()
dt_dir = "benchmark/"
import stlearn as st
#%% import data 
libd_data_dir = dataset_dir+"benchmark/spatiallibd/st_data/spatialLIBD/"
libd_st_obj = sc.read_visium(libd_data_dir+"151673/")
libd_st_meta = pd.read_table(libd_data_dir+"151673/"+"metadata.tsv")
libd_st_obj.obs.loc[:,"layer"]=pd.Categorical(libd_st_meta.layer_guess_reordered)
libd_x = libd_st_obj.obs.array_row
libd_y = libd_st_obj.obs.array_col
libd_st_array = libd_st_obj.X.toarray()
#%% stagate on libd
t1 = time.time()
adata = libd_st_obj.copy()
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
STAGATE.Cal_Spatial_Net(adata, rad_cutoff=40)
STAGATE.Stats_Spatial_Net(adata)
adata = STAGATE.train_STAGATE(adata, alpha=0)
sc.pp.neighbors(adata, use_rep='STAGATE')
sc.tl.umap(adata)
adata = STAGATE.mclust_R(adata, used_obsm='STAGATE', num_cluster=7)
t2=time.time()
print(t2-t1)
stagate_result = pd.DataFrame(adata.obs.mclust)