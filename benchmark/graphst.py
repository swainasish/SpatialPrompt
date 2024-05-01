#%% import
import scanpy as sc
import numpy as np
import pandas as pd
from GraphST import GraphST
import polars as pl
from sklearn.metrics import mean_squared_error,roc_auc_score,normalized_mutual_info_score,adjusted_rand_score
import time
#%% reading st data 
dataset_dir = ""
libd_data_dir = dataset_dir+"benchmark/spatiallibd/st_data/spatialLIBD/"
adata=sc.read_visium(libd_data_dir+"151673/")
libd_st_meta = pd.read_table(libd_data_dir+"151673/"+"metadata.tsv")
adata.obs.loc[:,"layer"]=pd.Categorical(libd_st_meta.layer_guess_reordered)
adata.var_names_make_unique()
#%%  process st data
GraphST.preprocess(adata)
# build graph
GraphST.construct_interaction(adata)
GraphST.add_contrastive_label(adata)
#%% reading sc ref
libd_sc_df = pl.read_csv("libd_sc_df2.csv")
libd_sc_df = libd_sc_df.to_pandas()
libd_sc_df = libd_sc_df.iloc[:,1:]
adata_sc = sc.AnnData(libd_sc_df)
adata_sc.var_names_make_unique()
#%% sc labels and preprocess
libd_sc_labels = pd.read_csv("libd_sc_labels.csv",)
adata_sc.obs.loc[:,"cell_type"]=pd.Categorical(libd_sc_labels.iloc[:,1])
GraphST.preprocess(adata_sc)
#%%# find overlap genes
from GraphST.preprocess import filter_with_overlap_gene
adata, adata_sc = filter_with_overlap_gene(adata, adata_sc)
#%% in action
# get features
GraphST.get_feature(adata)
import torch
# Run device, by default, the package is implemented on 'cpu'. We recommend using GPU.
device = "cpu"
# Train model
model = GraphST.GraphST(adata, adata_sc, epochs=1200, random_seed=50, device=device, deconvolution=True)
adata, adata_sc = model.train_map()
#%% # Project cells into spatial space
from GraphST.utils import project_cell_to_spot
project_cell_to_spot(adata, adata_sc, retain_percent=0.15)
decoonv_result=adata.obs.iloc[:,3:]