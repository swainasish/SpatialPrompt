#%%
import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
import stlearn as st
import argparse
dataset_dir = "spatial_prompt/"
dt_dir = "benchmark/"
from scipy.sparse import issparse
import random, torch
import warnings
warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import SpaGCN as spg
import time
#In order to read in image data, we need to install some package. Here we recommend package "opencv"
#inatll opencv in python
#!pip3 install opencv-python
import cv2
dataset_dir = "spatial_prompt/"

#%% data load 
from scanpy import read_10x_h5
spagcn_dir = dataset_dir+"benchmark/spatiallibd/st_data/spatialLIBD/spagcn/151673/"
adata = read_10x_h5(spagcn_dir+"expression_matrix.h5")
spatial=pd.read_csv(spagcn_dir+"positions.txt",sep=",",header=None,na_filter=False,index_col=0) 
adata.obs["x1"]=spatial[1]
adata.obs["x2"]=spatial[2]
adata.obs["x3"]=spatial[3]
adata.obs["x4"]=spatial[4]
adata.obs["x5"]=spatial[5]
#Select captured samples
adata=adata[adata.obs["x1"]==1]
adata.var_names=[i.upper() for i in list(adata.var_names)]
adata.var["genename"]=adata.var.index.astype("str")

#Read in hitology image
img=cv2.imread(spagcn_dir+"histology.tif")
#%% run spagcn 
#Set coordinates
adata.obs["x_array"]=adata.obs["x2"]
adata.obs["y_array"]=adata.obs["x3"]
adata.obs["x_pixel"]=adata.obs["x4"]
adata.obs["y_pixel"]=adata.obs["x5"]
x_array=adata.obs["x_array"].tolist()
y_array=adata.obs["y_array"].tolist()
x_pixel=adata.obs["x_pixel"].tolist()
y_pixel=adata.obs["y_pixel"].tolist()
#Run SpaGCN
t1 = time.time()
adata.obs["pred"]= spg.detect_spatial_domains_ez_mode(adata, img, x_array, y_array, x_pixel, y_pixel, n_clusters=7, histology=True, s=1, b=49, p=0.5, r_seed=100, t_seed=100, n_seed=100)
adata.obs["pred"]=adata.obs["pred"].astype('category')
t2 = time.time()
print(t2-t1)
adata.obs["refined_pred"]=spg.spatial_domains_refinement_ez_mode(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), x_array=x_array, y_array=y_array, shape="hexagon")
adata.obs["refined_pred"]=adata.obs["refined_pred"].astype('category')
spagecn_libd = pd.DataFrame(adata.obs["refined_pred"])