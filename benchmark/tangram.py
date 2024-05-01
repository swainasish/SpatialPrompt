#%% import libs
import pandas as pd
import numpy as np
import tangram as tg
import scanpy as sc
import time
dt_dir = "benchmark/"
#%% import ref sc data 
sc_data_ref = pd.read_csv(dt_dir+"hippo_t1_t4/sc_ref/sc_ref_hippo.csv",index_col=0)
sc_data_labels = pd.read_csv(dt_dir+"hippo_t1_t4/sc_ref/sc_label_hippo.csv",index_col=0)
ad_sc = sc.AnnData(sc_data_ref)
ad_sc.obs["celltype_key"] = pd.Categorical(sc_data_labels.iloc[:,0])
#%% tangram 5k 
st_dt_5k = pd.read_csv(dt_dir+"hippo_t1_t4/5k/hippo_st_5k.csv",index_col=0)
ad_sp = sc.AnnData(st_dt_5k)
#%% run 5k 
start = time.time()
sc.pp.normalize_total(ad_sc)
celltype_counts = ad_sc.obs["celltype_key"].value_counts()
celltype_drop = celltype_counts.index[celltype_counts < 2]
print(f'Drop celltype {list(celltype_drop)} contain less 2 sample')
ad_sc = ad_sc[~ad_sc.obs["celltype_key"].isin(celltype_drop),].copy()
sc.tl.rank_genes_groups(ad_sc, groupby="celltype_key", use_raw=False)
markers_df = pd.DataFrame(ad_sc.uns["rank_genes_groups"]["names"]).iloc[0:200, :]
print(markers_df)
genes_sc = np.unique(markers_df.melt().value.values)
print(genes_sc)
genes_st = ad_sp.var_names.values
genes = list(set(genes_sc).intersection(set(genes_st)))

tg.pp_adatas(ad_sc, ad_sp, genes=genes)

ad_map = tg.map_cells_to_space(
                   ad_sc,
                   ad_sp,
                   mode='clusters',
                   cluster_label="celltype_key")

tg.project_cell_annotations(ad_map, ad_sp, annotation="celltype_key")

celltype_density = ad_sp.obsm['tangram_ct_pred']
celltype_density = (celltype_density.T/celltype_density.sum(axis=1)).T

celltype_density.to_csv(dt_dir+"hippo_t1_t4/results_outputs/tangram/5k.csv")

end = time.time()
print(end - start)