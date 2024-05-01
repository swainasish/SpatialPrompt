import pandas as pd
import numpy as np
import scanpy as sc
import time
dt_dir = ""
import cell2location
import matplotlib.pyplot as plt
sc.set_figure_params(dpi_save = 600,figsize=(15,15),fontsize=45)
import datatable as dt
#%% import ref sc data 
sc_data_ref = pd.read_csv(dt_dir+"sc_ref_hippo.csv",index_col=0)
sc_data_labels = pd.read_csv(dt_dir+"sc_label_hippo.csv",index_col=0)
adata_ref = sc.AnnData(sc_data_ref)
adata_ref.obs.loc[:,"Subset"] = pd.Categorical(sc_data_labels.iloc[:,0])
adata_ref.var_names_make_unique()
adata_ref.var['SYMBOL'] = adata_ref.var.index 
adata_ref.var.set_index('SYMBOL', drop=True, inplace=True)
#%% import spatial 5k
st_dt_5k = pd.read_csv(dt_dir+"hippo_st_5k.csv",index_col=0)
adata_vis = sc.AnnData(st_dt_5k)
adata_vis.var['SYMBOL'] = adata_vis.var_names
adata_vis.var.set_index('SYMBOL', drop=True, inplace=True)
adata_vis.var_names_make_unique()
adata_vis.obs['sample'] = "hippo5k"
#%% run 
t1 = time.time()
from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

# filter the object
adata_ref = adata_ref[:, selected].copy()


# prepare anndata for the regression model
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref, 
                        labels_key='Subset' )

# create the regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref) 

# view anndata_setup as a sanity check
mod.view_anndata_setup()

mod.train(max_epochs=250, use_gpu=True)   #change epoch  to  250

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' 
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}' 
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]

intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")

mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver, 
    N_cells_per_location=30,
    detection_alpha=20) 

mod.train(max_epochs=2000, ###########change this 30k
          # train using full data (batch_size=None)
          batch_size=2000, 
          # use all data points in training because 
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=True,
         )

# plot ELBO loss history during training, removing first 100 epochs from the plot
mod.plot_history(1)
plt.legend(labels=['full data training']);

adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': 5000, 'use_gpu': True}
)

adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
t2  = time.time()
print(t2-t1)
#%%
cell2loc_result = adata_vis.obs.loc[:,np.unique(sc_data_labels)]
cell2loc_row_sum = cell2loc_result.sum(axis=1)
cell2loc_result_1 = cell2loc_result.div(cell2loc_row_sum,axis=0)
cell2loc_result_1.to_csv(dt_dir+"5k.csv")