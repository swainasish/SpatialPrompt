# SpatialPrompt: Scalable and efficient algorithm for cell type deconvolution and clustering in spatial transcriptomic
## Installation 
### Using conda
#### Create a conda environment having Python version 3.10
    conda create -n spatialpromptENV python=3.10
#### activate the environment 
    conda activate spatialpromptENV
#### Install using pip
    pip install spatialprompt

### Without conda, using pip can be installed directly 
###### Platforms : Ubuntu: 22.04/20.04, Windows: 10/11, Mac-OS: Ventura - Python: 3.10 Preferable
    pip install spatialprompt

## Spatial spot-deconvolution 
This program performs spot deconvolution in spatial data using scRNA-seq data reference.
```python
import spatialprompt as sp

# Create an instance of Spatialdeconvolution 
deconv_model = sp.SpatialDeconvolution()

# Example call to predict_cell_prop
result = deconv_model.predict_cell_prop(sc_array, st_array, 
	                          sc_genes,  st_genes, 
	                          sc_labels, 
	                          x_cord, y_cord)
```
## Spatial clustering
This program performs spatial clustering for spatial data.
```python
import spatialprompt as sp

# Create an instance of Spatialclustering
clus_model = sp.SpatialCluster()

# Example call to predict_cell_prop
cortex_clus_annotations = clus_model.fit_predict(cortex_st_mat,
                                            x_cor = cortex_x,
                                            y_cor = cortex_y,
                                            n_cluster=20)
```
