# SpatialPrompt : Spatially aware scalable and accurate tool for spot deconvolution and clustering in spatial transcriptomics

Article url: [abc.com] \
SpatialPromptDB database url: [abc.com] \
Documentation & tutorial url: [abc.com] 
### Graphical overview 
![alt text](tutorials/images/graphical_abstract.png)
SpatialPrompt integrates gene expression, spatial location, and reference single cell RNA sequencing data for spatial deconvolution.\ For spatial clustering, tool requires only gene expression with spatial locations.
### Installation using PIP 
Platforms: Ubuntu: 22.04/20.04, Windows: 10/11, Mac-OS: Ventura - Python: 3.10 Preferable
```
pip install spatialprompt
```
## Qucik start
### Spatial spot deconvolution 
For spot deconvolutiom SpatialPrompt requires:
- `sc_array`: Matrix of Single-cell data, where rows are the cells and columns are the genes.
- `st_array`: Matrix of Spatial data, where rows are the cells and columns are the genes.
- `sc_genes`: Gene names of the `sc_array` matrix.
- `st_genes`: Gene names of the `st_array` matrix.
- `sc_labels`: Cell type annotations of `sc_array`.
- `x_cord`: X coordinate array of spatial data.
- `y_cord`: Y coordinate array of spatial data.
##### Import the library and download the data
```
import scanpy as sc
import spatialprompt as sp
import gdown
gdown.download("https://drive.google.com/uc?id=1n3ACaWjfjXJ8P6IJhxSUax_vBagIUPLV","sc_m_cortex.h5ad", quiet=False)
gdown.download("https://drive.google.com/uc?id=1h7E5nPs2ga1ixOBDjLJKK7V8xJq93ez5","st_m_cortex.h5ad", quiet=False)
```