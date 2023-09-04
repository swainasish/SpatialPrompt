## Spatial spot-deconvolution 

```
spatialprompt.SpatialDeconvolution.predict_cell_prop(sc_array, st_array,sc_genes, st_genes, sc_labels, x_cord, y_cord,
n_hvgs=1000, min_cell=10, max_cell=15, return_prop=True,
spot_ratio=[0.33, 0.33, 0.33], n_neighbor=45, n_itr=3)
```

### Description
This program performs spot deconvolution in spatial data using scRNA-seq data reference.

### Parameters
- `sc_array`: Matrix of Single-cell data, where rows are the cells and columns are the genes.
- `st_array`: Matrix of Spatial data, where rows are the cells and columns are the genes.
- `sc_genes`: Gene names of the `sc_array` matrix.
- `st_genes`: Gene names of the `st_array` matrix.
- `sc_labels`: Cell type annotations of `sc_array`.
- `x_cord`: X coordinate array of spatial data.
- `y_cord`: Y coordinate array of spatial data.
- `n_hvgs` (default=1000): Number of high variance genes to consider for analysis.
- `min_cell` (default=10): Minimum number of cells to simulate the spatial spot.
- `max_cell` (default=15): Maximum number of cells to simulate the spatial spot.
- `return_prop` (default=True): Return proportions of cell types if true. Else return the cell type having a higher proportion.
- `spot_ratio` (default=[0.33, 0.33, 0.33]): Ratio of proportions of spots to be simulated using criteria 1/2/3 mentioned in the paper. If the labels are ambiguous cell types (e.g., EX_L3_4_5 have cell types of L3 AND L4), then `spot_ratio` should be provided as a list, e.g., `[0, 0, 1]`.
- `n_neighbor` (default=45): Number of neighbors to consider for weighted mean expression calculation.
- `n_itr` (default=3): Number of iterations message passing layer pull information from neighbors.

### Usage
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