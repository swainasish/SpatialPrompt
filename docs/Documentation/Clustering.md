## Spatial clustering

```
spatialprompt.SpatialCluster.fit_predict(self,st_array,x_cord,y_cord,
	n_neighbor=20,n_itr=3,n_cluster="auto",n_hvgs=1000)
```

### Description
This program perform spatial clustering for spatial data .
        
### Parameters
- `st_array`: Matrix of Spatial data, where rows are the cells and columns are the genes.
- `x_cord`: X coordinate array of spatial data.
- `y_cord`: Y coordinate array of spatial data.
- `n_hvgs` (default=1000): Number of high variance genes to consider for analysis.
- `n_neighbor` (default=45): Number of neighbors to consider for weighted mean expression calculation.
- `n_itr` (default=3): Number of iterations message passing layer pull information from neighbors.
- `n_cluster` (default:"auto"): Number of clusters needed to perform. 
### Usage
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