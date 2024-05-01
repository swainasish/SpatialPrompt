```{r}
library(CARD)
library(MuSiC)
library(Seurat)
library(spacexr)
```
load sc_ref_data
```{r}
sc_data_1 = read.csv("sc_ref_hippo.csv",
                             header = T,row.names = 1)
sc_data_labels = read.csv("sc_label_hippo.csv",row.names = 1)
sc_df = t(sc_data_1)
st_data_1_5k = read.csv("hippo_st_5k.csv",
                             header = T,row.names = 1)
st_data_5k_coordinates = read.csv("hippo_coord_5k.csv",row.names =1)
st_df_5k = t(st_data_1_5k)
```
preprocess
```{r}
st_mat = as.matrix(st_df_5k)
sc_mat = as.matrix(sc_df)
sc_labels = data.frame(sc_data_labels$ref.meta.data.celltype)
rownames(sc_labels) = colnames(sc_mat)
st_coordinates = data.frame(st_data_5k_coordinates$imagecol,st_data_5k_coordinates$imagerow)
colnames(st_coordinates) <- c("x","y")
rownames(st_coordinates) = colnames(st_mat)
colnames(sc_labels) = c("celltype")
```
run RCTD
```{r}
t1 = Sys.time()
rctd_celltypes = factor(sc_labels$celltype)
names(rctd_celltypes) = colnames(sc_mat)
sc_reference=Reference(
    counts=sc_mat,
    cell_types=rctd_celltypes,require_int = F
)
st_data=SpatialRNA(
    counts=st_mat,
    coords=st_coordinates,
    require_int=FALSE
)
myRCTD <- create.RCTD(st_data, sc_reference, max_cores = 1, CELL_MIN_INSTANCE = 1)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
t2 = Sys.time()
print(t2-t1)
weights=myRCTD@results$weights
norm_weights=data.frame(normalize_weights(weights))