
```{r}
library("Seurat")
library("BayesSpace")
```
load spatialLIBD
```{r}
st_data_libd = read.csv("libd_df.csv",
                              header = T,row.names = 1)
st_data_coordinates = data.frame(read.csv("meta_data.csv",row.names =1))
st_df_libd = t(st_data_libd)
```
analysis
```{r}
library(SpatialExperiment)
col_data_libd = data.frame(barcode = colnames(st_df_libd),imagerow = st_data_coordinates$pixel_row, 
                         imagecol=st_data_coordinates$pixel_col,
                        row= st_data_coordinates$array_row,
                        col= st_data_coordinates$array_col)


sce <- SingleCellExperiment(assays=list(counts=as(st_df_libd, "dgCMatrix")),
                            colData=col_data_libd)
counts <- assay(sce, "counts")
libsizes <- colSums(counts)
size.factors <- libsizes/mean(libsizes)
logcounts(sce) <- log2(t(t(counts)/size.factors) + 1)
t1 = Sys.time()
sce <- spatialPreprocess(sce, platform="Visium", 
                              n.PCs=7, n.HVGs=2000, log.normalize=F)
q <- 7  # Number of clusters
d <- 15  # Number of PCs

sce<- spatialCluster(sce, q=q, d=d, platform='Visium',
                        nrep=10000, gamma=3, save.chain=TRUE)

## We recorded the cluster labels to match the expected brain layers
t2 = Sys.time()
labelslibd <- data.frame(sce$spatial.cluster)
write.csv(labelslibd,"libd_manual.csv")
```