```{r}
#.libPaths("/home/asish/R/x86_64-pc-linux-gnu-library/4.2")
library(MuSiC)
library(Seurat)
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater)
library(scran)
library(scuttle)
library(SPOTlight)
```
load sc_ref_data
```{r}
sc_data_1 = read.csv("sc_ref_hippo.csv",
                             header = T,row.names = 1)
sc_data_labels = read.csv("sc_label_hippo.csv",row.names = 1)
sc_df = t(sc_data_1)
```
load 5k
```{r}
st_data_1_5k = read.csv("/media/asish/New Volume/Datasets_swain_asish/spatial_prompt/benchmark/hippo_t1_t4/5k/hippo_st_5k.csv",
                              header = T,row.names = 1)
st_data_5k_coordinates = read.csv("hippo_coord_5k.csv",row.names =1)
st_df_5k = t(st_data_1_5k)
```
create sc and st object 
```{r}
sc_obj = SingleCellExperiment(assays = (counts=sc_df))
names(assays(sc_obj)) = "counts"
st_obj = SingleCellExperiment(assays = (counts=st_df_5k))
names(assays(st_obj)) = "counts"
sce <- logNormCounts(sc_obj)
sce$free_annotation = sc_data_labels$ref.meta.data.celltype
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, n = 3000)
colLabels(sce) <- colData(sce)$free_annotation
genes <- !grepl(pattern = "^Rp[l|s]|Mt", x = rownames(sce))
# Compute marker genes
mgs <- scoreMarkers(sce, subset.row = genes,groups = sce$free_annotation)
```
keep relevant gene
```{r}
mgs_fil <- lapply(names(mgs), function(i) {
    x <- mgs[[i]]
    # Filter and keep relevant marker genes, those with AUC > 0.8
    x <- x[x$mean.AUC > 0.8, ]
    # Sort the genes from highest to lowest weight
    x <- x[order(x$mean.AUC, decreasing = TRUE), ]
    # Add gene and cluster id to the dataframe
    x$gene <- rownames(x)
    x$cluster <- i
    data.frame(x)
})
mgs_df <- do.call(rbind, mgs_fil)
```
down_sample
```{r}
# split cell indices by identity
idx <- split(seq(ncol(sce)), sce$free_annotation)
n_cells <- 100
cs_keep <- lapply(idx, function(i) {
    n <- length(i)
    if (n < n_cells)
        n_cells <- n
    sample(i, n_cells)
})
sce <- sce[, unlist(cs_keep)]
```

deconv
```{r}
t1 = Sys.time()
res <- SPOTlight(
    x = sce,
    y = st_obj,
    groups = sce$free_annotation,
    mgs = mgs_df,
    hvg = hvg,
    weight_id = "mean.AUC",
    group_id = "cluster",
    gene_id = "gene")
t2 = Sys.time()
print(t2-t1)
mat5k <- res$mat
write.csv(mat5k,"5k.csv")
```