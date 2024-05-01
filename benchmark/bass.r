```{r}
library(BASS)
library(Seurat)
library(tidyverse)
```
load dlpfc dataset 
```{r}
st_data_libd = read.csv("libd_df.csv",
                              header = T,row.names = 1)
st_data_coordinates = data.frame(read.csv("meta_data.csv",row.names =1))
st_df_libd = t(st_data_libd)
st_data_coordinates = st_data_coordinates[,1:2]
colnames(st_data_coordinates) = c("x","y")
```
pre-process bass
```{r}
load("starmap_mpfc.RData")
cnts <- list(st_df_libd)
xys <- lapply(list(st_data_coordinates), function(info.i){
  as.matrix(info.i[, c("x", "y")])
}) # a list of spatial coordinates matrices
BASS <- createBASSObject(cnts, xys, C = 20, R = 7, beta_method = "SW",init_method = "mclust", nsample = 10000)

```
bass_on_action
```{r}
BASS <- BASS.preprocess(BASS, doLogNormalize = TRUE,
  geneSelect = "sparkx", nSE = 3000, doPCA = TRUE, 
  scaleFeature = FALSE, nPC = 20)
BASS <- BASS.run(BASS)
BASS <- BASS.postprocess(BASS)
zlabels <- BASS@results$z # spatial domain labelsBASS <- BASS.run(BASS)
res = data.frame(zlabels)
write.csv(res,"dlpfc_clus.csv")
```