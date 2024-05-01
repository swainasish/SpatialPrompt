```{r}
library(CARD)
library(MuSiC)
library(Seurat)
# library(data.table)
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
st_data_1_5k = data.frame(as.matrix(fread("hippo_st_5k.csv",
                             header = T),rownames=1))
st_data_5k_coordinates = read.csv("hippo_coord_5k.csv",row.names =1)
st_df_5k = t(st_data_1_5k)
```
preprocess 5k
```{r}
card_st_mat = as.matrix(st_df_5k)
card_sc_mat = as.matrix(sc_df)
# colnames(card_sc_mat) = 1:44976
sc_labels = data.frame(sc_data_labels$ref.meta.data.celltype)
rownames(sc_labels) = colnames(card_sc_mat)
st_coordinates = data.frame(st_data_5k_coordinates$imagecol,st_data_5k_coordinates$imagerow)
colnames(st_coordinates) <- c("x","y")
rownames(st_coordinates) = colnames(card_st_mat)
colnames(sc_labels) = c("celltype")
sc_labels$sampleInfo = "sample1"
```
card in action 
```{r}
t1 = Sys.time()
CARD_obj = createCARDObject(
	sc_count = card_sc_mat,
	sc_meta = sc_labels,
	spatial_count = card_st_mat,
	spatial_location = st_coordinates,
	ct.varname = "celltype",
	ct.select = unique(sc_labels$celltype),
	sample.varname ="sampleInfo",
	minCountGene = 0,
	minCountSpot = 0) 
CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
t2 = Sys.time()
```