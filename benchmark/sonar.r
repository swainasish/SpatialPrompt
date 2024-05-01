```{r}
library(SONAR)
library(here)
library(Matrix)
library(data.table)
library(Seurat)
library(matlabr)
library(R.matlab)
```
```{r}
libd_sc_df = fread("libd_sc_df2.csv")
libd_st_df = read.csv("libd_st_df.csv",row.names = 1)
libd_st_coord = read.csv("libd_st_coord.csv",row.names = 1)
libd_sc_labels = read.csv("libd_sc_labels.csv")
st_mat = as.matrix(t(libd_st_df))
sc_mat = t(as.matrix(libd_sc_df[,2:30063]))
colnames(sc_mat) = 1:44976
sc_labels = data.frame(libd_sc_labels$X0)
rownames(sc_labels) = colnames(sc_mat)
sc_labels$cellname = rownames(sc_labels)
sc_labels2 <- data.frame(sc_labels$cellname,sc_labels$libd_sc_labels.X0)
colnames(sc_labels2) <- c("cellname","celltype")
st_coordinates = data.frame(libd_st_coord$x,libd_st_coord$y)
colnames(st_coordinates) <- c("x","y")
```
```{r}
st_mat <- as.data.frame(st_mat)
sc_mat <- as.data.frame(sc_mat)
overlap_gene <- intersect(rownames(st_mat), rownames(sc_mat))
ref <- sc_mat[overlap_gene,]
spots <- st_mat[overlap_gene,]
#calculate the nUMI and nUMI_spot
rownames(st_coordinates)=as.character(1:3639)
colnames(spots)=rownames(st_coordinates)
nUMI <- colSums(ref)
names(nUMI) <- colnames(ref)
nUMI_spot <- colSums(spots)
names(nUMI_spot) <- colnames(spots)
typeanno <- sc_labels2$celltype
names(typeanno) <- sc_labels2$cellname
typeanno <- as.factor(typeanno)
```
SONAR 
```{r}
processed_input<-SONAR.preprocess(sc_count=ref,sc_cell_type=typeanno,sc_nUMI=nUMI,
                                  sp_coords=st_coordinates,sp_count=spots,sp_nUMI=nUMI_spot,
                                  cores=8,spot_min_UMI =0)
#deliver the preprocessed data to SONAR
code_path <- file.path("sonar/")
trans_data<-SONAR.deliver(processed_data=processed_input,path=code_path)
temp<-dist(st_coordinates)
temp<-Matrix::Matrix(temp)
temp[temp==0] <- NA
mindist <- min(temp,na.rm = T)
h <- 1.2*mindist
SONAR.deconvolute(fname = paste0(code_path,"/SONAR_main.m"),path=code_path,h,wait = TRUE)

```
```{r}
SONAR.results <- readMat("SONAR_results.mat",sep = ",")
SONAR.results <- SONAR.results$JIE
u <- fread(paste0(code_path,"u.txt"))
u[,1] <- NULL
colnames(SONAR.results) <- colnames(u)
spot_name <- read.table(file=paste0(code_path,"coord.txt"),sep=",")
rownames(SONAR.results) <- rownames(spot_name)

#Complete! SONAR.results is the final results
result_path="sonar/"
write.table(SONAR.results,file = paste(result_path,'SONAR.results.txt',sep=""),sep = ",")
```