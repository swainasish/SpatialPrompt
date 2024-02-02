## Guideline to select single cell RNA sequencing reference dataset

#### When selecting a scRNA-seq reference dataset for spot deconvolution, it is crucial to ensure a comprehensive representation of cell types and expression profiles. Follow these guidelines for an effective dataset selection:

### 1. Diversity of Cell Types:
Choose a scRNA-seq dataset that encompasses a broad range of cell types relevant to the biological context under investigation. Ensure the dataset covers both major and minor cell populations to capture the full spectrum of cellular diversity.

### 2. Expression Profile Variation:
Opt for a dataset that exhibits diverse expression profiles, representing a wide array of gene expression patterns across different cell types. This variation is essential for robust spot deconvolution and accurate characterization of the spatial heterogeneity within the sample.

### 3. Species Compatibility:
Ensure that the scRNA-seq reference dataset is derived from the same species as the spatial transcriptomics data. Consistency in species ensures the biological relevance and accurate representation of cellular characteristics within the spatial microenvironment.


 Considering these guidelines, our database contains single-cell RNA sequencing datasets from both mouse and human, providing a valuable resource for selecting appropriate reference datasets based on the specific requirements of your spatial transcriptomics study.

## Steps to annotate the cell types in scRNA-seq datasets 

### 1.Filtering Criteria:

Cells with fewer than 500 reads and genes with a total read count of less than 1000 were excluded.
Cells with more than 20% mitochondrial genes were also excluded from the analysis.
### 2.Normalization and Transformation:

Counts per million normalization and logarithmic transformation were applied to each dataset to ensure comparability.
### 3.Cluster Identification:

In cases where cell annotations were not provided by the original study, clusters were determined using the pp.neighbors function available in Scanpy.
Identification of Cell Type or Cluster-Specific Markers:

The tl.rank_genes_groups function was employed to identify markers specific to each cell type or cluster.
### 4.Validation of Markers:

The identified markers were manually validated using external databases, namely CellMarker and PanglaoDB.
The top 20 markers with the highest log fold change values were selected for further validation.