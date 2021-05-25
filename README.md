# scRNAseq-BladderCancer
Collaboration with Bishoy Morris Faltas (bmf9003@med.cornell.edu) at Weill Cornell and David Mulholland (david.mulholland@mssm.edu) at Mount Sinai for Invasive Bladder Cancer project.

## R Shinyapp

<a href=https://weillcornellmed.shinyapps.io/Human_BladderCancer_ShinyCell/ > Human bladder cance scRNA-seq</a><br />
![](https://github.com/nyuhuyang/scRNAseq-BladderCancer/blob/master/doc/Rshiny_human_bladder.png)

<a href=https://weillcornellmed.shinyapps.io/Mouse_BladderCancer_ShinyCell/ > Human bladder cance scRNA-seq </a><br />
![](https://github.com/nyuhuyang/scRNAseq-BladderCancer/blob/master/doc/Rshiny_mouse_bladder.png)

## METHOD

Single-cell RNA-seq data is pre-processed with the scater R package. Data normalization, unsupervised cell clustering, and differential expression analysis were carried out by the Seurat R package. Reference-based cell type annotation was carried out using the SingleR R package.

## How to use this Script

#### Key Software Setup
R version 3.6.0 <br />
Seurat_3.0.3 <br />
MAST_1.10.0 <br />
scater_1.12.0 <br />
scran_1.12.1 <br />
SingleR_1.0.1 <br />

After pulling this repository, create folders **_data_** and **_output_** in the top working folder.
Move Cell Ranger analysis results into **_data_** folder.
Tree structure of directory:

### 1. scater.R
<a href="https://github.com/nyuhuyang/scRNAseq-BladderCancer/blob/master/R/Seurat2/Human/scater.R">scater.R</a> for human<br />
<a href="https://github.com/nyuhuyang/scRNAseq-BladderCancer/blob/master/R/Seurat3/Mouse/QC_scater.R">scater.R</a> for mouse<br />
Initial quality control and remove low quality cells.

After running these two scripts, `sce_list_Human_{date}.Rda` and `sce_list_Mouse_{date}.Rda` files will be generated inside

### 2. Seurat_setup.R
<a href="https://github.com/nyuhuyang/scRNAseq-BladderCancer/blob/master/R/Seurat2/Human/Seurat_setup.R">Seurat_setup.R</a> for human<br />
<a href="https://github.com/nyuhuyang/scRNAseq-BladderCancer/blob/master/R/Seurat3/Mouse/Seurat_setup.R">Seurat_setup.R</a> for mouse<br />

Cells with less than 800 genes or 1500 UMIs or more than 15% of mitochondria genes were excluded from the analysis. Gene expression raw counts were normalized following a global-scaling normalization method with a scale factor of 10,000 and natural log transformation, using the Seurat NormalizeData function. The top 2000 highly variable genes were selected using the expression and dispersion (variance/mean) of genes, followed by canonical correlation analysis (CCA) to identify common sources of variation between the patient and normal datasets. The first 20 CCA results were chosen to generate dimensional t-Distributed Stochastic Neighbor Embedding (tSNE) plots, and cell clustering by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm.

Need to modify the code according to the date. After running these two scripts, `BladderCancer_H2_{date}.Rda` and `BladderCancer_H2_{date}.Rda` files will be generated inside **_data_** folder.
 Do not modify any files in **_data_** folder.
 
### 3. SingleR.R
<a href="https://github.com/nyuhuyang/scRNAseq-BladderCancer/blob/master/R/Seurat2/Human/SingleR.R">SingleR.R</a> for Human<br />
<a href="https://github.com/nyuhuyang/scRNAseq-BladderCancer/blob/master/R/Seurat3/Mouse/SingleR.R">SingleR.R</a> for Mouse<br />
Cell types were identified by SingleR (Single-cell Recognition) package. SingleR is a novel computational method for unbiased cell type recognition of scRNA-seq. SingleR leverages reference transcriptomic datasets of pure cell types to infer the cell of origin of each of the single cells independently.

After running this script, `singler_BladderCancer_H2_{date}.RData` and `singler_BladderCancer_M2_{date}.RData` file will be generated inside **_output_** folder.

### 4. Identify_Cell_Types_Manually.R
<a href="https://github.com/nyuhuyang/scRNAseq-BladderCancer/blob/master/R/Seurat2/Human/Identify_Cell_Types_Manually.R">Identify_Cell_Types_Manually.R</a> for Human<br />
<a href="https://github.com/nyuhuyang/scRNAseq-BladderCancer/blob/master/R/Seurat3/Mouse/Identify_Cell_Types_Manually.R">Identify_Cell_Types_Manually.R</a> for Mouse<br />
All clusters are tested against marker genes and gene sets.

Multiple plots and table will be generated, save them when necessary. 

### 5. Differential_analysis.R
<a href="https://github.com/nyuhuyang/scRNAseq-BladderCancer/blob/master/R/Seurat2/Human/Differential_analysis.R">Differential_analysis.R</a> for Human<br />
<a href="https://github.com/nyuhuyang/scRNAseq-BladderCancer/blob/master/R/Seurat3/Mouse/Differential_analysis.R">Differential_analysis.R</a> for Mouse<br />
Modified FindAllMarkers() `FindAllMarkers.UMI()` will generate similar dataframe plus two extra columns `UMI.1` and `UMI.2` to record nUMI. `UMI.1` is average nUMI of current cluster, `UMI.2` is average nUMI of all rest of clusters.<br />
`FindAllMarkers(object, test.use = "MAST")` : MAST (Model-based Analysis of Single Cell Transcriptomics), a GLM-framework that treates cellular detection rate as a covariate (Finak et al, Genome Biology, 2015)<br />

Below is an example of a Differential analysis output file.

| gene |   p_val | avg_logFC |  pct.1 |  pct.2 | p_val_adj |  UMI.1 |  UMI.2 |  cluster
| -----    | ------  | -------- | ----  | ----- | ------- | ------- | ------| --- |
| Psca    | 0 | 3.9340 | 0.939 | 0.055 | 0 | 3.5565 | 0.0339 | 0
| Ppbp    | 0 | 2.9163 | 0.99 | 0.161 | 0 | 3.0622 | 0.1834 | 0
| Ltf    | 0 | 2.9105 | 0.959 | 0.042 | 0 | 2.6070 | 0.0365 | 0
| Ecm1    | 0 | 2.7729 | 0.965 | 0.072 | 0 | 2.6652 | 0.0931 | 0
| Gsto1 | 0 | 2.7221 | 0.995 | 0.035 | 0 | 2.7625 | 0.0496 | 0

The results data frame has the following columns :

gene: gene name.<br />
p_val: p_val is calculated using MAST (Model-based Analysis of Single Cell Transcriptomics, Finak et al., Genome Biology, 2015) <br />
avg_logFC: log fold-change of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group.<br />
pct.1: The percentage of cells where the gene is detected in the first group.<br />
pct.2: The percentage of cells where the gene is detected in the second group.<br />
p_val_adj: Adjusted p-value, based on Bonferroni correction<br />
UMI.1 is average nUMI of the current cluster.<br />
UMI.2 is average nUMI of rest of clusters.<br />
cluster: either cell type or corresponding cluster.