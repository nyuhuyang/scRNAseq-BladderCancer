########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
library(DropletUtils)
library(scater)
library(EnsDb.Hsapiens.v86)
source("../R/Seurat_functions.R")
source("../R/scatter_utils.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
########################################################################
#
#  0.1~0.5 scater
# 
# ######################################################################
# 0.1. Setting up the data
# 0.1.1 Reading in a sparse matrix
df_samples <- readxl::read_excel("doc/181008_Single_cell_sample list.xlsx")
sample_n <- df_samples$species == "Human"
samples <- df_samples$samples[sample_n]
projects <- df_samples$projects[sample_n]
conditions <- df_samples$conditions[sample_n]
df_samples %>% kable() %>% kable_styling()

sce_list <- list()
for(i in 1:length(samples)){
        fname <- paste0("./data/",samples[i],
                        "/outs/filtered_gene_bc_matrices/hg19")
        sce_list[[i]] <- read10xCounts.1(fname, col.names=TRUE,
                                         add.colnames = samples[i])
}
names(sce_list) <- samples
# 0.1.2 Annotating the rows
for(i in 1:length(samples)){
        rownames(sce_list[[i]]) <- uniquifyFeatureNames(rowData(sce_list[[i]])$ID,
                                                        rowData(sce_list[[i]])$Symbol)
        print(head(rownames(sce_list[[i]]),3))
        print(length(rownames(sce_list[[i]])))
}

# We also identify the chromosomal location for each gene. 
# The mitochondrial percentage is particularly useful for later quality control.
for(i in 1:length(samples)){
        location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce_list[[i]])$ID, 
                           column="SEQNAME", keytype="GENEID")
        rowData(sce_list[[i]])$CHR <- location
        print(summary(location=="MT"))
}

# 0.2 Calling cells from empty droplets#######################################
## --------------------------------------------------------------------------
# 0.2.1 Testing for deviations from ambient expression
## --------------------------------------------------------------------------
# 0.2.2 Examining cell-calling diagnostics
## --------------------------------------------------------------------------
# 0.4 Quality control on the cells#########################################################3
sce_list <- lapply(sce_list, function(x) calculateQCMetrics(x,compact = FALSE,
                        feature_controls=list(Mito=which(location=="MT"))))
########################################################################
## --------------------------------------------------------------------------
## ----qchist, Histograms of QC metric distributions in the dataset."----

# Ideally, we would remove cells with low library sizes or total number of expressed features as described previously.
# However, this would likely remove cell types with low RNA content,
# especially in a heterogeneous population with many different cell types.
# Thus, we use a more relaxed strategy and only remove cells with large mitochondrial proportions,
# using it as a proxy for cell damage. 
# (Keep in mind that droplet-based datasets usually do not have spike-in RNA.)
# Low-quality cells are defined as those with extreme values for these QC metrics and are removed.
for(i in 1:length(samples)){
        high.mito <- isOutlier(sce_list[[i]]$pct_counts_Mito, nmads=3, type="higher")
        low.lib <- isOutlier(sce_list[[i]]$log10_total_counts, type="lower", nmad=3)
        low.genes <- isOutlier(sce_list[[i]]$log10_total_features_by_counts, type="lower", nmad=3)
        discard <- high.mito | low.lib | low.genes
        data.frame(HighMito= sum(high.mito),LowLib=sum(low.lib), 
                   LowNgenes=sum(low.genes),Discard=sum(discard))
        sce_list[[i]] <- sce_list[[i]][,!discard]
        print(summary(!discard))
}
# 0.5 Examining gene expression ############################
## --------------------------------------------------------------------------
########################################################################
#
#  0.6 log tranform
# 
# ######################################################################
# Use natural Log transform to fit Seurat

for(i in 1:length(sce_list)){
        logcounts(sce_list[[i]]) <- as(log(assay(sce_list[[i]], "counts")+1),"dgCMatrix")
}
save(sce_list, file = "data/sce_list_Human_20181109.Rda")