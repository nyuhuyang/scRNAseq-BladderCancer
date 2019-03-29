########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("DropletUtils","dplyr","scater","kableExtra","Matrix"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("../R/Seurat_functions.R")
source("../R/scatter_utils.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
########################################################################
#
#  0.1~0.5 scater
# 
# ######################################################################
# 0.1. Setting up the data
# read sample summary list
# read sample summary list
args <- commandArgs(trailingOnly = TRUE)
df_samples <- readxl::read_excel(args[1])
print(df_samples)
df_samples = as.data.frame(df_samples)
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% c("control",paste0("test",1:10)))
#df_samples[sample_n,] %>% kable() %>% kable_styling()
table(df_samples$tests);nrow(df_samples)
list_samples <- lapply(colnames(df_samples), function(col) df_samples[,col])
names(list_samples) = colnames(df_samples)
keep = sapply(list_samples, function(n) length(n[!is.na(n)])>1)
list_samples =list_samples[keep]
samples <- list_samples$sample

sce_list <- list()
# select species
if(unique(list_samples$species) == "Homo_sapiens") species <- "hg19"
if(unique(list_samples$species) == "Mus_musculus") species <- "mm10"
if(unique(list_samples$species) == "Danio_rerio") species <- "danRer10"
if(species == "hg19") suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
if(species == "mm10") suppressPackageStartupMessages(library(EnsDb.Mmusculus.v79))

for(i in 1:length(list_samples$sample)){
        fname <- paste0("./data/",list_samples$sample.id[i],
                        "/outs/filtered_gene_bc_matrices/",species)
        sce_list[[i]] <- read10xCounts.1(fname, col.names=TRUE,
                                         add.colnames = list_samples$sample[i])
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
        if(species == "hg19") {
                location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce_list[[i]])$ID, 
                           column="SEQNAME", keytype="GENEID")
        }
        if(species == "mm10") {
                location <- mapIds(EnsDb.Mmusculus.v79, keys=rowData(sce_list[[i]])$ID, 
                                   column="SEQNAME", keytype="GENEID")
        }
        rowData(sce_list[[i]])$CHR <- location
        print(summary(location=="MT"))
}

# 0.4 Quality control on the cells#########################
# It is entirely possible for droplets to contain damaged or dying cells,
# which need to be removed prior to downstream analysis. 
# We compute some QC metrics using  calculateQCMetrics() (McCarthy et al. 2017) 
# and examine their distributions in Figure 2.
sce_list <- lapply(sce_list, function(x) calculateQCMetrics(x,compact = FALSE,
                        feature_controls=list(Mito=which(location=="MT"))))

for(i in 1:length(sce_list)){
        logcounts(sce_list[[i]]) <- as(log1p(assay(sce_list[[i]], "counts")),"dgCMatrix")
}

########################################################################

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

########################################################################
#
#  0.6 natural Log transform 
# 
# ######################################################################
# Use natural Log transform to fit Seurat

for(i in 1:length(sce_list)){
        logcounts(sce_list[[i]]) <- as(log1p(assay(sce_list[[i]], "counts")),"dgCMatrix")
}
save(sce_list, file = paste0("data/","sce_list_",species,length(sample_n),"_",gsub("-","",Sys.Date()),".Rda"))


