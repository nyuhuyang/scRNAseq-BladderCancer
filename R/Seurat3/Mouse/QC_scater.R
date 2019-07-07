########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("R.utils","Seurat","dplyr","kableExtra","readxl",
                   "scater","scran","BiocSingular"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
        }))
source("R/utils/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
if(!dir.exists("data")) dir.create("data")
if(!dir.exists("doc")) dir.create("doc")
########################################################################
#
#  1 Data preprocessing
# 
# ######################################################################
#======1.1 Load the data files and Set up Seurat object =========================
# read sample summary list
args <- "doc/190625_scRNAseq_info.xlsx"
#args <- commandArgs(trailingOnly = TRUE)
df_samples <- read_excel(args[1])
print(df_samples)
df_samples = as.data.frame(df_samples)
colnames(df_samples) <- colnames(df_samples) %>% tolower
keep = df_samples$tests %in% paste0("test",c(2,4,5))
df_samples = df_samples[keep,]
df_samples
# check missing data
current <- list.files("data")
(current <- current[!grepl(".Rda|RData",current)])
(missing_data <- df_samples$sample[!(df_samples$sample %in% current)])

# select species
if(unique(df_samples$species) %in% c("Homo_sapiens","Human")) species <- "hg19"
if(unique(df_samples$species) %in% c("Mus_musculus","Mouse")) species <- "mm10"
if(unique(df_samples$species) == "Danio_rerio") species <- "danRer10"
if(species == "hg19") suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
if(species == "mm10") suppressPackageStartupMessages(library(EnsDb.Mmusculus.v79))

message("Copying the datasets")
if(length(missing_data)>0){
        # Move files from Download to ./data and rename them
        for(missing_dat in missing_data){
                old.pth  <- paste("~/Downloads", missing_dat,
                                  "filtered_feature_bc_matrix",sep = "/")
                list.of.files <- list.files(old.pth)
                new.folder <- paste("data", missing_dat,
                                    "filtered_feature_bc_matrix",sep = "/")
                if(!dir.exists(new.folder)) dir.create(new.folder, recursive = T)
                # copy the files to the new folder
                file.copy(paste(old.pth, list.of.files, sep = "/"), new.folder)
                print(list_files <- list.files(new.folder))
        }
}

message("Loading the datasets")
## Load the dataset
Seurat_raw <- list()
Seurat_list <- list()
for(i in 1:length(df_samples$sample.id)){
        Seurat_raw[[i]] <- Read10X(data.dir = paste0("data/",df_samples$sample.id[i],
                                   "/filtered_feature_bc_matrix/"))
        colnames(Seurat_raw[[i]]) = paste0(df_samples$sample[i],"_",colnames(Seurat_raw[[i]]))
        Seurat_list[[i]] <- CreateSeuratObject(Seurat_raw[[i]],
                                               min.cells = 0,
                                               min.features = 0)
}
remove(Seurat_raw);GC()


#======1.1.2 QC before merge =========================
# if args 2 is passed
args[2] = as.character(args[2])
if(is.na(args[2])){
        message("Starting QC")
        cell.number <- sapply(Seurat_list, function(x) length(colnames(x)))
        QC_list <- lapply(Seurat_list, function(x) as.matrix(GetAssayData(x, slot = "counts")))
        median.nUMI <- sapply(QC_list, function(x) median(colSums(x)))
        median.nGene <- sapply(QC_list, function(x) median(apply(x,2,function(y) sum(length(y[y>0])))))
        
        min.nUMI <- sapply(QC_list, function(x) min(colSums(x)))
        min.nGene <- sapply(QC_list, function(x) min(apply(x,2,function(y) sum(length(y[y>0])))))
        
        QC.list <- cbind(df_samples,cell.number, median.nUMI, median.nGene, 
                         min.nUMI,min.nGene, row.names = df_samples$samples)
        write.csv(QC.list,paste0(path,"QC_list.csv"))
        #QC.list %>% kable() %>% kable_styling()
        remove(QC_list,median.nUMI,median.nGene,min.nUMI,min.nGene,QC.list);GC()
}

#========1.1.3 merge ===================================
object <- Reduce(function(x, y) merge(x, y, do.normalize = F), Seurat_list)
remove(Seurat_list);GC()

# read and select mitochondial genes
if(species == "hg19") mito = "^MT-"
if(species == "mm10") mito = "^mt-" # not Mt-
if(species == "danRer10") mito = "^mt-"
message("mito.genes:")
(mito.features <- grep(pattern = mito, x = rownames(object_raw), value = TRUE))

Idents(object) = "orig.ident"
object[["percent.mt"]] <- PercentageFeatureSet(object = object, pattern = mito)
Idents(object) = factor(Idents(object),levels = df_samples$sample)
g1 <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mt"), function(features){
        VlnPlot(object = object, features = features, ncol = 3, pt.size = 0.01)+
                theme(axis.text.x = element_text(size=15),legend.position="none")
})
save(g1,file= paste0(path,"g1","_",length(df_samples$sample),"_",gsub("-","",Sys.Date()),".Rda"))

#========1.2 scatter ===============================
Seurat_list <- lapply(df_samples$sample, function(x) subset(object, idents = x))
remove(object);GC()

for(i in 1:length(df_samples$sample)){
        high.mito <- isOutlier(Seurat_list[[i]]$percent.mt, nmads=3, type="higher")
        low.lib <- isOutlier(log10(Seurat_list[[i]]$nCount_RNA), type="lower", nmad=3)
        low.genes <- isOutlier(log10(Seurat_list[[i]]$nFeature_RNA), type="lower", nmad=3)
        discard <- high.mito | low.lib | low.genes
        print(data.frame(HighMito= sum(high.mito),LowLib=sum(low.lib), 
                         LowNgenes=sum(low.genes),Discard=sum(discard)))
        Seurat_list[[i]] <- Seurat_list[[i]][,!discard]
        #print(summary(!discard))
        print(i)
}

sce_list <- lapply(Seurat_list, as.SingleCellExperiment)
names(sce_list) = df_samples$sample
remove(Seurat_list);GC()

# cluster
set.seed(1000)
clusters_list <- lapply(sce_list,function(x){
        quickCluster(x, use.ranks=FALSE, BSPARAM=IrlbaParam())
}) 
sapply(clusters_list, table)

# computeSumFactors and normalize
sce_list <- mapply(function(x,y){
        computeSumFactors(x, min.mean=0.1, cluster=y)},
        x=sce_list,y=clusters_list)
lapply(sce_list,function(x) summary(sizeFactors(x)))
remove(clusters_list);GC()
#plot(sce_list[[1]]$nCount_RNA, sizeFactors(sce_list[[1]]), log="xy")
sce_list <- lapply(sce_list, normalize)

save(sce_list, file = paste0("data/","sce_",species,"_",length(df_samples$samples),"_",gsub("-","",Sys.Date()),".Rda"))
