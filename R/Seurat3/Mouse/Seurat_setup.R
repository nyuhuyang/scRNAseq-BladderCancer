########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(cowplot)
library(kableExtra)
library(magrittr)
library(harmony)
library(scran)
source("R/utils/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
########################################################################
#
#  1 Seurat setup
# 
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# Load the mouse.eyes dataset

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here
df_samples <- readxl::read_excel("doc/190625_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
keep = df_samples$tests %in% paste0("test",c(2,4,5))
df_samples = df_samples[keep,]
(samples = df_samples$sample)

#======1.2 load  SingleCellExperiment =========================
(load(file = "data/sce_mm10_6_20190705.Rda"))
names(sce_list)
object_list <- lapply(sce_list, as.Seurat)

for(i in 1:length(samples)){
    object_list[[i]]$orig.ident <- df_samples$sample[i]
    object_list[[i]]$conditions <- df_samples$conditions[i]
    object_list[[i]]$groups <- df_samples$groups[i]
    object_list[[i]]$notes <- df_samples$notes[i]
    object_list[[i]]$tests <- df_samples$tests[i]
    }

#========1.3 merge ===================================
object <- Reduce(function(x, y) merge(x, y, do.normalize = F), object_list)
object@assays$RNA@data = object@assays$RNA@data *log(2) # change to natural log
remove(sce_list,object_list);GC()
save(object, file = paste0("data/BladderCancer_mm10_",length(df_samples$sample),"_",gsub("-","",Sys.Date()),".Rda"))

#======1.2 QC, pre-processing and normalizing the data=========================
# store mitochondrial percentage in object meta data
object <- PercentageFeatureSet(object = object, pattern = "^mt-", col.name = "percent.mt")
Idents(object) = "orig.ident"
Idents(object) %<>% factor(levels = samples)
(load(file = paste0(path, "g1_6_20190705.Rda"))

object %<>% subset(subset = nFeature_RNA > 800 & nCount_RNA > 1000 & percent.mt < 15)
# FilterCellsgenerate Vlnplot before and after filteration
g2 <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mt"), function(features){
    VlnPlot(object = object, features = features, ncol = 3, pt.size = 0.01)+
        theme(axis.text.x = element_text(size=15),legend.position="none")
})

save(g2,file= paste0(path,"g2_6_20190705.Rda"))
jpeg(paste0(path,"S1_nGene.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[1]]+ggtitle("nFeature_RNA before filteration")+
                    scale_y_log10(limits = c(100,10000)),
                g2[[1]]+ggtitle("nFeature_RNA after filteration")+
                    scale_y_log10(limits = c(100,10000))))
dev.off()
jpeg(paste0(path,"S1_nUMI.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[2]]+ggtitle("nCount_RNA before filteration")+
                    scale_y_log10(limits = c(500,100000)),
                g2[[2]]+ggtitle("nCount_RNA after filteration")+ 
                    scale_y_log10(limits = c(500,100000))))
dev.off()
jpeg(paste0(path,"S1_mito.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[3]]+ggtitle("mito % before filteration")+
                    ylim(c(0,50)),
                g2[[3]]+ggtitle("mito % after filteration")+ 
                    ylim(c(0,50))))
dev.off()

######################################
(load("data/BladderCancer_mm10_6_20190706.Rda"))
# After removing unwanted cells from the dataset, the next step is to normalize the data.
#BladderCancer <- NormalizeData(object = object, normalization.method = "LogNormalize", 
#                      scale.factor = 10000)
object <- FindVariableFeatures(object = object, selection.method = "vst",
                            num.bin = 20,
                            mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(object), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
jpeg(paste0(path,"VariableFeaturePlot.jpeg"), units="in", width=10, height=7,res=600)
print(plot2)
dev.off()
#======1.3 1st run of pca-tsne  =========================
DefaultAssay(object) <- "RNA"
object <- ScaleData(object = object,features = VariableFeatures(object))
object <- RunPCA(object, features = VariableFeatures(object),verbose =F,npcs = 120)
object <- JackStraw(object, num.replicate = 20,dims = 120)
object <- ScoreJackStraw(object, dims = 1:120)
jpeg(paste0(path,"JackStrawPlot.jpeg"), units="in", width=10, height=7,res=600)
JackStrawPlot(object, dims = 100:110)
dev.off()
npcs =105
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(reduction = "pca",resolution = 0.6,
                       dims.use = 1:npcs,print.output = FALSE)
object %<>% RunTSNE(reduction = "pca", dims = 1:npcs)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)

p0 <- TSNEPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,
                 label.size = 4, repel = T,title = "Original tsne plot")
p1 <- UMAPPlot(object, group.by="orig.ident",pt.size = 1,label = F,
               label.size = 4, repel = T)+ggtitle("Original umap plot")+
    theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"))

#======1.4 Performing CCA integration =========================
set.seed(100)
Idents(object) = "orig.ident"
object_list <- lapply(df_samples$sample,function(x) subset(object,idents=x))
anchors <- FindIntegrationAnchors(object.list = object_list, dims = 1:npcs)
object <- IntegrateData(anchorset = anchors, dims = 1:npcs)
remove(anchors,object_list);GC()
DefaultAssay(object) <- "integrated"
object %<>% ScaleData(verbose = FALSE)
object %<>% RunPCA(npcs = npcs, features = VariableFeatures(object),verbose = FALSE)
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(reduction = "pca",resolution = 0.6,
                         dims.use = 1:npcs,print.output = FALSE)
object %<>% RunTSNE(reduction = "pca", dims = 1:npcs)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)

p2 <- TSNEPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,
                 label.size = 4, repel = T,title = "CCA tsne plot")
p3 <- UMAPPlot(object, group.by="orig.ident",pt.size = 1,label = F,
               label.size = 4, repel = T)+ggtitle("CCA umap plot")+
    theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"))
#======1.6 Apply sctransform normalization =========================
DefaultAssay(object) <- "RNA"
object <- SCTransform(object, verbose = T)
(load(file="data/BladderCancer_mm10_6_20190706.Rda"))
DefaultAssay(object) <- "SCT"
object %<>% RunPCA(npcs = npcs, features = VariableFeatures(object),verbose = FALSE)
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(reduction = "pca",resolution = 0.6,
                         dims.use = 1:npcs,print.output = FALSE)
object %<>% RunTSNE(reduction = "pca", dims = 1:npcs)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)

p4 <- TSNEPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,
                 label.size = 4, repel = T,title = "SCT tsne plot")
p5 <- UMAPPlot(object, group.by="orig.ident",pt.size = 1,label = F,
               label.size = 4, repel = T)+ggtitle("SCT umap plot")+
    theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"))
#======1.7 Performing MNN-based correction =========================
#https://bioconductor.org/packages/3.8/workflows/vignettes/simpleSingleCell/inst/doc/work-5-mnn.html#4_performing_mnn-based_correction
set.seed(100)
original <- lapply(samples, function(x) object@assays$RNA@scale.data[VariableFeatures(object), 
                                                 (object$orig.ident %in% x)])
system.time(mnn.out <- do.call(fastMNN, c(original, list(k=20, d=npcs, auto.order=T,
                                             approximate=TRUE))))
dim(mnn.out$corrected)
rownames(mnn.out$rotation) = VariableFeatures(object)
colnames(mnn.out$rotation) = paste0("MNN_", 1:npcs)
rownames(mnn.out$corrected) = colnames(object)
colnames(mnn.out$corrected) = paste0("MNN_",1:ncol(mnn.out$corrected))
#Storing a new MNN
object[['mnn']] <- CreateDimReducObject(embeddings = mnn.out$corrected,
                                       assay = "RNA",key = "MNN_")
remove(original);GC()
Idents(object) <- "orig.ident"

DimPlot(object = object, reduction = "mnn", pt.size = 0.5)

object <- FindNeighbors(object, reduction = "mnn",dims = 1:npcs)
object <- FindClusters(object, reduction = "mnn",resolution = 0.6,
                       dims.use = 1:npcs,print.output = FALSE) %>%
    RunTSNE(reduction = "mnn", dims = 1:npcs)

p6 <- TSNEPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,
                 label.size = 4, repel = T,title = "MNN tsne plot")
p7 <- UMAPPlot(object, group.by="orig.ident",pt.size = 1,label = F,
               label.size = 4, repel = T)+ggtitle("MNN umap plot")+
    theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"))


#======1.8 RunHarmony=======================
jpeg(paste0(path,"RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
system.time(object %<>% RunHarmony(group.by.vars= "orig.ident", dims.use = 1:npcs,
                                   theta = NULL, plot_convergence = TRUE,epsilon.harmony = -Inf,
                                   nclust = 100, max.iter.cluster = 100))
dev.off()

object %<>% FindNeighbors(reduction = "harmony",dims = 1:npcs)
object %<>% FindClusters(reduction = "harmony",resolution = 0.6,
                                  dims.use = 1:npcs,print.output = FALSE)
object %<>% RunTSNE(reduction = "harmony", dims = 1:npcs)
object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs)

p6 <- TSNEPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,
                 label.size = 4, repel = T,title = "harmony tsne plot")
p7 <- UMAPPlot(object, group.by="orig.ident",pt.size = 1,label = F,
               label.size = 4, repel = T)+ggtitle("harmony umap plot")+
    theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"))

Idents(object) = "integrated_snn_res.0.6"
TSNEPlot.1(object, group.by="integrated_snn_res.0.6",pt.size = 1,label = F,
           label.size = 4, repel = T,title = "All cluster in tSNE plot",do.print = T)

jpeg(paste0(path,"tsne_umap.jpeg"), units="in", width=10, height=7,res=600)
UMAPPlot(object, group.by="integrated_snn_res.0.6",pt.size = 1,label = F,
               label.size = 4, repel = T)+ggtitle("All cluster in UMAP plot")+
    theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"))
dev.off()

object@assays$integrated@scale.data = matrix(0,0,0)
save(object, file = paste0("data/BladderCancer_mm10_",length(df_samples$sample),"_",gsub("-","",Sys.Date()),".Rda"))