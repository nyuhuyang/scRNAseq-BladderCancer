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
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
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
(load(file = "output/20190706/Quality control/g1_6_20190705.Rda"))

object %<>% subset(subset = nFeature_RNA > 400 & nCount_RNA > 1000 & percent.mt < 15)
# FilterCellsgenerate Vlnplot before and after filteration
g2 <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mt"), function(features){
    VlnPlot(object = object, features = features, ncol = 3, pt.size = 0.01)+
        theme(axis.text.x = element_text(size=15),legend.position="none")
})

save(g2,file= paste0(path,"g2_6_20190726.Rda"))
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
(load("data/BladderCancer_mm10_6_20190726.Rda"))
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
#p1 <- UMAPPlot(object, group.by="orig.ident",pt.size = 1,label = F,
#               label.size = 4, repel = T)+ggtitle("Original umap plot")+
#    theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"))

#======1.4 Performing SCTransform and integration =========================
set.seed(100)
object_list <- SplitObject(object, split.by = "groups")
object_list %<>% lapply(SCTransform)
object.features <- SelectIntegrationFeatures(object_list, nfeatures = 3000)
options(future.globals.maxSize= 8388608000)
object_list <- PrepSCTIntegration(object.list = object_list, anchor.features = object.features, 
                                    verbose = FALSE)
anchors <- FindIntegrationAnchors(object_list, normalization.method = "SCT", 
                                           anchor.features = object.features)
object <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
DefaultAssay(object) = "integrated"

remove(object.anchors,object_list);GC()
object  %<>% ScaleData()
object %<>% RunPCA(npcs = 105, verbose = FALSE)
object <- JackStraw(object, num.replicate = 20,dims = 120)
object <- ScoreJackStraw(object, dims = 1:120)
jpeg(paste0(path,"JackStrawPlot.jpeg"), units="in", width=10, height=7,res=600)
JackStrawPlot(object, dims = 100:105)+NoLegend()
dev.off()
npcs=105
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
resolutions = c(seq(0.1,1.5, by = 0.1),seq(2,5, by = 1))
for(i in 1:length(resolutions)){
        object %<>% FindClusters(resolution = resolutions[i], algorithm = 1,verbose = F)
        Progress(i, length(resolutions))
}
object$epi_subtypes <- plyr::mapvalues(object$integrated_snn_res.0.3,
                                       from = c("1","2","3","5","12"),
                                       from = c("1","2","3","5","Epcam+ Fibroblasts"),)
object %<>% RunTSNE(reduction = "pca", dims = 1:npcs)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)

object[["cca.umap"]] <- CreateDimReducObject(embeddings = object@reductions[["umap"]]@cell.embeddings,
                                                 key = "ccaUMAP_", assay = DefaultAssay(object))
p2 <- TSNEPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,
                 label.size = 4, repel = T,title = "CCA tsne plot")

#======1.4 UMAP from SCT =========================
DefaultAssay(object) <- "SCT"
object <- FindVariableFeatures(object, nfeatures = 3000)
object <- ScaleData(object = object,features = VariableFeatures(object))
object <- RunPCA(object, features = VariableFeatures(object),verbose =F,npcs = npcs)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)

object %<>% FindClusters(reduction = "pca",resolution = 0.6,
                         dims.use = 1:npcs,print.output = FALSE)
object %<>% FindClusters(reduction = "pca",resolution = 1.2,
                         dims.use = 1:npcs,print.output = FALSE)
#======1.8 UMAP from harmony =========================
DefaultAssay(object) = "SCT"

npcs = 105
jpeg(paste0(path,"S1_RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
system.time(object %<>% RunHarmony(group.by = "orig.ident", dims.use = 1:npcs,
                                     theta = 1, plot_convergence = TRUE,do_pca = FALSE,
                                     nclust = 50, max.iter.cluster = 100))
dev.off()


object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs)
#system.time(object %<>% RunTSNE(reduction = "harmony", dims = 1:npcs))

object[["harmony.umap"]] <- CreateDimReducObject(embeddings = object@reductions[["umap"]]@cell.embeddings,
                                                 key = "harmonyUMAP_", assay = DefaultAssay(object))

#=======1.9 summary =======================================
jpeg(paste0(path,"S1_TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p0+ggtitle("Clustering without integration")+
              theme(plot.title = element_text(hjust = 0.5,size = 18)),
          p2+ggtitle("Clustering with integration")+
              theme(plot.title = element_text(hjust = 0.5,size = 18)))
dev.off()

object@meta.data$conditions %<>% as.factor
object@meta.data$conditions %<>% factor(levels = c("CD45-positive", "CD45-negative"))

TSNEPlot.1(object, group.by="integrated_snn_res.0.6",pt.size = 1,label = T,no.legend = T,
           label.repel = T, alpha = 1,border = T,
           label.size = 4, repel = T,title = "All cluster in tSNE plot resolution = 0.6",do.print = T)

TSNEPlot.1(object, group.by="integrated_snn_res.0.6",pt.size = 1,label = T,no.legend = T,
           label.repel = T, alpha = 1,border = T,split.by = "conditions",
           label.size = 4, repel = T,title = NULL,do.print = T)

p3 <- TSNEPlot.1(object, group.by="integrated_snn_res.0.6",pt.size = 1,label = T,no.legend = T,
                 label.repel = T, alpha = 1,border = T,
                 label.size = 4, repel = T,title = "Total Clusters",do.print = T)
p4 <- TSNEPlot.1(object, group.by="integrated_snn_res.0.6",pt.size = 1,label = T,no.legend = T,
                 label.repel = T, alpha = 1,border = T, split.by = "conditions",
                 label.size = 4, repel = T,title = NULL, do.print = T)

jpeg(paste0(path,"S1_split_TSNEPlot_all.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p3, p4, align = "h")
dev.off()

Idents(object) = "groups"
for (sample in c("4950", "8524", "8525")) {
    subset_object <- subset(object, idents = sample)
    p5 <- TSNEPlot.1(subset_object, group.by="integrated_snn_res.0.6",pt.size = 1,label = T,no.legend = T,
                     label.repel = T, alpha = 1,label.size = 4, repel = T,border = T, 
                     title = "Total Clusters",do.print = F)
    
    p6 <- TSNEPlot.1(subset_object, group.by="integrated_snn_res.0.6",pt.size = 1,label = T,no.legend = T,
                     label.repel = T, alpha = 1,border = T, split.by = "conditions",
                     label.size = 4, repel = T,title = NULL, do.print = F)
    
    jpeg(paste0(path,"S1_split_TSNEPlot_",sample,".jpeg"), units="in", width=10, height=7,res=600)
    print(plot_grid(p5, p6, align = "h")+
              ggtitle(sample)+
              theme(plot.title = element_text(hjust = 0.5,size = 18)))
    dev.off()
}

object@assays$integrated@scale.data = matrix(0,0,0)
object@assays$RNA@scale.data = matrix(0,0,0)
object@assays$SCT@scale.data = matrix(0,0,0)

format(object.size(object),unit = "GB")

save(object, file = "data/BladderCancer_mm10_6_20190726.Rda")
saveRDS(object, file = "data/BladderCancer_mm10_6_20190726.RDS")
object_data = object@assays$SCT@data
save(object_data, file = "data/object_data_mm10_6_20190726.Rda")

(load(file = "data/BladderCancer_mm10_6_20190726.Rda"))
object %<>% UpdateSeuratObject()
Idents(object) = "cell.types"
UMAPPlot.1(object, do.print = T)
DefaultAssay(object) <- "SCT"
object$mouse <- gsub("N$|P$","",object$orig.ident)
s.genes <- cc.genes$s.genes %>% gsub("MLF1IP","CENPU",.) %>% tolower  %>% Hmisc::capitalize() 
g2m.genes <- cc.genes$g2m.genes %>% tolower  %>% Hmisc::capitalize()
object <- CellCycleScoring(object, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
colnames(object@meta.data) %<>% sub("Phase","cell.cycle.phase",.)
save(object, file = "data/BladderCancer_mm10_6_20190726.Rda")

