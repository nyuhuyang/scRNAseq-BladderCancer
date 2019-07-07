library(SingleR)
library(Seurat)
library(reshape2)
library(pheatmap)
library(kableExtra)
library(dplyr)
source("R/utils/Seurat3_functions.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file = "data/BladderCancer_mm10_6_20190706.Rda"))
(load(file = "output/singler_F_BladderCancer_mm10_6_20190706.Rda"))
# if singler didn't find all cell labels
length(singler$singler[[1]]$SingleR.single$labels) == ncol(object@assays$RNA@data)
if(length(singler$singler[[1]]$SingleR.single$labels) < ncol(object@assays$RNA@data)){
        all.cell = colnames(object);length(all.cell)
        know.cell = rownames(singler$singler[[1]]$SingleR.single$labels);length(know.cell)
        object = object[,know.cell]
}


table(rownames(singler$singler[[1]]$SingleR.single$labels) %in% colnames(object))
singler$meta.data$orig.ident = object@meta.data$orig.ident # the original identities, if not supplied in 'annot'
singler$meta.data$xy = object@reductions$tsne@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = Idents(object) # the Seurat clusters (if 'clusters' not provided)
save(singler,file="output/singler_F_BladderCancer_mm10_6_20190706.Rda")

##############################
# add singleR label to Seurat
###############################
singlerDF = data.frame("singler1sub" = singler$singler[[1]]$SingleR.single$labels,
                       "singler1main" = singler$singler[[1]]$SingleR.single.main$labels,
                       "orig.ident" = object@meta.data$orig.ident,
                       "UMAP_1" = object@reductions$umap@cell.embeddings[,"UMAP_1"],
                       "UMAP_2" = object@reductions$umap@cell.embeddings[,"UMAP_2"],
                       row.names = rownames(singler$singler[[1]]$SingleR.single$labels))

table(rownames(singlerDF) %in% colnames(object))
head(singlerDF)
apply(singlerDF,2,function(x) length(unique(x)))

##############################
# check the spearman correlation
###############################
#Or by all cell types (showing the top 50 cell types):
jpeg(paste0(path,"DrawHeatmap_sub1.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single, top.n = 50,normalize = F))
dev.off()
jpeg(paste0(path,"DrawHeatmap_sub1_N.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single,top.n = 50,normalize = T))
dev.off()

#Finally, we can also view the labeling as a table compared to the original identities:

kable(table(singlerDF$singler1sub, singlerDF$orig.ident)) %>%
        kable_styling()
singlerDF$orig.ident %>% table() %>% kable() %>% kable_styling()
singlerDF$singler1sub %>% table() %>% kable() %>% kable_styling()

singlerDF$singler1sub = gsub("Adipocytes|Ependymal|Microglia activated|NPCs|qNSCs",
                             "other cell types",singlerDF$singler1sub)
singlerDF$singler1sub = gsub("Fibroblasts|Fibroblasts activated|Fibroblasts senescent",
                             "Stromal cells",singlerDF$singler1sub)
singlerDF$singler1sub %<>% plyr::mapvalues(from = c("Hepatocytes","Macrophages activated"),
                                           to = c("Epithelial cells","Macrophages"))
##############################
# process color scheme
##############################

singler_colors <- readxl::read_excel("doc/singler.colors.xlsx")
singler_colors1 = as.vector(singler_colors$singler.color1[!is.na(singler_colors$singler.color1)])
singler_colors1[duplicated(singler_colors1)]
length(singler_colors1)
apply(singlerDF[,c("singler1sub","singler1main")],2,function(x) length(unique(x)))
singlerDF[,c("singler1sub")] %>% table() %>% kable() %>% kable_styling()
object <- AddMetaData(object = object,metadata = singlerDF)
object <- AddMetaColor(object = object, label= "singler1sub", colors = singler_colors1)
Idents(object) <- "singler1sub"
object %<>% sortIdent()
TSNEPlot.1(object, cols = ExtractMetaColor(object),label = F,pt.size = 1,
         label.size = 5, repel = T,do.print = T,title = "All cell types in tSNE plot")

jpeg(paste0(path,"tsne_umap_cell-type.jpeg"), units="in", width=10, height=7,res=600)
UMAPPlot(object, group.by="singler1sub",pt.size = 1,label = F,
         cols = ExtractMetaColor(object),
         label.size = 4, repel = T)+ggtitle("All cell types in UMAP plot")+
        theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"))
dev.off()

save(object,file="data/BladderCancer_mm10_6_20190706.Rda")
##############################
# draw tsne plot
##############################
object@meta.data$orig.ident %<>% factor(levels = c("4950N","8524N","8525N",
                                                   "4950P","8524P","8525P"))
Idents(object) <- "singler1sub"
TSNEPlot.1(object, cols = ExtractMetaColor(object),label = F,pt.size = 1,
         split.by = "orig.ident", group.by = "singler1sub",label.size = 4, repel = T, 
         no.legend = T, do.print = T,
         ncol=3,title = "tSNE plots for all samples")

UMAPPlot.1(object, cols = ExtractMetaColor(object),label = F,pt.size = 1,
           split.by = "orig.ident", label.size = 4, repel = T, 
           no.legend = T, do.print = T,unique.name =T,
           ncol=3,title = "UMAP plots for all samples")
