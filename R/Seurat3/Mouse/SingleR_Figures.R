library(SingleR)
library(Seurat)
library(reshape2)
library(pheatmap)
library(kableExtra)
library(dplyr)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file = "data/BladderCancer_mm10_6_20190726.Rda"))
(load(file = "output/singler_T_BladderCancer_mm10_6_20190726.Rda"))
# if singler didn't find all cell labels`
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
save(singler,file="output/singler_T_BladderCancer_mm10_6_20190726.Rda")

##############################
# add singleR label to Seurat
###############################
singlerDF = data.frame("singler1sub" = singler$singler[[1]]$SingleR.single$labels,
                       "singler1main" = singler$singler[[1]]$SingleR.single.main$labels,
                       "orig.ident" = object@meta.data$orig.ident,
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
singlerDF$singler1main %>% table() %>% kable() %>% kable_styling()

singlerDF$singler1main = gsub("Hepatocytes","Epithelial cells",singlerDF$singler1main)


##############################
# process color scheme
##############################

singler_colors <- readxl::read_excel("doc/singler.colors.xlsx")
singler_colors1 = as.vector(singler_colors$singler.color1[!is.na(singler_colors$singler.color1)])
singler_colors2 = as.vector(singler_colors$singler.color2[!is.na(singler_colors$singler.color2)])
singler_colors1[duplicated(singler_colors1)];singler_colors2[duplicated(singler_colors2)]
length(singler_colors1);length(singler_colors2)
apply(singlerDF[,c("singler1sub","singler1main")],2,function(x) length(unique(x)))
singlerDF[,c("singler1main")] %>% table() %>% kable() %>% kable_styling()
object <- AddMetaData(object = object,metadata = singlerDF)
object <- AddMetaColor(object = object, label= "singler1main", colors = singler_colors1)
Idents(object) <- "singler1main"
object %<>% sortIdent()
TSNEPlot.1(object, cols = ExtractMetaColor(object),split.by = "conditions",
           label = F,pt.size = 0.3,no.legend = F,label.repel = F,legend.size = 10,
           repel = T,do.return= T,do.print = T,alpha = 0.9,border = T,
         title = "All cell types identified by Mouse RNA-seq reference database")

save(object,file="data/BladderCancer_mm10_6_20190726.Rda")
##############################
# Adjust cell type manually
##############################
Idents(object) = "integrated_snn_res.0.6"
object %<>% RenameIdents('0' = 'Fibroblasts',
                         '1' = 'A',
                         '2' = 'Fibroblasts',
                         '3' = 'Epithelial cells',
                         '4' = 'A',
                         '5' = 'Fibroblasts',
                         '6' = 'A',
                         '7' = 'Fibroblasts',
                         '8' = 'Epithelial cells',
                         '9' = 'A',
                         '10' = 'Endothelial cells',
                         '11' = 'Fibroblasts',
                         '12' = 'T cells',
                         '13' = 'A',
                         '14' = 'B cells',
                         '15' = 'A',
                         '16' = 'A',
                         '17' = 'T cells',
                         '18' = 'Mast cells',
                         '19' = 'Stromal cells',
                         '20' = 'Lymphatic endothelial cells',
                         '21' = 'A')
object[["cell.types"]] <- as.character(Idents(object))
rename_cells <- rownames(object@meta.data)[grep("A",object$cell.types)]
object@meta.data[rename_cells,"cell.types"] = as.character(object@meta.data[rename_cells,"singler1main"])

cluster_5_4950_8525 <- rownames(object@meta.data)[object$integrated_snn_res.0.6 %in% 5 & 
                                                      object$groups %in% c("4950","8525")]
object@meta.data[cluster_5_4950_8525,"cell.types"] = "Epithelial cells"

table(object[["cell.types"]])
object[["cell.types"]] %>% table() %>% kable() %>% kable_styling()
object[["cell.types"]] %>% table() %>% prop.table() %>% kable() %>% kable_styling()

object <- AddMetaColor(object = object, label= "cell.types", colors = singler_colors2)
Idents(object) <- "cell.types"
object %<>% sortIdent()
TSNEPlot.1(object, group.by="cell.types",cols = ExtractMetaColor(object),
           label = T,pt.size = 1,no.legend = F,label.repel = T,
           label.size = 4, repel = T,do.return= T,do.print = T,alpha = 0.9,
           title = "All cell types in tSNE plot")

save(object,file="data/BladderCancer_mm10_6_20190726.Rda")
##############################
# split tsne plot
##############################
TSNEPlot.1(object, cols = ExtractMetaColor(object),label = F,pt.size = 1,
         split.by = "groups", group.by = "cell.types",label.size = 4, repel = T, 
         no.legend = T, do.print = T,border = T,
         ncol=3,title = "Compare cell types in each sample")
table(object$cell.types,object@meta.data[,'groups']) %>% kable() %>% kable_styling()

Idents(object) = "groups"
for (sample in c("4950", "8524", "8525")) {
        subset_object <- subset(object, idents = sample)
        Idents(subset_object) = "cell.types"
        p1 <- TSNEPlot.1(subset_object, group.by="cell.types",pt.size = 1,
                         label = F,no.legend = T,cols = ExtractMetaColor(subset_object),
                         label.repel = F, alpha = 1,label.size = 4, repel = T,border = T, 
                         title = "Total Clusters",do.print = F)
        
        p2 <- TSNEPlot.1(subset_object, group.by="cell.types",pt.size = 1,
                         label = F,no.legend = T,cols = ExtractMetaColor(subset_object),
                         label.repel = F, alpha = 1,border = T, split.by = "conditions",
                         label.size = 4, repel = T,title = NULL, do.print = F)
        
        jpeg(paste0(path,"S1_split_TSNEPlot_",sample,"_cell.types.jpeg"), units="in", width=10, height=7,res=600)
        print(plot_grid(p1, p2, align = "h")+
                      ggtitle(sample)+
                      theme(plot.title = element_text(hjust = 0.5,size = 18)))
        dev.off()
}

##############################
# split tsne plot and FeaturePlot
##############################
DefaultAssay(object) = "RNA"
Idents(object) = "cell.types"
object %<>% sortIdent()
TSNEPlot.1(object, group.by = "cell.types",split.by = "conditions",
           cols = ExtractMetaColor(object),label = T,label.size = 2, ncol = 1,
           width=3, height=7,do.print = T,border = T, no.legend = T,alpha = 1)

object@meta.data$project = "Mouse_tumor"
Idents(object) = "project"
FeaturePlot.1(object, features= c("Basal","Luminal","EMT_and_claudin",
                                "P53_like","neuroendocrine","Neuronal_differentiation"),
              ncol = 3,do.print = T,width=8, height=7, border = T)


Idents(object) = "groups"
for (sample in c("4950", "8524", "8525")) {
        subset_object <- subset(object, idents = sample)
        Idents(subset_object) = "cell.types"
        subset_object %<>% sortIdent()
        TSNEPlot.1(subset_object, group.by = "cell.types",split.by = "conditions",
                   cols = ExtractMetaColor(subset_object),label = T,label.size = 2, ncol = 1,
                   width=3, height=7,do.print = T,border = T, no.legend = T,alpha = 1,
                   unique.name = T)

        FeaturePlot.1(subset_object, features= c("Basal","Luminal","EMT_and_claudin",
                                          "P53_like","neuroendocrine","Neuronal_differentiation"),
                      ncol = 3,do.print = T,width=8, height=7, border = T,
                      unique.name = T)
}
