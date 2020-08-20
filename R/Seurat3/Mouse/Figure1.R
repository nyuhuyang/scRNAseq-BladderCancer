library(Seurat)
library(reshape2)
library(pheatmap)
library(kableExtra)
library(dplyr)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# load data
(load(file = "data/BladderCancer_mm10_6_20190726.Rda"))

# Fig 1-A
Idents(object) = "cell.types"
TSNEPlot.1(object, group.by="cell.types",cols = ExtractMetaColor(object),
           label = T,pt.size = 1,no.legend = T,label.repel = T,
           label.size = 7, repel = T,do.return= F,do.print = T,alpha = 0.85,
           border = T,title = NULL,width=7, height=7)
# Fig 1-B
TSNEPlot.1(object, group.by="cell.types",split.by = "conditions",ncol=1,
           cols = ExtractMetaColor(object),
           label = T,pt.size = 1,no.legend = T,label.repel = T,
           strip.text.size = 30,
           label.size = 8, repel = T,do.return= F,do.print = T,alpha = 0.85,
           border = T,title = NULL,width=7, height=14)

colnames(object@meta.data)[grep("EMT_and_claudin",names(object@meta.data))] = "emt_claudin"
colnames(object@meta.data)[grep("neuronal\ndifferentiation",names(object@meta.data))] = "neuronal differentiation"
colnames(object@meta.data)[grep("P53_like",names(object@meta.data))] = "P53 like"

DefaultAssay(object) = "SCT"
FeaturePlot.1(object, features = c("Basal","Luminal","emt_claudin"),
              border = T, strip.text.size = 30, width=21, height=7, no.legend = T,
              threshold = 0.5,do.print = T,do.return = F,pt.size = 1, unique.name = F)
FeaturePlot.1(object, features = c("P53 like", "neuroendocrine",
                                   "neuronal differentiation"),
              border = T, strip.text.size = 30, width=21, height=7, no.legend = T,
              threshold = 0.05,do.print = T,do.return = F,pt.size = 1, unique.name = F)
colnames(object@meta.data)
# scale bar
jpeg(paste0(path,"scale_bar.jpeg"), units="in", width=10, height=7,res=600)
FeaturePlot(object, features = c("Basal","Luminal","emt_claudin",
                                 "P53 like", "neuroendocrine",
                                 "Neuronal_differentiation"),
            cols = c("lightgrey", "red"),ncol = 3)
dev.off()

