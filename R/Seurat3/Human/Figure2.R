library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
marker_path <- paste0(path,"markers/")
if(!dir.exists(marker_path))dir.create(marker_path, recursive = T)

(load(file="data/BladderCancer_H2_20181109.Rda"))
object <- UpdateSeuratObject(BladderCancer)
remove(BladderCancer);GC()
DefaultAssay(object) <- "RNA"

GeneSets = read.delim("data/seurat_resources/gene_set_Bladder_Cancer.txt",row.names =1,header = F,
                      stringsAsFactors = F)
GeneSets %>% kable() %>% kable_styling()
GeneSets.df <- as.data.frame(t(GeneSets))
GeneSetNames <- c("Luminal","EMT_and_smooth_muscle","EMT_and_claudin",
                  "Basal","Squamous","Immune",
                  "Neuronal_differentiation","Down_regulated_in_CIS",
                  "Up_regulated_in_CIS")

GeneSets.df <- GeneSets.df[,GeneSetNames]
GeneSets.list <- df2list(GeneSets.df)
GeneSets_list <- lapply(GeneSets.list, function(x) FilterGenes(object,x))
GeneSets_list %>% list2df %>% t %>% kable() %>% kable_styling()
colnames(object@meta.data) = gsub("_markers","",colnames(object@meta.data))
object@meta.data = object@meta.data[,-which(colnames(object@meta.data) %in% GeneSetNames)]

DefaultAssay(object) = "RNA"
for(i in 1:length(GeneSets_list)){
        object %<>% AddModuleScore(features = GeneSets_list[i], name = names(GeneSets_list[i]))
        svMisc::progress(i/length(GeneSets_list)*100)
}

colnames(object@meta.data)[colnames(object@meta.data) %in% paste0(colnames(GeneSets.df),"1")] =
        colnames(GeneSets.df)
colnames(object@meta.data) = gsub("_markers","",colnames(object@meta.data))
GeneSetNames = c("Luminal", "EMT_and_smooth_muscle","EMT_and_claudin","Basal")
cols <- list(c("#c3d2e7", "#214069"),#blue
              c("#f7e5b3","#8a6601"), #orange
              c("#c1e2bf", "#1e601a"),#green
              c("#f4a3a4","#880f10"))#red

for(i in seq_along(GeneSetNames)){
        p <- FeaturePlot.1(object = object, feature = GeneSetNames[i],
                           pt.size = 0.5,split.by = "orig.ident",
                           col = cols[[i]],
                         reduction="tsne", label = F)+
                ggtitle(GeneSetNames[i])+
                theme(plot.title = element_text(hjust = 0.5,size = 15,face = "plain"),
                      panel.border = element_rect(colour = "black"))
        
        jpeg(paste0(path,GeneSetNames[i],".jpeg"),units="in", width=10, height=7,res=600)
        print(p)
        dev.off()
        svMisc::progress(i/length(GeneSetNames)*100)
}
