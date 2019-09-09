library(Seurat)
library(SingleR)
library(dplyr)
library(tidyr)
library(kableExtra)
source("../R/Seurat3_functions.R")
source("R/utils/SingleR_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# ======== 2.1 =========== test with known markers==================
(load(file="data/BladderCancer_mm10_6_20190726.Rda"))
DefaultAssay(object) <- "RNA"
df_markers <- readxl::read_excel("doc/Celltype Markers.xlsx")

#df_markers <- readxl::read_excel("../seurat_resources/bio-rad-markers.xlsx",sheet = "Human.sub")
colnames(df_markers) = gsub(" ","_",colnames(df_markers))
colnames(df_markers) = gsub(":|\\/","_",colnames(df_markers))
colnames(df_markers) = gsub("\\+","",colnames(df_markers))
markers = df_markers[,-grep("Alias",colnames(df_markers))]
marker.list <- df2list(markers)

GeneSets = read.delim("data/seurat_resources/gene_set_Bladder_Cancer.txt",row.names =1,header = F,
                      stringsAsFactors = F)
GeneSets %>% kable() %>% kable_styling()
marker.list <- df2list(t(GeneSets))

marker.list %<>% lapply(function(x) x[1:16]) %>% 
    lapply(function(x) FilterGenes(object,x)) %>% 
    lapply(function(x) x[!is.na(x)]) %>% 
    lapply(function(x) x[1:min(length(x),9)])
marker.list <- marker.list[!is.na(marker.list)]
marker.list %>% list2df %>% t %>% kable() %>% kable_styling()

Idents(object) = "integrated_snn_res.0.6"
for(i in 1:length(marker.list)){
    if(length(marker.list[[i]]) == 0) next
    p <- lapply(marker.list[[i]], function(marker) {
        FeaturePlot(object = object, feature = marker,pt.size = 0.5,
                    reduction="tsne", label = T)+
            NoLegend()+
            ggtitle(paste0(marker,Alias(df = df_markers,gene = marker)))+
            theme(plot.title = element_text(hjust = 0.5,size = 15,face = "plain"))
    })
    jpeg(paste0(path,"markers/",names(marker.list)[i],".jpeg"),units="in", width=10, height=7,res=600)
    print(do.call(cowplot::plot_grid, p)+ ggtitle(paste(names(marker.list)[i],"markers"))+
              theme(plot.title = element_text(hjust = 0.5,size = 20)))
    dev.off()
    print(paste0(i,":",length(marker.list)))
}

object_data <- object@assays$RNA@data
save(object_data, file = "data/object_data_mm10_6_20190706.Rda")

#====== 2.2 pathway analysis ==========================================
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

DefaultAssay(object) = "RNA"
for(i in 1:length(GeneSets_list)){
    object %<>% AddModuleScore(features = GeneSets_list[i], name = names(GeneSets_list[i]))
    svMisc::progress(i/length(GeneSets_list)*100)
}

colnames(object@meta.data)[colnames(object@meta.data) %in% paste0(colnames(GeneSets.df),"1")] =
    colnames(GeneSets.df)
colnames(object@meta.data) = gsub("_markers","",colnames(object@meta.data))

#=============== multiple color in the single Featureplot===================================
features.plot.list = list(c("Basal","Luminal"), 
                          c("Basal","EMT_and_claudin"),
                          c("Luminal","EMT_and_claudin"),
                          c("Basal","ITGA6"))
cols.use.list = list(c("#F8766D", "#00B0F6","#E31A1C"),
                     c("#F8766D", "#A3A500","#E31A1C"),
                     c("#00B0F6", "#A3A500","#E31A1C"),
                     c("#F8766D", "#00BF7D","#E31A1C"))
split.object <- SplitSeurat(object)
for(i in 1:4){
    jpeg(paste0(path,"Human_tsne_359_",i,".jpeg"), units="in", width=10, height=7,res=600)
    FeaturePlot(object = split.object[[2]], features.plot = features.plot.list[[i]], 
                cols.use = c("grey",cols.use.list[[i]]), overlay = TRUE, no.legend = FALSE)
    dev.off()
}

features.plot = c("Basal","EMT_and_claudin","Luminal")
p <- list()
for(i in 1:length(features.plot)){
    p[[i]] <- SingleFeaturePlot.1(object = object,
                                  feature = features.plot[i],
                                  threshold = 0.2)
}
jpeg(paste0(path,"Human_tsne.jpeg"), units="in", width=10, height=7,res=600)
print(do.call(plot_grid, p))
dev.off()

#=============== single color per Featureplot for gene sets ===================================
object_list <- SplitSeurat(object = object, split.by = "orig.ident")
gradient_list <- list(c("#c3d2e7", "#214069"),#blue
                      c("#f7e5b3","#8a6601"), #orange
                      c("#c1e2bf", "#1e601a"),#green
                      c("#f4a3a4","#880f10"),#red
                      c("#e7cdbe","#6a3518"), #brown
                      c("#f999cb","#90014c")) #purple
(levels <- object_list[[length(object_list)]])
for(j in 1:5){
    p <- list()
    for(i in 1:length(levels)){
        p[[i]] <- SingleFeaturePlot.1(object = object_list[[i]],
                                      feature = names(GeneSets_list[j]),
                                      gradient.use = gradient_list[[j]],
                                      title=levels[i], threshold = 0.2)
    }
    jpeg(paste0(path,names(GeneSets_list[j]),".jpeg"), units="in", width=10, height=7,res=600)
    print(do.call(plot_grid, p))
    dev.off()
}

for(j in 1:5/5){
    p <- list()
    for(i in 1:length(levels)){
        p[[i]] <- SingleFeaturePlot.1(object = object_list[[i]],
                                      feature = "Luminal",
                                      title=levels[i], threshold = j)
    }
    jpeg(paste0(path,"Luminal_",j,".jpeg"), units="in", width=10, height=7,res=600)
    print(do.call(plot_grid, p)+
              ggtitle(paste0("threshold = ",j))+
              theme(text = element_text(size=15),
                    plot.title = element_text(hjust = 0.5,size = 15, face = "bold")))
    dev.off()
}
#====== 2.2 marker gene analysis ==========================================
Blueprint_encode = read.csv("R/seurat_resources/Hpca_Blueprint_encode_main.csv",row.names =1,header = T,
                      stringsAsFactors = F)
marker.list <- df2list(Blueprint_encode)
marker.list <- lapply(marker.list, function(x) HumanGenes(object,x))

marker.list %>% list2df %>% head(15) %>% kable() %>% kable_styling()

.FeaturePlot <- function(x){
    p <- FeaturePlot(object = object, 
                     reduction.use = "tsne",
                     features.plot = x, min.cutoff = NA, do.return =T,
                     cols.use = c("lightgrey","blue"), pt.size = 2)
    return(p)
}

marker.list = "ITGA6"
for(i in 1:length(marker.list)){
    p <- .FeaturePlot(x = marker.list[[i]][1:9])+ ggtitle(names(marker.list)[i])+
        theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"))
    jpeg(paste0(path,names(marker.list)[i],".jpeg"), units="in", width=10, height=7,res=600)
    print(do.call(plot_grid, p))
    dev.off()
}

Epithelium <- HumanGenes(object,c("Epcam","KRT19","KRT5","Cdh1",
                                         "MUC1","Msln"))
p <- .FeaturePlot(x = Epithelium)
jpeg(paste0(path,"Epithelium",".jpeg"), units="in", width=10, height=7,res=600)
print(do.call(plot_grid, p))
dev.off()


jpeg(paste0(path,"ITGA6",".jpeg"), units="in", width=10, height=7,res=600)
print(.FeaturePlot(x = "ITGA6"))
dev.off()

HumanGenes(object,c("PTPRC"))
df_markers = data.frame("gene" = "PTPRC", "gene.alias" = "CD45")
SplitSingleFeaturePlot(object, alias = df_markers, 
                       group.by = "ident",split.by = "orig.ident",
                       no.legend = T,label.size=3,do.print =T,nrow = 1,
                       markers = "PTPRC", threshold = NULL)

hmarkers <- HumanGenes(object, c("APOBEC3A","APOBEC3B","APOBEC3D","APOBEC3C",
                                        "APOBEC3G","APOBEC3F","APOBEC3H","AICDA"))
df_markers = data.frame("gene" = "AICDA", "gene.alias" = "AID")
SplitSingleFeaturePlot(object, alias = df_markers, 
                       group.by = "ident",split.by = "orig.ident",
                       no.legend = T,label.size=3,do.print =T,nrow = 1,
                       markers = hmarkers, threshold = 0.1)
