library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
source("R/utils/Seurat3_functions.R")
source("R/utils/FeaturePlot.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
marker_path <- paste0(path,"markers/")
if(!dir.exists(marker_path))dir.create(marker_path, recursive = T)

(load(file="data/BladderCancer_H2_20181109.Rda"))
object <- UpdateSeuratObject(BladderCancer)
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

Idents(object) = "integrated_snn_res.0.6"
for(i in 1:length(GeneSets_list)){
    p <- FeaturePlot(object = object, feature = GeneSetNames[i],pt.size = 0.5,
                    reduction="tsne", label = T)+
        ggtitle(GeneSetNames[i])+
        theme(plot.title = element_text(hjust = 0.5,size = 15,face = "plain"),
              panel.border = element_rect(colour = "black"))
    
    jpeg(paste0(path,"markers/",GeneSetNames[i],".jpeg"),units="in", width=10, height=7,res=600)
    print(p)
    dev.off()
    print(paste0(i,":",length(GeneSets_list)))
}

#=============== multiple color in the single Featureplot===================================
features.list = list(c("Basal","Luminal"),
                     c("Luminal","EMT_and_claudin"),
                     c("Basal","EMT_and_claudin"))
cols.use.list = list(c("#b88801", "#2c568c","#E31A1C"),
                     c("#2c568c", "#b88801","#E31A1C"),
                     c("#b88801", "#54307b","#E31A1C"))

for(i in 1:length(features.list)){
    g <- FeaturePlot.2(object = object, features = features.list[[i]],do.return = T,
                       overlay = T,cols = c("#d8d8d8",cols.use.list[[i]]),pt.size = 0.8)
    jpeg(paste0(path,"Mouse_tsne_",paste(features.list[[i]],collapse = "_"),".jpeg"), 
         units="in", width=10, height=7,res=600)
    print(g+theme(plot.title = element_text(hjust = 0.5,size = 20),
                  legend.position="bottom",
                  legend.text=element_text(size=15))+
              guides(colour = guide_legend(override.aes = list(size=5)), 
                     shape = guide_legend(override.aes = list(size=5))))
    dev.off()
}


#========== GenePlot ============
feature1 <- FilterGenes(object,c("KRT5","KRT14","CLDN3","CLDN4", "CLDN7"))
feature2 <- c("KRT8","KRT8","VIM","VIM","VIM")

for(i in 1:2){
    jpeg(paste0(path,feature1[i],"+",feature2[i],".jpeg"), 
         units="in", width=10, height=7,res=600)
    print(FeatureScatter(object, feature1 = feature1[i], feature2 = feature2[i]))
    dev.off()
    
}

KRT5_KRT8_p <- subset(object, subset = KRT5 > 0 & KRT8 > 0)
FeatureScatter(KRT5_KRT8_p, feature1 = feature1[i], feature2 = feature2[i])

KRT5_KRT8_n <- subset(object, subset = KRT5 == 0 & KRT8 == 0)
FeatureScatter(KRT5_KRT8_n, feature1 = feature1[i], feature2 = feature2[i])

for(i in 3:5){
    jpeg(paste0(path,feature1[i],"+",feature2[i]," KRT5_KRT8_p.jpeg"), 
         units="in", width=10, height=7,res=600)
    print(FeatureScatter(KRT5_KRT8_p, feature1 = feature1[i], feature2 = feature2[i]))
    dev.off()
    
}

for(i in 3:5){
    jpeg(paste0(path,feature1[i],"+",feature2[i]," KRT5_KRT8_n.jpeg"), 
         units="in", width=10, height=7,res=600)
    print(FeatureScatter(KRT5_KRT8_n, feature1 = feature1[i], feature2 = feature2[i]))
    dev.off()
    
}
