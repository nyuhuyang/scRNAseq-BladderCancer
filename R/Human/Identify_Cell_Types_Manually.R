library(Seurat)
library(SingleR)
library(dplyr)
library(tidyr)
library(kableExtra)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
#====== 2.1 pathway analysis ==========================================
(lnames = load(file="./output/H_BladderCancer_CCA_20181109.Rda"))

gene_set = read.delim("./doc/gene_set_Bladder_Cancer.txt",row.names =1,header = F,
                   stringsAsFactors = F)
gene_set %>% kable() %>% kable_styling()
gene_set.df <- as.data.frame(t(gene_set))

gene_set.list <- df2list(gene_set.df)
gene_set_list <- lapply(gene_set.list, function(x) HumanGenes(BladderCancer,x))

for(i in 1:length(gene_set_list)){
    BladderCancer <- .AddModuleScore(BladderCancer, genes.list = gene_set_list[i],
                                   ctrl.size = 5,enrich.name = names(gene_set_list[i]),
                                   only.pos = F)
}
#=============== multiple color in the single Featureplot===================================
features.plot.list = list(c("Basal_markers","Luminal_markers"), 
                          c("Basal_markers","EMT_and_claudin_markers"),
                          c("Luminal_markers","EMT_and_claudin_markers"),
                          c("Basal_markers","ITGA6"))
cols.use.list = list(c("#F8766D", "#00B0F6","#E31A1C"),
                     c("#F8766D", "#A3A500","#E31A1C"),
                     c("#00B0F6", "#A3A500","#E31A1C"),
                     c("#F8766D", "#00BF7D","#E31A1C"))
for(i in 1:4){
    jpeg(paste0(path,"Human_tsne_",i,".jpeg"), units="in", width=10, height=7,res=600)
    FeaturePlot(object = BladderCancer, features.plot = features.plot.list[[i]], 
                cols.use = c("grey",cols.use.list[[i]]), overlay = TRUE, no.legend = FALSE)
    dev.off()
}

features.plot = c("Basal_markers","EMT_and_claudin_markers","Luminal_markers")
p <- list()
for(i in 1:length(features.plot)){
    p[[i]] <- SingleFeaturePlot.1(object = BladderCancer,
                                  feature = features.plot[i],
                                  threshold = 0.2)
}
jpeg(paste0(path,"Human_tsne.jpeg"), units="in", width=10, height=7,res=600)
print(do.call(plot_grid, p))
dev.off()

#=============== single color per Featureplot for gene sets ===================================
BladderCancer_list <- SplitSeurat(object = BladderCancer, split.by = "orig.ident")
gradient_list <- list(c("#c3d2e7", "#214069"),#blue
                      c("#f7e5b3","#8a6601"), #orange
                      c("#c1e2bf", "#1e601a"),#green
                      c("#f4a3a4","#880f10"),#red
                      c("#e7cdbe","#6a3518"), #brown
                      c("#f999cb","#90014c")) #purple
(levels <- BladderCancer_list[[length(BladderCancer_list)]])
for(j in 1:5){
    p <- list()
    for(i in 1:length(levels)){
        p[[i]] <- SingleFeaturePlot.1(object = BladderCancer_list[[i]],
                                      feature = names(gene_set_list[j]),
                                      gradient.use = gradient_list[[j]],
                                      title=levels[i], threshold = 0.2)
    }
    jpeg(paste0(path,names(gene_set_list[j]),".jpeg"), units="in", width=10, height=7,res=600)
    print(do.call(plot_grid, p))
    dev.off()
}

for(j in 1:5/5){
    p <- list()
    for(i in 1:length(levels)){
        p[[i]] <- SingleFeaturePlot.1(object = BladderCancer_list[[i]],
                                      feature = "Luminal_markers",
                                      title=levels[i], threshold = j)
    }
    jpeg(paste0(path,"Luminal_markers_",j,".jpeg"), units="in", width=10, height=7,res=600)
    print(do.call(plot_grid, p)+
              ggtitle(paste0("threshold = ",j))+
              theme(text = element_text(size=15),
                    plot.title = element_text(hjust = 0.5,size = 15, face = "bold")))
    dev.off()
}
#====== 2.2 marker gene analysis ==========================================
Blueprint_encode = read.csv("../SingleR/output/Hpca_Blueprint_encode_main.csv",row.names =1,header = T,
                      stringsAsFactors = F)
marker.list <- df2list(Blueprint_encode)
marker.list <- lapply(marker.list, function(x) HumanGenes(BladderCancer,x))

marker.list %>% list2df %>% head(15) %>% kable() %>% kable_styling()

.FeaturePlot <- function(x){
    p <- FeaturePlot(object = BladderCancer, 
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

Epithelium <- HumanGenes(BladderCancer,c("Epcam","KRT19","KRT5","Cdh1",
                                         "MUC1","Msln"))
p <- .FeaturePlot(x = Epithelium)
jpeg(paste0(path,"Epithelium",".jpeg"), units="in", width=10, height=7,res=600)
print(do.call(plot_grid, p))
dev.off()


jpeg(paste0(path,"ITGA6",".jpeg"), units="in", width=10, height=7,res=600)
print(.FeaturePlot(x = "ITGA6"))
dev.off()
