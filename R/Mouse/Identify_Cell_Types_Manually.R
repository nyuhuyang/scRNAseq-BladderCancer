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
(lnames = load(file="./output/BladderCancer_20181026.RData"))

gene_set = read.delim("./doc/gene_set_Bladder_Cancer.txt",row.names =1,header = F,
                   stringsAsFactors = F)
gene_set %>% kable() %>% kable_styling()
gene_set.df <- as.data.frame(t(gene_set))

gene_set.list <- df2list(gene_set.df)
gene_set_list <- lapply(gene_set.list, function(x) MouseGenes(BladderCancer,x))

for(i in 1:length(gene_set_list)){
    BladderCancer <- .AddModuleScore(BladderCancer, genes.list = gene_set_list[i],
                                   ctrl.size = 5,enrich.name = names(gene_set_list[i]))
}

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
                                      title=levels[i], threshold = 0.1)
    }
    jpeg(paste0(path,names(gene_set_list[j]),".jpeg"), units="in", width=10, height=7,res=600)
    print(do.call(plot_grid, p))
    dev.off()
}

#====== 2.2 marker gene analysis ==========================================
immgen_main = read.csv("../SingleR/output/immgen_main.csv",row.names =1,header = T,
                      stringsAsFactors = F)
marker.list <- df2list(immgen_main)
marker.list <- lapply(marker.list, function(x) MouseGenes(BladderCancer,x))

marker.list %>% list2df %>% head(15) %>% kable() %>% kable_styling()

.FeaturePlot <- function(x){
    p <- FeaturePlot(object = BladderCancer, 
                     reduction.use = "tsne",
                     features.plot = x, min.cutoff = NA, do.return =T,
                     cols.use = c("lightgrey","blue"), pt.size = 0.5)
    return(p)
}

for(i in 1:length(marker.list)){
    p <- .FeaturePlot(x = marker.list[[i]][1:9])
    jpeg(paste0(path,names(marker.list)[i],".jpeg"), units="in", width=10, height=7,res=600)
    print(do.call(plot_grid, p))
    dev.off()
}

Epithelium <- MouseGenes(BladderCancer,c("Epcam","KRT19","KRT5","Cdh1",
                                         "MUC1","Msln"))
p <- .FeaturePlot(x = Epithelium)
jpeg(paste0(path,"Epithelium",".jpeg"), units="in", width=10, height=7,res=600)
print(do.call(plot_grid, p))
dev.off()
