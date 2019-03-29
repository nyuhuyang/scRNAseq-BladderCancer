library(Seurat)
library(SingleR)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
source("R/utils/Seurat_functions.R")
source("R/utils/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
#====== 2.1 pathway analysis ==========================================
(load(file="./data/BladderCancer_M2_20181026.Rda")

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
#=============== multiple color in the single Featureplot===================================
features.plot.list = list(c("Basal_markers","Luminal_markers"), 
                          c("Basal_markers","EMT_and_claudin_markers"),
                          c("Luminal_markers","EMT_and_claudin_markers"),
                          c("Basal_markers","Itga6"))
cols.use.list = list(c("#F8766D", "#00B0F6","#E31A1C"),
                     c("#F8766D", "#A3A500","#E31A1C"),
                     c("#00B0F6", "#A3A500","#E31A1C"),
                     c("#F8766D", "#00BF7D","#E31A1C"))
for(i in 1:4){
    jpeg(paste0(path,"Mouse_tsne_",i,".jpeg"), units="in", width=10, height=7,res=600)
    FeaturePlot(object = BladderCancer, features.plot = features.plot.list[[i]], 
                cols.use = c("grey",cols.use.list[[i]]), overlay = TRUE, no.legend = FALSE)
    dev.off()
}
#=============== single color per Featureplot===================================
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
immgen_main = read.csv("R/seurat_resources/immgen_main.csv",row.names =1,header = T,
                      stringsAsFactors = F)
marker.list <- df2list(immgen_main)
marker.list <- lapply(marker.list, function(x) MouseGenes(BladderCancer,x))

marker.list %>% list2df %>% head(15) %>% kable() %>% kable_styling()

.FeaturePlot <- function(x){
    p <- FeaturePlot(object = BladderCancer, 
                     reduction.use = "tsne",
                     features.plot = x, min.cutoff = NA, do.return =T,
                     cols.use = c("lightgrey","blue"), pt.size = 2)
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

jpeg(paste0(path,"m_ITGA6",".jpeg"), units="in", width=10, height=7,res=600)
print(.FeaturePlot(x = "Itga6"))
dev.off()

# pull marker genes from panel
mouse.rnaseq_main = read.csv("R/seurat_resources/mouse.rnaseq_main.csv",
                             row.names =1, header = T, stringsAsFactors = F)
immgen_main = read.csv("R/seurat_resources/immgen_main.csv",
                             row.names =1, header = T, stringsAsFactors = F) 

immgen_main %<>% lapply(function(x) Hmisc::capitalize(tolower(x)))

(intersect(mouse.rnaseq_main$T_cells, immgen_main$T_cells))[1:10]

#===== for Jan 15, 2019 email========
library(biomaRt)
human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
searchAttributes(mouse, "mgi_symbol")
searchAttributes(human, "hgnc_symbol")

hmarkers <-c("APOBEC3A","APOBEC3B","APOBEC3D","APOBEC3C",
            "APOBEC3G","APOBEC3F","APOBEC3H","AID")
(markers = biomaRt::getLDS(attributes = c("hgnc_symbol"), 
                           filters = "hgnc_symbol", 
                           values = hmarkers, mart = human, 
                           attributesL = c("mgi_symbol"), martL = mouse))


mmarkers <-c("Apobec1","Apobec3","Apobec2")
(df_markers = biomaRt::getLDS(attributes = c("mgi_symbol"), 
                           filters = "mgi_symbol", 
                           values = mmarkers, mart = mouse, 
                           attributesL = c("hgnc_symbol"), martL = human))
colnames(df_markers) = c("markers","markers.Alias")

SplitSingleFeaturePlot(BladderCancer, 
                       #alias = Alias(df = df_markers, gene = marker), 
                       group.by = "ident",split.by = "orig.ident",
                       no.legend = T,label.size=3,do.print =T,nrow = 1,
                       markers = df_markers$markers, threshold = NULL)

