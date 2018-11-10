library(SingleR)
library(Seurat)
library(reshape2)
library(pheatmap)
library(kableExtra)
library(dplyr)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
lname1 = load(file = "./data/BladderCancer_20181026.Rda");lname1
singler = CreateSinglerObject(BladderCancer@data, annot = NULL, project.name=BladderCancer@project.name,
                              min.genes = 500,technology = "10X", species = "Mouse", citation = "",
                              normalize.gene.length = F, variable.genes = "de",
                              fine.tune = F, do.signatures = F, clusters = NULL)
# singler didn't find all cell labels
length(singler$singler[[1]]$SingleR.single$labels) == ncol(BladderCancer@data)
singler$meta.data$orig.ident = BladderCancer@meta.data$orig.ident # the original identities, if not supplied in 'annot'
singler$meta.data$xy = BladderCancer@dr$tsne@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = BladderCancer@ident # the Seurat clusters (if 'clusters' not provided)
save(singler,file="./output/singler_BladderCancer_20181026.RData")
#====== 3.2 SingleR specifications ==========================================
# Step 1: Spearman coefficient
lnames = load(file = "./output/singler_BladderCancer_20181026.RData")
lnames = load(file = "../SingleR/data/immgen.rda") 
lnames = load(file = "../SingleR/data/mouse.rnaseq.rda") 

singler$seurat = BladderCancer
SingleR.DrawScatter(sc_data = singler$seurat@data,cell_id = 10, 
                    ref = Hpca, sample_id = 232)

# Step 2: Multiple correlation coefficients per cell types are aggregated 

# to provide a single value per cell type per single-cell. 
# In the examples below we use the 80% percentile of correlation values.
# for visualization purposes we only present a subset of cell types (defined in labels.use)
out = SingleR.DrawBoxPlot(sc_data = singler$seurat@data,cell_id = 10, 
                          ref = immgen,main_types = T,
                          labels.use=c('B cells','T cells','DC','Macrophages','Monocytes','NK cells',
                                       'Mast cells','Neutrophils','Fibroblasts','Endothelial cells'))
print(out$plot)

# Step 3: DrawHeatmap

SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single.main, top.n = Inf)
#Or by all cell types (showing the top 50 cell types):
jpeg(paste0(path,"DrawHeatmap.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single.main, top.n = 50))
dev.off()
jpeg(paste0(path,"DrawHeatmap_normF.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single.main,top.n = 50,normalize = F))
dev.off()
#Next, we can use the fine-tuned labels to color the t-SNE plot:

# Step 4: DrawHeatmap
# ImmueGene-------
ncolor = length(unique(singler$singler[[1]]$SingleR.single$labels));ncolor
out1 = SingleR.PlotTsne.1(singler$singler[[1]]$SingleR.single,
                          singler$meta.data$xy,do.label=F,
                          do.letters = F,labels = singler$singler[[1]]$SingleR.single$labels,
                          label.size = 5, dot.size = 2,do.legend = T,alpha = 1,
                          label.repel = F,force=2,
                          title = "Supervised sub cell type labeling by HPCA",
                          colors = singler.colors[1:ncolor])
jpeg(paste0(path,"/PlotTsne_sub1.jpeg"), units="in", width=10, height=7,
     res=600)
out1
dev.off()
# Mouse.rnaseq-------
out2 = SingleR.PlotTsne.1(singler$singler[[2]]$SingleR.single,
                         singler$meta.data$xy,do.label=T,
                         do.letters = F,labels = singler$singler[[2]]$SingleR.single$labels,
                         label.size = 5, dot.size = 2,do.legend = F,alpha = 1,
                         label.repel = T,force=4,
                         title = "Supervised sub cell type labeling by Blueprint+ENCOD",)
#label.repel = T,force=2)
jpeg(paste0(path,"/PlotTsne_sub2.jpeg"), units="in", width=10, height=7,
     res=600)
out2
dev.off()

#Finally, we can also view the labeling as a table compared to the original identities:

kable(table(singler$singler[[2]]$SingleR.single$labels,
            singler$meta.data$orig.ident)) %>%
        kable_styling()
kable(table(singler$meta.data$orig.ident,
            singler$singler[[1]]$SingleR.single.main$labels)) %>%
        kable_styling()

##############################
# add singleR label to Seurat
###############################
singler_colors <- readxl::read_excel("./doc/singler.colors.xlsx")
singler_colors1 = as.vector(singler_colors$singler.color1[!is.na(singler_colors$singler.color1)])
singler_colors2 = as.vector(singler_colors$singler.color2[!is.na(singler_colors$singler.color2)])
length(singler_colors1);length(singler_colors2)

singlerDF = data.frame("singler1sub"=singler$singler[[1]]$SingleR.single$labels,
                       "singler1main"=singler$singler[[1]]$SingleR.single.main$labels,
                       "singler2sub"=singler$singler[[2]]$SingleR.single$labels,
                       "singler2main"=singler$singler[[2]]$SingleR.single.main$labels,
                       "cell.names" = rownames(singler$singler[[1]]$SingleR.single$labels))
knowDF = data.frame("cell.names"= BladderCancer@cell.names)
ident.DF = full_join(singlerDF,knowDF, by="cell.names")
ident.DF<- apply(ident.DF,2,as.character)
rownames(ident.DF) = ident.DF[,"cell.names"]
ident.DF = ident.DF[,-which(colnames(ident.DF) == "cell.names")]
apply(ident.DF,2,function(x) length(unique(x)))

BladderCancer <- AddMetaData(object = BladderCancer,
                   metadata = as.data.frame(ident.DF))
BladderCancer <- SetAllIdent(object = BladderCancer, id = "singler1main")
##############################
# process singler.color
##############################

#' AddMetaColor: prepare meta.data data frame to store color code
#' @object Seurat object
#' @colors vector of hex colors
#' @label ident label
#' @example 
# singlerDF = AddMetaColor(object = BladderCancer, colors = singler_colors1[1:16], 
#                         label = "singler2sub")
AddMetaColor <- function(object, colors, label=NULL){
        df_ident <- data.frame("index" =  as.numeric(as.factor(as.vector(object@ident))),
                               "ident" =  as.vector(object@ident),
                               "cell.names" = object@cell.names,
                               row.names = object@cell.names)
        df_colors = data.frame("colors" = colors[1:length(unique(df_ident$index))],
                               "index" = 1:length(unique(df_ident$index)))
        ident_colors <- full_join(df_ident, df_colors, by = "index")
        ident_colors = ident_colors[!is.na(ident_colors$cell.names),]
        rownames(ident_colors) = ident_colors$cell.names

        DF_colors <- data.frame("colors" = ident_colors$colors,
                                row.names = ident_colors$cell.names,
                                stringsAsFactors = F)
        if(!is.null(label)) colnames(DF_colors)[1] = paste(colnames(DF_colors)[1],label,
                                                            sep = ".")
        return(DF_colors)
}

singlerDF1 = AddMetaColor(object = BladderCancer, colors = singler_colors1[1:16], 
                         label = "singler1main")
BladderCancer <- AddMetaData(object = BladderCancer, metadata = singlerDF1)

##############################
# draw tsne plot
##############################

#' ExtractMetaColor: extract color code from meta.data
#' @object seurat object, meta.data slot must have "colors.singlerxxx"
#' @example 
# singlerDF = factor2color(mat = singler$singler[[2]]$SingleR.single$labels,
#                       colors = singler_colors[3:26], col.name = "singler2sub")
ExtractMetaColor <- function(object){
        meta.data = object@meta.data
        colors.id = grep("colors", colnames(meta.data),value = T)
        id = sub("colors.","",colors.id)
        meta.data = meta.data[,c(id,colors.id)]
        meta.data$index <- as.numeric(as.factor(meta.data[,id]))
        df_colors = meta.data[!duplicated(meta.data$index),]
        df_colors = df_colors[order(df_colors$index),]

        return(as.character(df_colors$colors))
}

df_colors <- ExtractMetaColor(object = BladderCancer)
p3 <- DimPlot.1(object = BladderCancer, reduction.use = "tsne", dim.1 = 1, dim.2 = 2, 
          group.by = "ident", do.return = TRUE,
          cols.use = df_colors, no.legend = T,
          do.label =T,label.size=4, label.repel = T,force=3)+
        ggtitle("Supervised main cell type labeling by ImmueGene")+
        theme(text = element_text(size=15),
              plot.title = element_text(hjust = 0.5))

jpeg(paste0(path,"PlotTsne_main1~.jpeg"), units="in", width=10, height=7,
     res=600)
print(p3)
dev.off()

save(BladderCancer,file="./output/BladderCancer_20181026.RData")
##############################
# subset Seurat
###############################
table(BladderCancer@meta.data$orig.ident)
table(BladderCancer@ident)

p4 <- SplitTSNEPlot(BladderCancer,group.by = "ident", split.by = "orig.ident",
                    no.legend = T,do.label =T,label.size=3,
                    return.plots =T, label.repel = T,force=1)
jpeg(paste0(path,"splitetsne.jpeg"), units="in", width=10, height=7,
     res=600)
do.call(plot_grid, p4)
dev.off()