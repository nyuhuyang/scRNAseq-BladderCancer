library(SingleR)
library(Seurat)
library(reshape2)
library(pheatmap)
library(kableExtra)
library(dplyr)
source("R/utils/Seurat_functions.R")
source("R/utils/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
lname1 = load(file = "./data/BladderCancer_H2_20181109.Rda");lname1
singler = CreateSinglerObject(BladderCancer@data, annot = NULL, project.name=BladderCancer@project.name,
                              min.genes = 500,technology = "10X", species = "Human", citation = "",
                              normalize.gene.length = F, variable.genes = "de",
                              fine.tune = F, do.signatures = F, clusters = NULL)
# singler didn't find all cell labels
length(singler$singler[[1]]$SingleR.single$labels) == ncol(BladderCancer@data)
singler$meta.data$orig.ident = BladderCancer@meta.data$orig.ident # the original identities, if not supplied in 'annot'
singler$meta.data$xy = BladderCancer@dr$tsne@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = BladderCancer@ident # the Seurat clusters (if 'clusters' not provided)
save(singler,file="./output/singler_BladderCancer_H2_20181109.RData")
#====== 3.2 SingleR specifications ==========================================
# Step 1: Spearman coefficient
lnames = load(file = "./output/singler_BladderCancer_H2_20181109.RData")
attach(immgen)
attach(mouse.rnaseq)

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

#Or by all cell types (showing the top 50 cell types):
jpeg(paste0(path,"DrawHeatmap_norm_sub1.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single, top.n = 50,normalize = F))
dev.off()

jpeg(paste0(path,"DrawHeatmap_norm_sub2.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single,top.n = 50,normalize = F))
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

kable(table(singler$singler[[1]]$SingleR.single$labels,
            singler$meta.data$orig.ident)) %>%
        kable_styling()
kable(table(singler$singler[[2]]$SingleR.single$labels,
            singler$meta.data$orig.ident)) %>%
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
BladderCancer <- SetAllIdent(object = BladderCancer, id = "singler2sub")

BladderCancer = AddMetaColor(object = BladderCancer, colors = singler_colors1[1:19], 
                         label = "singler2sub")

##############################
# draw tsne plot
##############################
p5 <- DimPlot.1(object = BladderCancer, reduction.use = "tsne", dim.1 = 1, dim.2 = 2, 
          group.by = "ident", do.return = TRUE,
          cols.use = ExtractMetaColor(object = BladderCancer), no.legend = T,
          do.label =T,label.size= 4, label.repel = T,force=1)+
        ggtitle("Supervised cell type labeling by Blueprint+ENCOD")+
        theme(text = element_text(size=15),
              plot.title = element_text(hjust = 0.5))

jpeg(paste0(path,"PlotTsne_sub2.jpeg"), units="in", width=10, height=7,
     res=600)
print(p5)
dev.off()

save(BladderCancer,file = "./data/BladderCancer_H2_20181109.Rda")
##############################
# subset Seurat
###############################
table(BladderCancer@meta.data$orig.ident)
table(BladderCancer@ident)

p4 <- SplitTSNEPlot(object = BladderCancer, split.by = "orig.ident",
                    do.label = T, group.by = "ident", 
                    do.return = TRUE, no.legend = T, do.print =T,
                    pt.size = 1,label.size = 3, label.repel = T,force=1)