########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(kableExtra)
source("R/utils/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
########################################################################
#
#  1 Seurat Alignment 
# 
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# Load the mouse.eyes dataset

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here
df_samples <- readxl::read_excel("doc/181008_Single_cell_sample list.xlsx")
sample_n <- df_samples$species == "Human"
samples <- df_samples$samples[sample_n]
projects <- df_samples$projects[sample_n]
conditions <- df_samples$conditions[sample_n]
species <- df_samples$species[sample_n]
df_samples %>% kable() %>% kable_styling()

BladderCancer_raw <- list()
BladderCancer_Seurat <- list()
for(i in 1:length(samples)){
    BladderCancer_raw[[i]] <- Read10X(data.dir = paste0("./data/",
                                             samples[i],"/outs/filtered_gene_bc_matrices/hg19/"))
    colnames(BladderCancer_raw[[i]]) <- paste0(samples[i],
                                            "_",colnames(BladderCancer_raw[[i]]))
    BladderCancer_Seurat[[i]] <- CreateSeuratObject(BladderCancer_raw[[i]],
                                                 min.cells = 3,
                                                 min.genes = 0,
                                                 project = projects[i],
                                                 names.delim = "_")
    BladderCancer_Seurat[[i]]@meta.data$conditions <- conditions[i]
}
#======1.1.2 QC before merge =========================
cell.number <- sapply(BladderCancer_Seurat, function(x) length(x@cell.names))
QC_list <- lapply(BladderCancer_Seurat, function(x) as.matrix(x = x@raw.data))
median.nUMI <- sapply(QC_list, function(x) median(colSums(x)))
median.nGene <- sapply(QC_list, function(x) median(apply(x,2,function(y) sum(length(y[y>0])))))

min.nUMI <- sapply(QC_list, function(x) min(colSums(x)))
min.nGene <- sapply(QC_list, function(x) min(apply(x,2,function(y) sum(length(y[y>0])))))

QC.list <- cbind(df_samples[sample_n,],cell.number, median.nUMI,median.nGene,min.nUMI,min.nGene,
                 row.names = samples)
write.csv(QC.list,paste0(path,"QC_list.csv"))
QC.list %>% kable() %>% kable_styling()
remove(QC_list,QC.list);GC()

#========1.1.3 merge ===================================
BladderCancer <- Reduce(function(x, y) MergeSeurat(x, y, do.normalize = F), BladderCancer_Seurat)
remove(BladderCancer_raw,BladderCancer_Seurat);GC()

mito.genes <- grep(pattern = "^MT-", x = rownames(x = BladderCancer@data), value = TRUE)
percent.mito <- Matrix::colSums(BladderCancer@raw.data[mito.genes, ])/Matrix::colSums(BladderCancer@raw.data)
BladderCancer <- AddMetaData(object = BladderCancer, metadata = percent.mito, col.name = "percent.mito")

BladderCancer@ident = factor(BladderCancer@ident,levels = samples)

g1 <- VlnPlot(object = BladderCancer, features.plot = c("nGene", "nUMI", "percent.mito"), 
              nCol = 1,point.size.use = 0.2,group.by= "orig.ident",
              x.lab.rot = T, do.return = T,return.plotlist =T)

remove(BladderCancer);GC()
#======1.2 load  SingleCellExperiment and normalizing the data=========================
# 1.2.1 Calculate median UMI per cell
Iname = load(file = "data/sce_list_Human_20181109.Rda")
seurat_list <- lapply(sce_list, as.seurat) %>%
                lapply(NormalizeData) %>%
                lapply(ScaleData) %>%
                lapply(FindVariableGenes, do.plot = FALSE)

for(i in 1:length(samples)){
    seurat_list[[i]]@meta.data$conditions <- conditions[i]
    seurat_list[[i]]@meta.data$species <- species[i]
}
# we will take the union of the top 1k variable genes in each dataset for
# alignment note that we use 1k genes in the manuscript examples, you can
# try this here with negligible changes to the overall results
g <- lapply(seurat_list, function(x) head(rownames(x@hvg.info), 2000))
genes.use <- unique(unlist(g))
for(i in 1:length(conditions)){
    genes.use <- intersect(genes.use, rownames(seurat_list[[i]]@data))
}
length(genes.use)

#  Perform a canonical correlation analysis (CCA) =========================
# run a canonical correlation analysis to identify common sources
# of variation between the two datasets.
remove(seurat_list,sce_list)
GC();GC();GC();GC();GC();GC();GC();GC();GC();GC()
BladderCancer <- RunCCA(seurat_list[[1]],seurat_list[[2]], 
                   genes.use = genes.use,
                   num.cc = 30)
save(BladderCancer, file = "./data/BladderCancer_H2_20181109.Rda")
# 1.2.3 calculate mitochondria percentage
mito.genes <- grep(pattern = "^MT-", x = rownames(x = BladderCancer@data), value = TRUE)
percent.mito <- Matrix::colSums(BladderCancer@raw.data[mito.genes, ])/Matrix::colSums(BladderCancer@raw.data)
BladderCancer <- AddMetaData(object = BladderCancer, metadata = percent.mito, col.name = "percent.mito")

BladderCancer@ident = factor(BladderCancer@ident,levels = samples)

g1 <- VlnPlot(object = BladderCancer, features.plot = c("nGene", "nUMI", "percent.mito"), 
              nCol = 1,point.size.use = 0.2,
              x.lab.rot = T, do.return = T,return.plotlist =T)

BladderCancer <- FilterCells(object = BladderCancer, subset.names = c("nGene","nUMI","percent.mito"),
                    low.thresholds = c(800,1500, -Inf), 
                    high.thresholds = c(Inf,Inf, 0.15))

g2 <- VlnPlot(object = BladderCancer, features.plot = c("nGene", "nUMI", "percent.mito"), 
              nCol = 1,point.size.use = 0.2,
              x.lab.rot = T, do.return = T,return.plotlist =T)

jpeg(paste0(path,"/S1_nGene.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[1]]+ggtitle("nGene in raw data")+ 
                    scale_y_log10(limits = c(200,10000)),#+ylim(c(0,1000)),
                g2[[1]]+ggtitle("nGene after filteration")+ 
                    scale_y_log10(limits = c(200,10000))))
dev.off()
jpeg(paste0(path,"/S1_nUMI.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[2]]+ggtitle("nUMI in raw data")+ 
                    scale_y_log10(limits = c(1000,100000)),#+ylim(c(0,1000)),
                g2[[2]]+ggtitle("nUMI after filteration")+ 
                    scale_y_log10(limits = c(1000,100000))))
dev.off()
jpeg(paste0(path,"/S1_mito.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[3]]+ggtitle("mito % in raw data")+ 
                    ylim(c(0,0.4)),
                g2[[3]]+ggtitle("mito % after filteration")+ 
                    ylim(c(0,0.4))))
dev.off()
######################################

# After removing unwanted cells from the dataset, the next step is to normalize the data.
BladderCancer <- NormalizeData(object = BladderCancer, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
BladderCancer <- FindVariableGenes(object = BladderCancer, mean.function = ExpMean, 
                          dispersion.function = LogVMR, do.plot = FALSE, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(BladderCancer@var.genes)
#======1.3 1st run of pca-tsne  =========================
BladderCancer <- ScaleData(object = BladderCancer) %>%
    RunPCA() %>%
    FindClusters(dims.use = 1:20, force.recalc = T, print.output = FALSE) %>%
    RunTSNE()

p1 <- TSNEPlot(object = BladderCancer, do.label = F, group.by = "orig.ident", 
         do.return = TRUE, no.legend = F, #colors.use = singler.colors,
         pt.size = 1,label.size = 8 )+
    ggtitle("Before batch correction")+
    theme(text = element_text(size=15),							
          plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

save(BladderCancer, file = "./data/BladderCancer_H2_20181109.Rda")



#======1.4 align seurat objects =========================
#Now we can run a single integrated analysis on all cells!
#Now we align the CCA subspaces, which returns a new dimensional reduction called cca.aligned
p3 <- MetageneBicorPlot(BladderCancer, grouping.var = "orig.ident", dims.eval = 1:30, 
                        display.progress = FALSE)

jpeg(paste0(path,"/CCA_MetageneBicorPlot.jpeg"), units="in", width=10, height=7,res=600)
p3
dev.off()

set.seed(42)
BladderCancer <- AlignSubspace(object = BladderCancer, reduction.type = "cca", grouping.var = "orig.ident", 
                     dims.align = 1:20)

BladderCancer <- RunTSNE(object = BladderCancer, reduction.use = "cca.aligned", dims.use = 1:20, 
               do.fast = TRUE)
BladderCancer <- FindClusters(object = BladderCancer, reduction.type = "cca.aligned", dims.use = 1:20, 
                    resolution = 0.6, force.recalc = T, save.SNN = TRUE)

BladderCancer <- RunPCA(object = BladderCancer, pcs.compute = 30, do.print = TRUE, 
              pcs.print = 1:5, genes.print = 5)

p2 <- TSNEPlot(object = BladderCancer, do.label = F, group.by = "orig.ident", 
               do.return = TRUE, no.legend = T, #colors.use = singler.colors,
               pt.size = 1,label.size = 8 )+
    ggtitle("After batch correction")+
    theme(text = element_text(size=15),							
          plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 


jpeg(paste0(path,"/TSNEplot_CCA_alignment.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p1 +theme(legend.position="none"), p2)
dev.off()

jpeg(paste0(path,"/TSNEplot.jpeg"), units="in", width=10, height=7,res=600)
p3
dev.off()

table(BladderCancer@meta.data$orig.ident)


save(BladderCancer, file = "./data/BladderCancer_H2_20181109.Rda")

#======1.4 Add Cell-cycle score =========================
# Read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "R/seurat_resources/regev_lab_cell_cycle_genes.txt")
# We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- FilterGenes(BladderCancer,cc.genes[1:43])
g2m.genes <- FilterGenes(BladderCancer,cc.genes[44:97])
# Assign Cell-Cycle Scores
BladderCancer <- CellCycleScoring(object = BladderCancer, s.genes = s.genes, g2m.genes = g2m.genes, 
                                  set.ident = TRUE)
# Visualize the distribution of cell cycle markers across
RidgePlot(object = BladderCancer, features.plot = FilterGenes(BladderCancer,c("CCND1","CDK4","CCND2","CDK6","CCND3","RB1")), 
          nCol = 2)
# regressing out the difference between the G2M and S phase scores
BladderCancer@meta.data$CC.Difference <- BladderCancer@meta.data$S.Score - BladderCancer@meta.data$G2M.Score
# view cell cycle scores and phase assignments
head(x = BladderCancer@meta.data)

#======1.5 Performing MNN-based correction =========================
#https://bioconductor.org/packages/3.8/workflows/vignettes/simpleSingleCell/inst/doc/work-5-mnn.html#4_performing_mnn-based_correction
set.seed(100)
original <- lapply(samples, function(x) BladderCancer@scale.data[BladderCancer@var.genes, 
                                                                 (BladderCancer@meta.data$orig.ident %in% x)])
mnn.out <- do.call(fastMNN, c(original, list(k=20, d=50, auto.order=T,
                                             approximate=TRUE)))
dim(mnn.out$corrected)
rownames(mnn.out$corrected) = BladderCancer@cell.names
colnames(mnn.out$corrected) = paste0("MNN_",1:ncol(mnn.out$corrected))
#Storing a new MNN
BladderCancer <- SetDimReduction(object = BladderCancer, reduction.type = "MNN", slot = "cell.embeddings",
                                 new.data = mnn.out$corrected)
BladderCancer <- SetDimReduction(object = BladderCancer, reduction.type = "MNN", slot = "key", 
                                 new.data = "MNN_")
remove(original);GC()
BladderCancer <- SetAllIdent(BladderCancer,id = "orig.ident")
DimPlot(object = BladderCancer, reduction.use = "MNN", pt.size = 0.5)

#======1.7 unsupervised clustering based on MNN =========================
BladderCancer <- RunTSNE(object = BladderCancer, reduction.use = "MNN", dims.use = 1:50, 
                do.fast = TRUE, perplexity= 30)

BladderCancer <- FindClusters(object = BladderCancer, reduction.type = "MNN", 
                    dims.use = 1:50, resolution = 0.6, 
                     k.param = 30,force.recalc = T,
                     save.SNN = TRUE, n.start = 100, nn.eps = 0, print.output = FALSE)

p4 <- TSNEPlot.1(object = BladderCancer, do.label = F, group.by = "orig.ident", 
           do.return = TRUE, no.legend = T, 
           pt.size = 1,label.size = 4 )+
    ggtitle("After batch correction")+
    theme(text = element_text(size=15),							
          plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"/TSNEplot_MNN_alignment.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p1 +theme(legend.position="none"), p4)
dev.off()


p3 <- TSNEPlot.1(object = BladderCancer, do.label = T, group.by = "ident", 
                 do.return = TRUE, no.legend = T, 
                 pt.size = 1,label.size = 6 )+
    ggtitle("Tsne plot for all clusters")+
    theme(text = element_text(size=15),							
          plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"tsneplot.jpeg"), units="in", width=10, height=7,res=600)
p3
dev.off()


p4 <- SplitTSNEPlot(object = BladderCancer, split.by = "orig.ident",
                    do.label = T, group.by = "ident", 
                 do.return = TRUE, no.legend = T, do.print =T,
                 pt.size = 1,label.size = 6 )

jpeg(paste0(path,"split_tsneplot.jpeg"), units="in", width=10, height=7,res=600)
print(do.call(plot_grid,p4))
dev.off()

save(BladderCancer, file = "./data/BladderCancer_H2_20181109.Rda")