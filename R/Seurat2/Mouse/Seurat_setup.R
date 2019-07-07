########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(scran)
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
samples <- df_samples$samples
projects <- df_samples$projects
conditions <- df_samples$conditions
df_samples %>% kable() %>% kable_styling()

BladderCancer_raw <- list()
BladderCancer_Seurat <- list()
for(i in 1:length(samples)){
    BladderCancer_raw[[i]] <- Read10X(data.dir = paste0("./data/",
                                             samples[i],"/outs/filtered_gene_bc_matrices/mm10/"))
    colnames(BladderCancer_raw[[i]]) <- paste0(samples[i],
                                            "_",colnames(BladderCancer_raw[[i]]))
    BladderCancer_Seurat[[i]] <- CreateSeuratObject(BladderCancer_raw[[i]],
                                                 min.cells = 3,
                                                 min.genes = 200,
                                                 project = projects[i],
                                                 names.delim = "_")
    BladderCancer_Seurat[[i]]@meta.data$conditions <- conditions[i]
}
#---QC----
QC_list <- lapply(BladderCancer_Seurat, function(x) as.matrix(x = x@raw.data))
lapply(QC_list, function(x) median(colSums(x))) # Median nUMI
lapply(QC_list, function(x) median(apply(x,2,function(y) sum(length(y[y>0]))))) # Median nGene

lapply(QC_list, function(x) min(colSums(x))) # Median nUMI
lapply(QC_list, function(x) min(apply(x,2,function(y) sum(length(y[y>0]))))) # Median nGene

BladderCancer <- Reduce(function(x, y) MergeSeurat(x, y, do.normalize = F), BladderCancer_Seurat)
remove(BladderCancer_raw,BladderCancer_Seurat);GC()
BladderCancer <- FilterCells(BladderCancer, subset.names = "nGene",
                    low.thresholds = 200,
                    high.thresholds = Inf) %>%
    NormalizeData() %>%
    ScaleData(display.progress = FALSE) %>%
    FindVariableGenes(do.plot = FALSE, display.progress = FALSE)
save(BladderCancer, file = "./data/BladderCancer_M2_20181026.Rda")

#======1.2 QC, pre-processing and normalizing the data=========================
# 1.2.1 Calculate median UMI per cell
Iname = load(file = "./data/BladderCancer_M2_20181026.Rda")
# 1.2.3 calculate mitochondria percentage
mito.genes <- grep(pattern = "^mt-", x = rownames(x = BladderCancer@data), value = TRUE)
percent.mito <- Matrix::colSums(BladderCancer@raw.data[mito.genes, ])/Matrix::colSums(BladderCancer@raw.data)
BladderCancer <- AddMetaData(object = BladderCancer, metadata = percent.mito, col.name = "percent.mito")

BladderCancer@ident = factor(BladderCancer@ident,levels = samples)

g1 <- VlnPlot(object = BladderCancer, features.plot = c("nGene", "nUMI", "percent.mito"), 
              nCol = 1,point.size.use = 0.2,
              x.lab.rot = T, do.return = T,return.plotlist =T)

BladderCancer <- FilterCells(object = BladderCancer, subset.names = c("nGene","nUMI","percent.mito"),
                    low.thresholds = c(400,1500, -Inf), 
                    high.thresholds = c(Inf,Inf, 0.1))

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
    ggtitle("Original")+
    theme(text = element_text(size=15),							
          plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

save(BladderCancer, file = "./data/BladderCancer_M2_20181026.Rda")
Iname = load("./data/BladderCancer_M2_20181026.Rda")

#======1.4 Add Cell-cycle score =========================
# Read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "R/seurat_resources/regev_lab_cell_cycle_genes.txt")
# We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- MouseGenes(BladderCancer,cc.genes[1:43])
g2m.genes <- MouseGenes(BladderCancer,cc.genes[44:97])
# Assign Cell-Cycle Scores
BladderCancer <- CellCycleScoring(object = BladderCancer, s.genes = s.genes, g2m.genes = g2m.genes, 
                         set.ident = TRUE)
# Visualize the distribution of cell cycle markers across
RidgePlot(object = BladderCancer, features.plot = MouseGenes(BladderCancer,c("CCND1","CDK4","CCND2","CDK6","CCND3","RB1")), 
          nCol = 2)
# regressing out the difference between the G2M and S phase scores
BladderCancer@meta.data$CC.Difference <- BladderCancer@meta.data$S.Score - BladderCancer@meta.data$G2M.Score
# view cell cycle scores and phase assignments
head(x = BladderCancer@meta.data)

#======1.5 Add project id =========================
batchname = BladderCancer@meta.data$orig.ident
batch.effect = as.numeric(factor(batchname,levels = samples))
names(batch.effect) = rownames(BladderCancer@meta.data)
BladderCancer <- AddMetaData(object = BladderCancer, metadata = batch.effect, col.name = "batch.effect")
table(BladderCancer@meta.data$batch.effect)
head(x = BladderCancer@meta.data)

#======1.6 vars.to.regress ScaleData =========================
features_threshold <- data.frame(c("nUMI","nGene","batch.effect","percent.mito","CC.Difference"),
                        c(10000,2000,1.0,0.05,0.05))
for(i in 1:nrow(features_threshold)){
    jpeg(paste0(path,"S2_",features_threshold[i,1],".jpeg"), units="in", width=10, height=7,res=600)
    P <- SingleFeaturePlot.1(BladderCancer, feature = features_threshold[i,1],
                        threshold= features_threshold[i,2])
    print(P)
    dev.off()
}
BladderCancer <- ScaleData(object = BladderCancer, 
                  model.use = "linear", do.par=T, do.center = T, do.scale = T,
                  #vars.to.regress = c("nUMI","percent.mito","batch.effect","CC.Difference"),
                  display.progress = T)
#======1.7 Performing MNN-based correction =========================
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
BladderCancer <- RunPCA(object = BladderCancer, pc.genes = BladderCancer@var.genes, pcs.compute = 100, 
               do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PCAPlot(object = BladderCancer)
PCElbowPlot(object = BladderCancer, num.pc = 100)
PCHeatmap(BladderCancer, pc.use = c(1:3, 25:30), cells.use = 500, do.balanced = TRUE)

DimElbowPlot.1(object = BladderCancer, reduction.type = "MNN", 
             dims.plot = 50,slot = "cell.embeddings")

BladderCancer <- RunTSNE(object = BladderCancer, reduction.use = "MNN", dims.use = 1:50, 
                do.fast = TRUE, perplexity= 30)

BladderCancer <- FindClusters(object = BladderCancer, reduction.type = "MNN", 
                    dims.use = 1:50, resolution = 0.6, 
                     k.param = 30,force.recalc = T,
                     save.SNN = TRUE, n.start = 100, nn.eps = 0, print.output = FALSE)

p2 <- TSNEPlot.1(object = BladderCancer, do.label = F, group.by = "orig.ident", 
           do.return = TRUE, no.legend = T, 
           pt.size = 1,label.size = 4 )+
    ggtitle("Corrected")+
    theme(text = element_text(size=15),							
          plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"remove_batch.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p1 +theme(legend.position="none"),p2)
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
                 do.return = TRUE, no.legend = T, 
                 pt.size = 1,label.size = 6 )

jpeg(paste0(path,"split_tsneplot.jpeg"), units="in", width=10, height=7,res=600)
print(do.call(plot_grid,p4))
dev.off()

save(BladderCancer, file = "./data/BladderCancer_M2_20181026.Rda")
