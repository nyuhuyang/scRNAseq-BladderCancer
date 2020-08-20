########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(reshape2)
library(gplots)
library(MAST)
library(cowplot)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 3.1.1 load data
# Rename ident
(load(file = "data/BladderCancer_mm10_6_20190726.Rda"))
# markers by clusters =========
Idents(object) <-"integrated_snn_res.0.6"
object <- sortIdent(object,numeric = T)
TSNEPlot(object)
DefaultAssay(object) <- "RNA"
BladderCancer.markers <- FindAllMarkers.UMI(object = object, only.pos = F, logfc.threshold = 1,
                                        test.use = "MAST")
write.csv(BladderCancer.markers,paste0(path,"BladderCancer_Mouse_markers_clusters_logfc0.25.csv"))

top <-  BladderCancer.markers %>% group_by(cluster) %>% top_n(Top_n, avg_logFC)
add.genes = unique(as.character(top$gene))
object %<>% ScaleData(features = add.genes)
DoHeatmap.1(object, features = add.genes, Top_n = Top_n, do.print=T, angle = 0,
            group.bar = T, title.size = 20, no.legend = F,size=5,hjust = 0.5,
            label=F, cex.row=6, legend.size = NULL,width=10, height=7,
            pal_gsea = FALSE, unique.name = T,
            title = paste("Top",Top_n,"markers in each cluster"))

DE_files <- c("BladderCancer_Mouse_markers_clusters_logfc0.25.csv",
              "BladderCancer_Human_357-Bladder_logfc0.25.csv",
              "BladderCancer_Human_359-Bladder_logfc0.25.csv")
Top_n = 10
for(file in DE_files){
        BladderCancer.markers = read.csv(paste0(path, file), row.names = 1)
        top <-  BladderCancer.markers %>% group_by(cluster) %>% top_n(Top_n, avg_logFC)
        bottom <-  BladderCancer.markers %>% group_by(cluster) %>% top_n(-Top_n, avg_logFC)
        top_bottom_10 <- rbind(top,bottom)
        top_bottom_10 = top_bottom_10[order(top_bottom_10$cluster),]
        write.csv(top_bottom_10,paste0(path,file))
        
}

# Gene heatmap arrange by gene set ======================
GeneSets = read.delim("data/seurat_resources/gene_set_Bladder_Cancer.txt",row.names =1,header = F,
                      stringsAsFactors = F)
GeneSets %>% kable() %>% kable_styling()
GeneSets.df <- as.data.frame(t(GeneSets))
GeneSetNames <- c("Luminal","EMT_and_smooth_muscle","EMT_and_claudin",
                         "Basal","Squamous","Immune",
                         "Down_regulated_in_CIS","Up_regulated_in_CIS")
GeneSets.df <- GeneSets.df[,GeneSetNames]

GeneSets.list <- df2list(GeneSets.df)
#GeneSets_list <- lapply(GeneSets.list, function(x) Human2Mouse(x))
GeneSets_list <- lapply(GeneSets.list, function(x) FilterGenes(object,x))
heatmap_genes <- unlist(GeneSets_list,use.names = F)

object %<>% ScaleData(features = heatmap_genes)
Idents(object) %<>% factor(levels =c(3,8,5,11,7,0,2,10,20,17,12,15,14,9,19,21,18,16,4,1,6,13))
DoHeatmap.1(object, features = heatmap_genes, do.print=T, angle = 0,
            group.bar = F, title.size = 20, no.legend = F,size=5,hjust = 0.5,
            group.bar.height = 0,label=F, cex.row= 500/length(heatmap_genes), legend.size = 0,
            width=8, height=6.5, unique.name = F, pal_gsea = F,
            title = "mRNA Expression Subtypes (Mouse)")


marker_df <- gather(list2df(GeneSets_list))
colnames(marker_df) = c("cluster","gene")
marker_df = marker_df[nchar(marker_df$gene) >0 & !is.na(nchar(marker_df$gene)),]
marker_df = marker_df[!duplicated(marker_df$gene),]
marker_df$cluster[marker_df$cluster %in% c("Down_regulated_in_CIS", 
                                           "Up_regulated_in_CIS")] = "Down_up_regulated_in_CIS"
marker_df$cluster %<>% as.factor
marker_df$cluster %<>% factor(levels = unique(marker_df$cluster))
GeneSetColors <- c("#386CB0","#ff0000","#6A3D9A","#E6AB02","#F0027F","#BF5B17","#4DAF4A") 
MakeCorlorBar(object, marker_df, color = GeneSetColors)


Idents(object) = "groups"
for (sample in c("4950", "8524", "8525")) {
        subset_object <- subset(object, idents = sample)
        Idents(subset_object) <-"integrated_snn_res.0.6"
        subset_object %<>% ScaleData(features = heatmap_genes)
        Idents(subset_object) %<>% factor(levels =c(3,8,5,11,7,0,2,10,20,17,12,15,14,9,19,21,18,16,4,1,6,13))
        DoHeatmap.1(subset_object, features = heatmap_genes, do.print=T, angle = 0,
                    group.bar = T, title.size = 20, no.legend = T,size=5,hjust = 0.5,
                    group.bar.height = 0.02,label=T, cex.row= 500/length(heatmap_genes), 
                    legend.size = 0,
                    width=8, height=6.5, unique.name = T, pal_gsea = F,
                    title = paste("mRNA Expression Subtypes in",sample,"(Mouse)"))
        }


# gene set heatmap  ===========
Idents(object) <-"integrated_snn_res.0.6"
integrated_snn_res.0.6 <- AverageExpression(object,assays = "RNA", slot = "scale.data",
                                            features = heatmap_genes)

y = as.matrix(integrated_snn_res.0.6$RNA)
## Column clustering (adjust here distance/linkage methods to what you need!)
hc <- hclust(as.dist(1-cor(y, method="pearson")))

jpeg(paste0(path,"Mouse_Heatmap_geneSet_zscore.jpeg"), units="in", width=10, height=7,res=600)
heatmap.2(y,
          Colv = as.dendrogram(hc), 
          Rowv= FALSE,
          trace ="none",
          dendrogram = "column",
          key.xlab = "Row z-score")
dev.off()

# geom_density  ===========
BladderCancer_subset <- SplitObject(object,split.by = "groups")
(samples <- sapply(BladderCancer_subset,function(x) unique(x@meta.data$groups)))
BladderCancer_subset[[4]] = object
(samples = c(samples,"all"))
GeneSetNames <- c("Basal","EMT_and_claudin",
                  "EMT_and_smooth_muscle","Luminal","Squamous")

g <- list()
for(i in 1:length(samples)){
        data.use <- BladderCancer_subset[[i]]@meta.data[,GeneSetNames] %>% t() %>% #scale(center = F) %>%
                t() %>% as.data.frame() %>% gather(key = Subtypes.markers, value = ave.expr)
        g[[i]] <- ggplot(data.use, aes(x = ave.expr, fill = Subtypes.markers)) +
                geom_density(alpha = .5) + scale_y_sqrt() +
                theme(legend.position="none")+
                xlab("Average expression (log nUMI)")+
                ggtitle(paste("density plot for",samples[i]))+
                theme_bw() +
                theme(text = element_text(size=15),
                      legend.position=c(0.5,0.8),
                      plot.title = element_text(hjust = 0.5,size = 20, face = "plain"),
                      panel.border = element_blank(), 
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), 
                      axis.line = element_line(colour = "black"))
        jpeg(paste0(path,"DensityPlot_",samples[i],".jpeg"), units="in", width=10, height=7,res=600)
        print(g[[i]])
        dev.off()
        }

# histgram  ===========
BladderCancer_subset <- SplitSeurat(BladderCancer)
samples <- BladderCancer_subset[[length(BladderCancer_subset)]]

g1 <- list()
for(i in 1:length(samples)){
        data.use <- BladderCancer_subset[[i]]@meta.data[,GeneSets] %>% t() %>% #scale(center = F) %>%
                t() %>% as.data.frame() %>% gather(key = Subtypes.markers, value = ave.expr)
        g1[[i]] <- ggplot(data.use, aes(x = ave.expr, fill = Subtypes.markers)) +
                geom_histogram(binwidth=0.03, alpha=0.5, position="identity") + scale_y_sqrt() +
                theme(legend.position="none")+
                ggtitle(samples[i])+
                theme(text = element_text(size=15),
                      legend.position=c(0.35,0.8),
                      plot.title = element_text(hjust = 0.5,size = 15, face = "bold"))
}
jpeg(paste0(path,"Mouse_histogram.jpeg"), units="in", width=10, height=7,res=600)
do.call(plot_grid,g1)+
        ggtitle("Muscle-invasive bladder cancer lineage scores in mouse samples")+
        theme(text = element_text(size=15),							
              plot.title = element_text(hjust = 0.5,size = 15, face = "bold"))
dev.off()
