########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(gplots)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 3.1.1 load data
# Rename ident
(load(file = "data/BladderCancer_mm10_6_20190706.Rda"))
Idents(object) <-"RNA_snn_res.0.6"
TSNEPlot(object)
BladderCancer.markers <- FindAllMarkers.UMI(object = object, only.pos = F, 
                                        test.use = "MAST")
write.csv(BladderCancer.markers,paste0(path,"BladderCancer_markers.csv"))
BladderCancer.markers =read.csv(file = paste0(path,"BladderCancer_markers.csv"),
                                row.names = 1, stringsAsFactors=F)

DoHeatmap.1(object,BladderCancer.markers,Top_n = 3, do.print=T,angle = 0,
            group.bar = T,title.size = 20, no.legend = F,size=4,label=T,
            title = "Top 3 markers in all clusters in 4950PN, PP, and B5011")

#split by samples
object$tumor = gsub('4950PN|4950PP','4950',object$orig.ident)
Idents(object) <- "RNA_snn_res.0.6"

TSNEPlot.1(object,group.by = "RNA_snn_res.0.6",split.by = "tumor",do.print=T)

Idents(object) <- "tumor"
B5011_4950.markers <- FindAllMarkers.UMI(object = object, only.pos = F, 
                                            test.use = "MAST")
write.csv(B5011_4950.markers,paste0(path,"B5011_4950.markers.csv"))
DoHeatmap.1(object,B5011_4950.markers,Top_n = 25, do.print=T,angle = 0,
            group.bar = T,title.size = 20, no.legend = F,size=4,label=T,
            title = "Top 25 markers in between 4950 and B5011")

# gene set heatmap  ===========
GeneSets <- c("Luminal_markers","EMT_and_smooth_muscle","EMT_and_claudin_markers",
              "Basal_markers","Squamous_markers")
object <- SetAllIdent(object,id = "orig.ident")

y = t(object@meta.data[,GeneSets]) %>% scale()
## Column clustering (adjust here distance/linkage methods to what you need!)
hc <- hclust(as.dist(1-cor(y, method="pearson")), method="complete")
cc = gsub("_.*","",hc$labels)
cc = gsub("4950PP","#E31A1C",cc)
cc = gsub("4950PN","#33A02C",cc)

jpeg(paste0(path,"/Mouse_Heatmap_geneSet_zscore.jpeg"), units="in", width=10, height=7,res=600)
heatmap.2(y,
          Colv = as.dendrogram(hc), Rowv= FALSE,
          ColSideColors = cc, trace ="none",labCol = FALSE,dendrogram = "column",#scale="row",
          adjRow = c(1, NA),offsetRow = -52, cexRow = 1.5,
          key.xlab = "Row z-score",
          col = bluered)
par(lend = 1)           # square line ends for the color legend
legend(0, 0.83,       # location of the legend on the heatmap plot
       legend = c("4950PP", "4950PN"), # category labels
       col = c("#E31A1C", "#33A02C"),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)
dev.off()

# geom_density  ===========
BladderCancer_subset <- SplitSeurat(BladderCancer)
samples <- BladderCancer_subset[[length(BladderCancer_subset)]]

g <- list()
for(i in 1:length(samples)){
        data.use <- BladderCancer_subset[[i]]@meta.data[,GeneSets] %>% t() %>% #scale(center = F) %>%
                t() %>% as.data.frame() %>% gather(key = Subtypes.markers, value = ave.expr)
        g[[i]] <- ggplot(data.use, aes(x = ave.expr, fill = Subtypes.markers)) +
                geom_density(alpha = .5) + scale_y_sqrt() +
                theme(legend.position="none")+
                xlab("Average expression (log nUMI)")+
                ggtitle(samples[i])+
                theme(text = element_text(size=15),
                      legend.position=c(0.35,0.8),
                      plot.title = element_text(hjust = 0.5,size = 15, face = "bold"))
}
jpeg(paste0(path,"Mouse_density.jpeg"), units="in", width=10, height=7,res=600)
do.call(plot_grid,g)+
        ggtitle("Muscle-invasive bladder cancer lineage scores in mouse samples")+
        theme(text = element_text(size=15),							
              plot.title = element_text(hjust = 0.5,size = 15, face = "bold"))
dev.off()


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
