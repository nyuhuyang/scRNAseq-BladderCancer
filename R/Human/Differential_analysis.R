########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(SingleR)
library(gplots)
source("R/utils/Seurat_functions.R")
source("R/utils/SingleR_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 3.1.1 load data
# Rename ident
(lnames = load(file="./data/BladderCancer_H2_20181109.Rda"))

BladderCancer <- SetAllIdent(BladderCancer,id = "res.0.6")
TSNEPlot(BladderCancer)
BladderCancer.markers <- FindAllMarkers.UMI(object = BladderCancer, only.pos = TRUE, 
                                        test.use = "MAST")
write.csv2(BladderCancer.markers,paste0(path,"BladderCancer_markers.csv"))
BladderCancer.markers %>% group_by(cluster) %>% 
        top_n(2, avg_logFC) %>% kable() %>% kable_styling()

clusters <- c(5,0,15,18,2,1,4,6,9,13,14,3,7,8,10:12,16:17)
top <- group_top_mutate(df = BladderCancer.markers, clusters)
BladderCancer.markers %>% head(20) %>% kable() %>% kable_styling()

DoHeatmap.1(BladderCancer,top,Top_n = 15, 
            group.order = clusters,ident.use = "all cell types",
            group.label.rot = T,cex.row = 5,remove.key =T)

# Gene heatmap arrange by gene set ======================
gene_set = read.delim("./doc/gene_set_Bladder_Cancer.txt",row.names =1,header = F,
                      stringsAsFactors = F)
gene_set %>% kable() %>% kable_styling()
gene_set.df <- as.data.frame(t(gene_set))
gene_set.df <- gene_set.df[,c("Luminal_markers","EMT_and_smooth_muscle","EMT_and_claudin_markers",
                              "Basal_markers","Squamous_markers","Immune_markers",
                              "Neuronal_differentiation","Down_regulated_in_CIS",
                              "Up_regulated_in_CIS")]
gene_set.list <- df2list(gene_set.df)
gene_set_list <- lapply(gene_set.list, function(x) MouseGenes(BladderCancer,x))
heatmap_genes <- unlist(gene_set_list,use.names = F)
heatmap_caterlog <- gsub("[[:digit:]]+","",names(unlist(gene_set_list)))

top1 <- heatmap_genes[(heatmap_genes %in% top$gene)]
vertical_bar <- heatmap_caterlog[(heatmap_genes %in% top$gene)]

jpeg(paste0(path,"/Heatmap~.jpeg"), units="in", width=10, height=7,res=600)
DoHeatmap(BladderCancer,genes.use = top1, #Top_n = 15,
            group.order = clusters,#ident.use = "all cell types",
            group.label.rot = F,cex.row = 8, group.cex = 8,slim.col.label = TRUE, 
          remove.key =T)
dev.off()

# gene set heatmap  ===========
GeneSets <- c("Luminal_markers","EMT_and_smooth_muscle","EMT_and_claudin_markers",
              "Basal_markers","Squamous_markers")
BladderCancer <- SetAllIdent(BladderCancer,id = "orig.ident")

y = t(BladderCancer@meta.data[,GeneSets]) #%>% scale()
## Column clustering (adjust here distance/linkage methods to what you need!)
hc <- hclust(as.dist(1-cor(y, method="pearson")), method="complete")
cc = gsub("_.*","",hc$labels)
cc = gsub("357-Bladder","#E31A1C",cc)
cc = gsub("359-Bladder","#33A02C",cc)

jpeg(paste0(path,"/Human_Heatmap_geneSet.jpeg"), units="in", width=10, height=7,res=600)
heatmap.2(y,
          Colv = as.dendrogram(hc), Rowv= FALSE,
          ColSideColors = cc, trace ="none",labCol = FALSE,dendrogram = "column",#scale="row",
          adjRow = c(1, NA),offsetRow = -52, cexRow = 1.5,
          key.xlab = "Average expression (log nUMI)",
          col = bluered)
par(lend = 1)           # square line ends for the color legend
legend(0, 0.83,       # location of the legend on the heatmap plot
       legend = c("357-Bladder", "359-Bladder"), # category labels
       col = c("#E31A1C", "#33A02C"),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)
dev.off()


# geom_density  ===========
GeneSets <- c("Luminal_markers","EMT_and_smooth_muscle","EMT_and_claudin_markers",
              "Basal_markers","Squamous_markers")
BladderCancer_subset <- SplitSeurat(BladderCancer)
(samples <- BladderCancer@meta.data$orig.ident %>% unique %>% sort)

g <-  BladderCancer_subset[[1]]@meta.data[,GeneSets] %>% as.data.frame %>% 
      gather(key = Subtypes.markers, value = ave.expr) %>%
      ggplot(aes(x = ave.expr, color = Subtypes.markers)) + 
          geom_density(size = 1) +
          scale_y_sqrt() + #ylim(0, 1)+
          xlab("lineage score")+
          #ggtitle(samples[1])+
          theme(text = element_text(size=15),
                #legend.position="none", 
                legend.position=c(0.35,0.85) ,
                plot.title = element_text(hjust = 0.5,size = 15, face = "bold"))+
        geom_vline(xintercept=0,size = 1.5) #Adding vertical line in plot ggplot

jpeg(paste0(path,"density_mouse.jpeg"), units="in", width=10, height=7,res=600)
print(g+ggtitle("Muscle-invasive bladder cancer lineage scores in CD45 negative mouse samples")+
        theme(text = element_text(size=15),							
              plot.title = element_text(hjust = 0.5,size = 15, face = "bold")))
dev.off()

# histgram  ===========
BladderCancer_subset <- SplitSeurat(BladderCancer)
(samples <- BladderCancer@meta.data$orig.ident %>% unique %>% sort)

g1 <- list()
for(i in 1:length(samples)){
      data.use <- BladderCancer_subset[[i]]@meta.data[,GeneSets] %>% t() %>% #scale(center = F) %>%
        t() %>% as.data.frame() %>% gather(key = Subtypes.markers, value = ave.expr)
      g1[[i]] <- ggplot(data.use, aes(x = ave.expr, fill = Subtypes.markers)) +
        geom_histogram(binwidth=0.03, alpha=0.5, position="identity") + scale_y_sqrt() +
        theme(legend.position="none")+
        ggtitle(samples[i])+
        theme(text = element_text(size=15),
              legend.position=c(0.4,0.8),
              plot.title = element_text(hjust = 0.5,size = 15, face = "bold"))
}
jpeg(paste0(path,"Human_histogram.jpeg"), units="in", width=10, height=7,res=600)
do.call(plot_grid,g1)+
  ggtitle("Muscle-invasive bladder cancer lineage scores in Human samples")+
  theme(text = element_text(size=15),							
        plot.title = element_text(hjust = 0.5,size = 15, face = "bold"))
dev.off()

# bar chart  ===========
#scheme one. Any cells with lineage score >0 will be assigned the correspondence subypte labels.
g <-  BladderCancer_subset[[1]]@meta.data[,GeneSets] %>%
             apply(2,function(y) length(y[y>0])) %>%
             c("Total"= nrow(BladderCancer_subset[[1]]@meta.data), .) %>%
      data.frame("Subtype" = gsub("_markers", "",names(.)),
              "cell numbers" = ., row.names = NULL)%>%
      ggplot(aes(x = reorder(Subtype, -cell.numbers), y = cell.numbers, fill = Subtype)) + 
          geom_bar(stat="identity")+theme_minimal()+
          theme(axis.text.x=element_text(angle=45, hjust=1))+
      xlab("Sub cell types")+
      scale_fill_manual(values=c(gg_color_hue(5),"#999999"))+
      ggtitle("Total and subtype bladder cancer cell numbers in CD45 negative mouse samples")+
      theme(text = element_text(size=15),							
            plot.title = element_text(hjust = 0.2,size = 15, face = "bold"))

jpeg(paste0(path,"Mouse_bar_chart.jpeg"), units="in", width=10, height=7,res=600)
g
dev.off()
