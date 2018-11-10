########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(kableExtra)
library(SingleR)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 3.1.1 load data
# Rename ident
lnames = load(file="./output/BladderCancer_20181026.RData")

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

# Make color scheme vector
# names must be ordered vector or factor
as.data.frame(table(vertical_bar))
color = c("#53B400","#00C094","#00B6EB","#A58AFF","#FB61D7","#F8766D","#C49A00")
MakeCorlorBar <- function(df, color =NULL, remove.legend =F){
        g <- ggplot(data = df, aes(x, major_cells, fill = time_points)) +
                geom_tile()+
                theme_bw() + 
                theme(panel.border = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.line = element_blank(), 
                      axis.title.x=element_blank(),
                      axis.text.x=element_blank(),
                      axis.ticks.x=element_blank(),
                      axis.title.y=element_blank(),
                      axis.text.y=element_blank(),
                      axis.ticks.y=element_blank())
        if(remove.legend) g = g + theme(legend.position="none")
        if(!is.null(color)) {
                colors_fill = color_scheme(color = color,
                                           names = unique(df$time_points))
                g = g + scale_fill_manual(values = colors_fill)
        }
        g
}

show_col(hue_pal()(7))
table(b_color_bar[b_color_bar$major_cells == "Spermatids","time_points"])
g_legend <- MakeCorlorBar(df = b_color_bar, cell_type = "Spermatogonia",
                          remove.legend = F, color = color)
g_Spermatogonia <- MakeCorlorBar(df = b_color_bar, cell_type = "Spermatogonia",
                                 color = color)

# gene set heatmap  ===========
GeneSets <- c("Luminal_markers","EMT_and_smooth_muscle","EMT_and_claudin_markers",
              "Basal_markers","Squamous_markers")
BladderCancer <- SetAllIdent(BladderCancer,id = "orig.ident")
jpeg(paste0(path,"/Heatmap.jpeg"), units="in", width=10, height=7,res=600)
DoHeatmap(BladderCancer,data.use = scale(t(BladderCancer@meta.data[,GeneSets]),center = F),
          col.low = "#07e007",col.mid = "#FFFFFF", col.high = "#e00707",
          group.label.rot = T,cex.row = 8, group.cex = 15,slim.col.label = TRUE, 
          remove.key =F,title="Muscle-invasive bladder cancer lineage scores in mouse samples")

dev.off()

# histgram  ===========
data.use <- BladderCancer@meta.data[,GeneSets] %>% t() %>% #scale(center = F) %>%
        t() %>% as.data.frame() %>% gather(key = Subtypes.markers, value = ave.expr)
jpeg(paste0(path,"/density.jpeg"), units="in", width=10, height=7,res=600)
ggplot(data.use, aes(x = ave.expr, fill = Subtypes.markers)) +
        geom_density(alpha = .5) + scale_y_sqrt() +
        ggtitle("Muscle-invasive bladder cancer lineage scores in mouse samples")+
        theme(text = element_text(size=15),							
              plot.title = element_text(hjust = 0.5,size = 15, face = "bold")) 
dev.off()

