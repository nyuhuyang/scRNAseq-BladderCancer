# conda activate r4.1.1
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
library(data.table)
library(RColorBrewer)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_differential_expression.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

object <- readRDS("data/Sfakianos_2020_Nat_Commun-GSE146137_mouse-seurat-annotated.RDS")
object %<>% AddMetaColor(label = "Celltypist_Cell_Type_Ontology_Name",
             colors = colorRampPalette(brewer.pal(n =12,"Paired"))(length(unique(object$Celltypist_Cell_Type_Ontology_ID_mouse))))
#meta.data <- readRDS(file = "output/BladderCancer_mm10_6_20190726_meta.data.rds")
#if(all(colnames(object) == rownames(meta.data))){
#        print("all cellID match!")
#        object@meta.data = meta.data
#}

#============== re-run UMAP ==============================
# Building trajectories with Monocle 3
object@reductions$umap = object@reductions$cca.umap
colnames(object@reductions$umap@cell.embeddings) = paste0("UMAP_",1:2)
object@reductions$umap@key = "UMAP_"
object@meta.data %<>% cbind(object[["umap"]]@cell.embeddings)

epi_object <- subset(object, subset = UMAP_1 > -1 & UMAP_2 > -5 )#& orig.ident %in% c("4950N","8524N","8525N"))
epi_object$Celltypist_Cell_Type_Ontology_Name %<>% droplevels()
UMAPPlot.1(epi_object,group.by = "Celltypist_Cell_Type_Ontology_Name",do.print = TRUE)

cds <- as.cell_data_set(epi_object)
cds <- cluster_cells(cds = cds, reduction_method = "UMAP",k = 13,verbose = TRUE)
cds <- learn_graph(cds, use_partition = TRUE, close_loop = TRUE)


get_earliest_principal_node <- function(cds, cell.type="basal cell"){
    cell_ids <- which(colData(cds)[, "cell_type"] == cell.type)
    
    closest_vertex <-
        cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    root_pr_nodes <-
        igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                  (which.max(table(closest_vertex[cell_ids,]))))]
    
    root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
#cds <- order_cells(cds)

epi_cds <- cds[,colData(cds)$Celltypist_Cell_Type_Ontology_Name %in% c("bladder urothelial cell","fibroblast",
                                                                          "myofibroblast cell","pericyte cell",
                                                                          "smooth muscle cell")]

colData(epi_cds)$Celltypist_Cell_Type_Ontology_Name %<>% droplevels()

epi_object <- subset(epi_object, subset = Celltypist_Cell_Type_Ontology_Name %in% c("bladder urothelial cell",
                                                                                    "fibroblast","myofibroblast cell",
                                                                                   "pericyte cell","smooth muscle cell") )
epi_object$Celltypist_Cell_Type_Ontology_Name %<>% droplevels()
#=============== prepare figures ===================

g <- plot_cells(epi_cds,color_cells_by = "cell.types",show_trajectory_graph = TRUE,group_label_size = 5,graph_label_size=5,cell_size = 0.5)

save_path <- paste0(path,"total/")
if(!dir.exists(save_path)) dir.create(save_path, recursive = T)

for(group.by in c("cell.types","Celltypist_Cell_Type_Ontology_Name")){
        print(group.by)
        save_sub_path <- paste0(save_path,group.by)
        if(!dir.exists(save_sub_path)) dir.create(paste0(save_sub_path), recursive = T)
        for(label in c("subtype_raw","subtype","subtype_label")){
                g <- plot_cells(
                        cds = epi_cds,
                        cell_size = 0.5,
                        label_cell_groups = switch(label,
                                                   "subtype_raw" = TRUE,
                                                   "subtype" = TRUE,
                                                   "subtype_label" = TRUE),
                        color_cells_by = group.by,
                        show_trajectory_graph = switch(label,
                                                       "subtype_raw" = FALSE,
                                                       TRUE),
                        group_label_size = switch(label,
                                                  "subtype_raw" = 5,
                                                  "subtype_label" = 5,
                                                  0),
                        graph_label_size=5) +
                        scale_color_manual(values = ExtractMetaColor(epi_object, group.by = group.by))
                jpeg(paste0(save_sub_path,"/Mouse_bladderCancer_",group.by,"_",label,".jpeg"), units="in", width=7, height=7,res=600)
                print(g)
                dev.off()
        }
}


for(label in c("pseudotime","pseudotime_label")){
    save_sub_path <- paste0(save_path,"pseudotime")
    if(!dir.exists(save_sub_path)) dir.create(paste0(save_sub_path), recursive = T)
    g1 <-     plot_cells(
                cds = epi_cds,
                cell_size = 0.5,
                color_cells_by = "pseudotime",
                show_trajectory_graph = switch(label,
                                               "pseudotime" = FALSE,
                                               "pseudotime_label" = TRUE
                                               ),
                scale_to_range = TRUE,
                group_label_size=5,
                graph_label_size=5
                )
    jpeg(paste0(save_sub_path,"/Mouse_bladderCancer_",label,".jpeg"), units="in", width=7, height=7,res=600)
    print(g1)
    dev.off()
}

### each tumor ###
samples = c("4950N","8524N","8525N")

for(sample in samples){
        print(sample)
        save_path <- paste0(path,"/",sample)
        if(!dir.exists(save_path)) dir.create(save_path, recursive = T)
        
        cds_subset <- epi_cds[,pData(epi_cds)$orig.ident == sample]
        cds_subset <- order_cells(cds_subset, root_pr_nodes=get_earliest_principal_node(cds_subset))
        #=============== prepare figures ===================

        
        
        for(group.by in c("cell.types","Celltypist_Cell_Type_Ontology_Name")){
            print(group.by)
            save_sub_path <- paste0(save_path,"/",group.by)
            if(!dir.exists(save_sub_path)) dir.create(paste0(save_sub_path), recursive = T)
            for(label in c("subtype_raw","subtype","subtype_label")){
                g <- plot_cells(
                    cds = cds_subset,
                    cell_size = 0.5,
                    label_cell_groups = switch(label,
                                               "subtype_raw" = TRUE,
                                               "subtype" = TRUE,
                                               "subtype_label" = TRUE),
                    color_cells_by = group.by,
                    show_trajectory_graph = switch(label,
                                                   "subtype_raw" = FALSE,
                                                   TRUE),
                    group_label_size = switch(label,
                                              "subtype_raw" = 5,
                                              "subtype_label" = 5,
                                              0),
                    graph_label_size=5) +
                    scale_color_manual(values = ExtractMetaColor(epi_object, group.by = group.by))
                jpeg(paste0(save_sub_path,"/",sample,"_",group.by,"_",label,".jpeg"), units="in", width=7, height=7,res=600)
                print(g)
                dev.off()
            }
        }
        
        for(label in c("pseudotime","pseudotime_label")){
            save_sub_path <- paste0(save_path,"/pseudotime")
            if(!dir.exists(save_sub_path)) dir.create(paste0(save_sub_path), recursive = T)
            g1 <-     plot_cells(
                cds = cds_subset,
                cell_size = 0.5,
                color_cells_by = "pseudotime",
                show_trajectory_graph = switch(label,
                                               "pseudotime" = FALSE,
                                               "pseudotime_label" = TRUE
                ),
                scale_to_range = TRUE,
                group_label_size=5,
                graph_label_size=5
            )
            jpeg(paste0(save_sub_path,"/",sample,"_",label,".jpeg"), units="in", width=7, height=7,res=600)
            print(g1)
            dev.off()
        }
}

epi_object
epi_object$Cluster <- epi_cds@clusters$UMAP$clusters
UMAPPlot.1(epi_object,do.print = TRUE, group.by = "Cluster")
UMAPPlot.1(epi_object,do.print = TRUE, group.by = "Cluster",label =T)

table(epi_object$orig.ident)
epi_object$orig.ident
Idents(epi_object) = "Cluster"
markers <- FindAllMarkers_UMI(epi_object,group.by ="Cluster",logfc.threshold = 0.1,
                   only.pos = F)
write.csv(markers,paste0(path,"DEGs_Clusters.csv"))
