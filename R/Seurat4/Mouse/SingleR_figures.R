# conda activate r4.1.3
library(Seurat)
library(magrittr)
library(pheatmap)
library(dplyr)
library(tidyr)
library(ggpubr)
library(S4Vectors)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 load data ==========================================
(load(file = "data/BladderCancer_mm10_6_20190726.Rda"))
meta.data <- object@meta.data
##############################
# create singleR data frame
###############################
pred <- readRDS(file = "output/BladderCancer_mm10_6_20190726_Tabula_2022_singleR_pred.rds")

singlerDF = data.frame("ontology_id_uberon_id" = pred$pruned.labels,
                   row.names = rownames(pred))
table(rownames(pred) == rownames(object@meta.data))
table(is.na(singlerDF$ontology_id_uberon_id))
singlerDF$ontology_id_uberon_id[is.na(singlerDF$ontology_id_uberon_id)]= "unknown"
singlerDF$ontology_id <- gsub("-UBERON.*","",singlerDF$ontology_id_uberon_id)
singlerDF$uberon_id <- gsub(".*-UBERON","UBERON",singlerDF$ontology_id_uberon_id)


FilePath = "../../../Rancho BioSciences/Projects/scds-pipeline"
Tabula_2022 = readRDS(file.path(FilePath,"output/Tabula_2022_Science-FigShare100973-pseudo-bulk_seurat.RDS"))
meta_data = Tabula_2022@meta.data
meta_data1 = meta_data[!duplicated(meta_data$ontology_id),]
meta_data2 = meta_data[!duplicated(meta_data$uberon_id),]

singlerDF$cell_type <- plyr::mapvalues(singlerDF$ontology_id,
                                      from = meta_data1$ontology_id,
                                      to =  meta_data1$cell_type)
singlerDF$cell_type_compartment <- plyr::mapvalues(singlerDF$ontology_id,
                                                  from = meta_data1$ontology_id,
                                                  to =  meta_data1$cell_type_compartment)
singlerDF$sample_tissue_uberon_name <- plyr::mapvalues(singlerDF$uberon_id,
                                                      from = meta_data2$uberon_id,
                                                      to =  meta_data2$sample_tissue_uberon_name)

if(all(colnames(object) == rownames(singlerDF))){
        print("all cellID match!")
        object@meta.data %<>% cbind(singlerDF)
}

saveRDS(object@meta.data, file = "output/BladderCancer_mm10_6_20190726_meta.data.rds")
