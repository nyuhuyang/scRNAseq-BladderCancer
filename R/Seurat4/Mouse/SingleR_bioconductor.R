#====== 3.1 Create Singler Object  ==========================================
# conda activate r4.0.3 linux
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(magrittr)

source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/SingleR_functions.R")

# ====== load single cell =============
(load(file = "data/BladderCancer_mm10_6_20190726.Rda"))

sce <- SummarizedExperiment(list(logcounts=object[["SCT"]]@data),
                                colData=DataFrame(object@meta.data))
rm(object);GC()


# ====== load reference =============
Tabula_2022 <- readRDS("../../../Rancho BioSciences/Projects/scds-pipeline/output/Tabula_2022_Science-FigShare100973-pseudo-bulk_seurat.RDS")
Tabula_2022_sce <- SummarizedExperiment(list(logcounts=Tabula_2022[["RNA"]]@data),
                             colData=DataFrame("ontology_id_uberon_id"=Tabula_2022$ontology_id_uberon_id))
rownames(Tabula_2022_sce) <- Hmisc::capitalize(tolower(rownames(Tabula_2022_sce)))
common <- Reduce(intersect, list(rownames(sce),
                                 rownames(Tabula_2022_sce)
))
length(common)

table(Tabula_2022_sce$ontology_id_uberon_id)
system.time(trained <- trainSingleR(ref = Tabula_2022_sce[common,],
                                    labels=Tabula_2022_sce$ontology_id_uberon_id))
system.time(pred <- classifySingleR(sce[common,], trained))
# elapsed 4872.846 sec
saveRDS(object = pred, file = "output/BladderCancer_mm10_6_20190726_Tabula_2022_singleR_pred.rds")
