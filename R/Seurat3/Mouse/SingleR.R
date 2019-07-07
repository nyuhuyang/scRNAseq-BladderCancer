library(SingleR)
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file = "data/object_data_mm10_6_20190706.Rda"))
attach(mouse.rnaseq)
singler = CreateSinglerObject(object_data, annot = NULL, project.name="Mouse_BladderCancer",
                              min.genes = 500,technology = "10X", species = "Mouse", citation = "",
                              ref.list = list(mouse.rnaseq),
                              normalize.gene.length = F, variable.genes = "de",
                              fine.tune = F, do.signatures = F, clusters = NULL)
save(singler,file="output/singler_F_BladderCancer_mm10_6_20190706.Rda")
  
