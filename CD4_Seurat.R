library(Seurat)

rds_directory <- "/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Data"
data_location <- c("Resting/Resting.rds", "Resting_Set2/resting_set2.rds", "2/2.rds", "2_Set2/2_set2.rds", "5/5.rds" , 
                     "5_Set2/5_set2.rds", "15/15.rds", "15_Set2/15_set2.rds")
object_names <- c("Resting", "Resting_2", "Two", "Two_2", "Five", "Five_2", "Fifteen", "Fifteen_2")

obj.list <- list()

for(i in 1:length(object_names)){
  obj.list[[object_names[i]]] <- readRDS(file.path(rds_directory, data_location[i]))
  obj.list[[i]][["percent.mt"]] <- PercentageFeatureSet(obj.list[[i]], pattern = "^MT-")
}

for(i in 1:length(obj.list)){
  pdf(paste(rds_directory,object_names[i],".pdf", sep=""))
  print(VlnPlot(obj.list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  dev.off()
}




