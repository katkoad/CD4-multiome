library(Seurat)
library(Signac)

Resting <- readRDS("Resting_subset.rds")
Two <- readRDS("2_subset.rds")
Five <- readRDS("5_subset.rds")
Fifteen <- readRDS("15_subset.rds")

Resting$dataset <- "Resting"
Two$dataset <- "Two"
Five$dataset <- "Five"
Fifteen$dataset <- "Fifteen"

merged <- merge(x=Resting, y=c(Two, Five, Fifteen))

saveRDS(merged, "merged.rds")
