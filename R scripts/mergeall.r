library(Seurat)
library(Signac)

Donor1 <- readRDS("Donor1/merged.rds")
Donor2 <- readRDS("Donor2/merged.rds")
Donor3 <- readRDS("Donor3/merged.rds")
Donor4 <- readRDS("Donor4/merged.rds")

x <- Donor1
y <- c(Donor2, Donor3, Donor4)

merged <- merge(x, y)
saveRDS(merged, "merge_all.rds")
