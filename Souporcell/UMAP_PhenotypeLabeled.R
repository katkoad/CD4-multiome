library(Seurat)

Seurat.obj <- readRDS("/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Data/Donor3/5/Seurat/Seurat.obj.rds")
barcodes <- read.table("/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/Souporcell_out/Donor3/15_commonvariants/cluster1_barcodes.txt")
barcodes <- as.matrix(barcodes)

seurat.obj <- NormalizeData(seurat.obj)
seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst")
seurat.obj <- ScaleData(seurat.obj)
seurat.obj <- RunPCA(seurat.obj, features = VariableFeatures(object = seurat.obj))
seurat.obj <- RunUMAP(seurat.obj, dims = 1:20)


tiff("Umap_Donor3Integrated.tiff", units="in", width=5, height=5, res=300)
DimPlot(seurat.obj, reduction = "umap", cells.highlight = barcodes, sizes.highlight = 0.1)
dev.off()