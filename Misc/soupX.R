library(SoupX)
library(Seurat)
library(ggplot2)

## Load Seurat Object and create meta data

seurat.obj <- readRDS("/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Data/Donor4/5/Seurat/Donor4_5.rds")
seurat.obj <- NormalizeData(seurat.obj)
seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst")
seurat.obj <- ScaleData(seurat.obj)
seurat.obj <- RunPCA(seurat.obj, features = VariableFeatures(object = seurat.obj))
seurat.obj <- FindNeighbors(seurat.obj, dims=1:30)
seurat.obj <- FindClusters(seurat.obj, resolution = 0.3)
seurat.obj <- RunUMAP(seurat.obj, dims = 1:30)

## Check number of clusters and change res if needed
head(Idents(seurat.obj))

tiff("Umap2.tiff", units="in", width=5, height=5, res=300)
DimPlot(seurat.obj, reduction = "umap")
dev.off()

## Create metadata
metaData = as.data.frame(seurat.obj$umap@cell.embeddings)
colnames(metaData) <- c("RD1", "RD2")
metaData$Cluster <- seurat.obj@meta.data$seurat_clusters

## Since dont have annotations yet, make it same as metadata
metaData$Annotation <- metaData$Cluster

## Create SoupX object
dir = "/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/soupX/Donor4/5"
sc = load10X(dir)
sc = setClusters(sc, setNames(metaData$Cluster, rownames(metaData)))
sc = setDR(sc, metaData[colnames(sc$toc), c("RD1", "RD2")])

tiff("SoupX_Umap2.tiff", units="in", width=5, height=5, res=300)
dd = metaData[colnames(sc$toc), ]
mids = aggregate(cbind(RD1, RD2) ~ Annotation, data = dd, FUN = mean)
gg = ggplot(dd, aes(RD1, RD2)) + geom_point(aes(colour = Annotation), size = 0.2) + 
    geom_label(data = mids, aes(label = Annotation)) + ggtitle("PBMC 4k Annotation") + 
    guides(colour = guide_legend(override.aes = list(size = 1)))
plot(gg)
dev.off()

tiff("SoupX_contam2.tiff", units="in", width=5, height=5, res=300)
sc = autoEstCont(sc)
dev.off()
