args = commandArgs(trailingOnly=TRUE)
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(dplyr))
suppressMessages(library(HGNChelper))

# Path to seurat object to cluster
path <- args[1]
seurat.obj <- readRDS(paste(norm, "integrated.rds", sep=""))
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

RNA_Dimensions <- args[2]
ATAC_Dimensions <- args[3]
integrate <- args[4]
norm <- args[5]
res <- args[6]


## Cluster with the Harmony integration assays
if(integrate == "Harmony"){
    seurat.obj <- FindMultiModalNeighbors(
        object = seurat.obj,
        reduction.list = list("harmony_rna", "harmony_atac"), 
        dims.list = list(RNA_Dimensions, ATAC_Dimensions),
        modality.weight.name = "RNA.weight",
        verbose = TRUE
    )
}

## Cluster with the Default integration assays
if(integrate == "Default")
seurat.obj <- FindMultiModalNeighbors(
  object = seurat.obj,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(RNA_Dimensions, ATAC_Dimensions),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
seurat.obj <- RunUMAP(seurat.obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
seurat.obj <- FindClusters(seurat.obj, graph.name = "wsnn", algorithm = 3, resolution = res, verbose = FALSE)

# Plot UMAP
tiff(paste(norm, "_", integrate, "_", "wnn_umap.tiff"), units="in", width=5, height=5, res=300)
DimPlot(seurat.obj, label = FALSE, reduction = "wnn.umap", group.by = "dataset")
dev.off()

markers <- FindAllMarkers(seurat.obj)
write.csv(markers, "cluster_markers.tsv", quote=F, sep="\t")

gs_list = gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", "Immune system")
cell_types = sctype_score(scRNAseqData = seurat.obj, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
View(cell_types)
