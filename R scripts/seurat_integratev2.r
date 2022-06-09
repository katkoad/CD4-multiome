args = commandArgs(trailingOnly=TRUE)
library(Seurat)
library(sctransform)
library(harmony)


# Arguments
path <- args[1]
norm <- args[2]
integration <- args[3]

# Load in seurat objects of the four timepoints
Resting <- readRDS(paste(path,"Resting_subset.rds",sep=""))
Two <- readRDS(paste(path, "2_subset.rds", sep=""))
Five  <- readRDS(paste(path,"5_subset.rds", sep=""))
Fifteen <- readRDS(paste(path, "15_subset.rds", sep=""))

# Change default assay to RNA in case it's currently ATAC
DefaultAssay(Resting) <- "RNA"
DefaultAssay(Two) <- "RNA"
DefaultAssay(Five) <- "RNA"
DefaultAssay(Fifteen) <- "RNA"

Resting$dataset = "Resting"
Two$dataset = "Two"
Five$dataset = "Five"
Fifteen$dataset = "Fifteen"

## Harmony RNA Integration
if(integration == "Harmony") {
    immune.combined <- readRDS("merged.rds")
    DefaultAssay(immune.combined) <- "RNA"
    # Default Normalization 
    if(norm == "Default"){
        immune.combined <- NormalizeData(immune.combined)
        immune.combined <- FindVariableFeatures(immune.combined)
        immune.combined <- ScaleData(immune.combined)
        immune.combined <- RunPCA(immune.combined, dims=1:30)
	immune.combined <- RunUMAP(immune.combined, dims=1:30)
	tiff(paste(norm, "_", integration, "_rna_merged_umap.tiff", sep=""),units="in", width=5, height=5, res=300)
	print(DimPlot(immune.combined, group.by="dataset", reduction = "umap"))
	dev.off()
	immune.combined <- RunHarmony(immune.combined, "dataset", reduction.save = "harmony_rna", assay.use="RNA")
    }
    # SCT Normalization
    if(norm == "SCT"){
        immune.combined <- SCTransform(immune.combined)
        immune.combined <- RunPCA(immune.combined, dims=1:30)
	immune.combined <- RunUMAP(immune.combined, dims=1:30)
	tiff(paste(norm, "_", integration, "_rna_merged_umap.tiff", sep=""),units="in", width=5, height=5, res=300)
        print(DimPlot(immune.combined, group.by="dataset", reduction = "umap"))
        dev.off()
	immune.combined <- RunHarmony(immune.combined, "dataset", reduction.save = "harmony_rna", assay.use="SCT")
    }
    immune.combined <- RunUMAP(immune.combined, dims=1:30, reduction = "harmony_rna")
    tiff(paste(norm, "_", integration, "_rna_integrated_umap.tiff", sep=""),units="in", width=5, height=5, res=300)
    print(DimPlot(immune.combined, group.by="dataset", reduction = "umap"))
    dev.off()
}

## Seurat RNA Integration
if(integration == "Default") {
    # Default Normalization
    if(norm == "Default"){
        Resting <- NormalizeData(Resting)
        Resting <- FindVariableFeatures(Resting, selection.method = "vst", nfeatures = 2000)
        Two <- NormalizeData(Two)
        Two <- FindVariableFeatures(Two, selection.method = "vst", nfeatures = 2000)
        Five <- NormalizeData(Five)
        Five <- FindVariableFeatures(Five, selection.method = "vst", nfeatures = 2000)
        Fifteen <- NormalizeData(Fifteen)
        Fifteen <- FindVariableFeatures(Fifteen, selection.method = "vst", nfeatures = 2000)
    }
    # SCT Normalization
    if(norm == "SCT"){
        Resting <- SCTransform(Resting)
        Two <- SCTransform(Two)
        Five <- SCTransform(Five)
        Fifteen <- SCTransform(Fifteen)
    }
    # SelectIntegrationFeatures
    features <- SelectIntegrationFeatures(c(Resting, Two, Five, Fifteen))
    # FindIntegrationAnchors
    immune.anchors <- FindIntegrationAnchors(c(Resting, Two, Five, Fifteen), anchor.features = features)
    # Integrate
    immune.combined <- IntegrateData(anchorset = immune.anchors)
    immune.combined <- ScaleData(immune.combined)
    immune.combined <- RunPCA(immune.combined, dims=1:30)
    immune.combined <- RunUMAP(immune.combined, dims=1:30)
    tiff(paste(norm, "_", integration, "_rna_merged_umap.tiff", sep=""),units="in", width=5, height=5, res=300)
    print(DimPlot(immune.combined, group.by="dataset", reduction = "umap"))
    dev.off()
}


DefaultAssay(Resting) <- "peaks"
DefaultAssay(Two) <- "peaks"
DefaultAssay(Five) <- "peaks"
DefaultAssay(Fifteen) <- "peaks"
# Get the grange object of the peaks
Resting_peaks <- granges(Resting)
Two_peaks <- granges(Two)
Five_peaks <- granges(Five)
Fifteen_peaks <- granges(Fifteen)
# Find overlapping peaks
combined.peaks <- reduce(c(Resting_peaks, Two_peaks, Five_peaks, Fifteen_peaks))
# Create counts matrix using set of overlapping peaks for each object
Resting.counts <- FeatureMatrix(
  fragments = Fragments(Resting),
  features = combined.peaks,
  cells = colnames(Resting)
)
Two.counts <- FeatureMatrix(
  fragments = Fragments(Two),
  features = combined.peaks,
  cells = colnames(Two)
)
Five.counts <- FeatureMatrix(
  fragments = Fragments(Five),
  features = combined.peaks,
  cells = colnames(Five)
)
Fifteen.counts <- FeatureMatrix(
  fragments = Fragments(Fifteen),
  features = combined.peaks,
  cells = colnames(Fifteen)
)
# Create assay and then an object for each set
Resting_peaks_integrated_assay <- CreateAssayObject(counts = Resting.counts, min.cells = 1)
Two_peaks_integrated_assay <- CreateAssayObject(counts = Two.counts, min.cells = 1)
Five_peaks_integrated_assay <- CreateAssayObject(counts = Five.counts, min.cells = 1)
Fifteen_peaks_integrated_assay <- CreateAssayObject(counts = Fifteen.counts, min.cells = 1)
Resting_peaks_integrated <- CreateSeuratObject(Resting_peaks_integrated_assay, assay = "ATAC")
Two_peaks_integrated <- CreateSeuratObject(Two_peaks_integrated_assay, assay = "ATAC")
Five_peaks_integrated <- CreateSeuratObject(Five_peaks_integrated_assay, assay = "ATAC")
Fifteen_peaks_integrated <- CreateSeuratObject(Fifteen_peaks_integrated_assay, assay = "ATAC")

# Add in dataset meta data for each of the new objects
Resting_peaks_integrated$dataset <- "Resting"
Two_peaks_integrated$dataset <- "Two"
Five_peaks_integrated$dataset <- "Five"
Fifteen_peaks_integrated$dataset <- "Fifteen"

# merge the objects
atac_merged <- merge(x=Resting_peaks_integrated, y=c(Two_peaks_integrated, Five_peaks_integrated, Fifteen_peaks_integrated))
# Normalizae the merged object and plot unintegrated umap
atac_merged <- RunTFIDF(atac_merged)
atac_merged <- FindTopFeatures(atac_merged, min.cutoff = 20)
atac_merged <- RunSVD(atac_merged)
atac_merged <- RunUMAP(atac_merged, dims = 2:50, reduction = 'lsi')
tiff(paste(norm, "_", integration, "_atac_merged_umap.tiff", sep=""),units="in", width=5, height=5, res=300)
print(DimPlot(atac_merged, group.by="dataset", reduction = "umap"))
dev.off()


if(integration == "Harmony"){
    atac_merged <- ScaleData(atac_merged)
    atac_merged <- RunHarmony(atac_merged, "dataset", reduction = "lsi", reduction.save = "harmony_atac", assay.use="ATAC")
    atac_merged <- RunUMAP(atac_merged, dims = 2:50, reduction = 'harmony_atac')
    tiff(paste(norm, "_", integration, "_atac_integrated_umap.tiff", sep=""),units="in", width=5, height=5, res=300)
    print(DimPlot(atac_merged, group.by="dataset", reduction = "umap"))
    dev.off()
    integrated.data <- atac_merged[["harmony_atac"]]
    immune.combined[["integrated_atac"]] <- integrated.data
}

## Default ATAC Integration
if(integration == "Default") {
    # Integrate
    integration.anchors <- FindIntegrationAnchors(
        object.list = list(Resting_peaks_integrated, Two_peaks_integrated, Five_peaks_integrated, Fifteen_peaks_integrated),
        anchor.features = rownames(Resting),
     )
    integrated <- IntegrateEmbeddings(
        anchorset = integration.anchors,
        reductions = atac_merged[["lsi"]],
        new.reduction.name = "integrated_lsi",
        dims.to.integrate = 1:30
    )
    integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)

    # Plot integrated umap
    tiff(paste(norm, "_", integration, "_atac_integrated_umap.tiff", sep=""),units="in", width=5, height=5, res=300)
    print(DimPlot(atac_merged, group.by="dataset", reduction = "umap"))
    dev.off()

    # Add the integrated atac assay to our main seurat object and save it
    integrated.data <- GetAssayData(object =  integrated[['ATAC']], slot = 'data')
    immune.combined[["integrated_atac"]] <- CreateAssayObject(data = integrated.data)
}
saveRDS(immune.combined, paste(norm, "_", integration, "_integrated.rds", sep=""))

