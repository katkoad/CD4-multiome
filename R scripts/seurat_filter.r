args = commandArgs(trailingOnly=TRUE)
library(Seurat)
library(dplyr)
library(ggplot2)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

path <- args[1]
nGenes_min <- args[2]
nUMI_min <- args[3]
mitoratio_max <- args[4]
blacklist_max <- args[5]
Nucleosomesignal_max <- args[6]
TSSenrichment_min <- args[7]
frag.path <- args[8]
save_name <- args[9]

Seurat.obj <- readRDS(path)

Seurat.obj$log10GenesPerUMI <- log10(Seurat.obj$nFeature_RNA) / log10(Seurat.obj$nCount_RNA)

Seurat.obj$mitoRatio <- PercentageFeatureSet(object = Seurat.obj, pattern = "^MT-")
Seurat.obj$mitoRatio <- Seurat.obj@meta.data$mitoRatio / 100

metadata <- Seurat.obj@meta.data
metadata$cells <- rownames(metadata)
metadata <- metadata %>%
        dplyr::rename(seq_folder = orig.ident,
                      nUMI = nCount_RNA,
                      nGene = nFeature_RNA)

Seurat.obj@meta.data <- metadata

tiff("qc_umi.tiff", units="in", width=5, height=5, res=300)
metadata %>% 
  	ggplot(aes(x=nUMI, color="none", fill="none")) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") +
  	geom_vline(xintercept = 60) + theme(legend.position="none")
dev.off()

tiff("qc_genes.tiff", units="in", width=5, height=5, res=300)
metadata %>% 
  	ggplot(aes(color="none", x=nGene, fill= "none")) + 
  	geom_density(alpha = 0.2) + 
  	theme_classic() +
  	scale_x_log10() + 
  	geom_vline(xintercept = 300) + theme(legend.position="none")
dev.off()

tiff("qc_umi_vs_genes.tiff", units="in", width=5, height=5, res=300)
metadata %>% 
  	ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
  	stat_smooth(method=lm) +
  	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 600) +
  	geom_hline(yintercept = 300) 
dev.off()

tiff("mito_genes.tiff", units="in", width=5, height=5, res=300)
metadata %>% 
  	ggplot(aes(color="none", x=mitoRatio, fill="none")) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 0.2) + theme(legend.position="none")
dev.off()

tiff("complexity.tiff", units="in", width=5, height=5, res=300)
metadata %>%
  	ggplot(aes(x=log10GenesPerUMI, color = "none", fill="none")) +
  	geom_density(alpha = 0.2) +
  	theme_classic() +
  	geom_vline(xintercept = 0.8) + theme(legend.position="none")
dev.off()



DefaultAssay(Seurat.obj) <- "ATAC"
fragments <- CreateFragmentObject(
  path = frag.path,
  cells = colnames(Seurat.obj), 
  validate.fragments = TRUE
)
Fragments(Seurat.obj) <- NULL
Fragments(Seurat.obj) <- fragments

Seurat.obj <- NucleosomeSignal(object = Seurat.obj)
Seurat.obj <- TSSEnrichment(object = Seurat.obj, fast = FALSE)

Seurat.obj$blacklist_fraction <- FractionCountsInRegion(
  object = Seurat.obj, 
  assay = 'ATAC',
  regions = blacklist_hg38
)

tiff("blacklist.tiff", units="in", width=5, height=5, res=300)
Seurat.obj@meta.data %>%
  	ggplot(aes(x=blacklist_fraction, color = "none", fill="none")) +
  	geom_density(alpha = 0.2) +
  	theme_classic() +
  	geom_vline(xintercept = 0.05) + theme(legend.position="none")
dev.off()

tiff("nucleosome.tiff", units="in", width=5, height=5, res=300)
Seurat.obj@meta.data %>%
  	ggplot(aes(x=nucleosome_signal, color = "none", fill="none")) +
  	geom_density(alpha = 0.2) +
  	theme_classic() +
  	geom_vline(xintercept = 2.5) + theme(legend.position="none")
dev.off()

tiff("tss.tiff", units="in", width=5, height=5, res=300)
Seurat.obj@meta.data %>%
  	ggplot(aes(x=TSS.enrichment, color = "none", fill="none")) +
  	geom_density(alpha = 0.2) +
  	theme_classic() +
  	geom_vline(xintercept = 2) + theme(legend.position="none")
dev.off()

annotations <- readRDS('/data/miraldiNB/wayman/databank/genes/hg38/annot_hg38_EnsDbHsapiens_v86_UCSC.rds')
Annotation(Seurat.obj) <- annotations
peaks <- CallPeaks(Seurat.obj, macs2.path='/usr/local/python/2.7.15/bin/macs2')
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# FeatureMatrix
macs2_counts <- FeatureMatrix(
  fragments = Fragments(Seurat.obj),
  features = peaks,
  cells = colnames(Seurat.obj)
)

# Create an assay for the seurat object containing the peaks
Seurat.obj[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = Fragments(Seurat.obj),
  annotation = Annotation(Seurat.obj),
  min.features=-1
)


Seurat.obj_subset <- subset(Seurat.obj, subset = nGene > 400 & nUMI > 750 & mitoRatio < 20 & blacklist_fraction < 0.05 & nucleosome_signal < 2.5 & TSS.enrichment > 2 )

saveRDS(Seurat.obj_subset, save_name)

