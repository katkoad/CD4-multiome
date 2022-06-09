library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

# path to filtered seurat object
path <- "/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Data/Donor4/15/Seurat/Donor4_15_subset.rds"

# Load in the pre-made annotations from miraldi drive
annotations <- readRDS('/data/miraldiNB/wayman/databank/genes/hg38/annot_hg38_EnsDbHsapiens_v86_UCSC.rds')
fragpath <- "/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Data/Donor4/15/Seurat/atac_fragments.tsv.gz"


# Load in seurat object, call peaks, and filter peaks
seurat.obj <- readRDS(path)
Annotation(seurat.obj) <- annotations
peaks <- CallPeaks(seurat.obj, macs2.path='/usr/local/python/2.7.15/bin/macs2')
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# FeatureMatrix
macs2_counts <- FeatureMatrix(
  fragments = Fragments(seurat.obj),
  features = peaks,
  cells = colnames(seurat.obj)
)

# Create an assay for the seurat object containing the peaks
seurat.obj[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = Annotation(seurat.obj)
)

saveRDS(seurat.obj, path)
