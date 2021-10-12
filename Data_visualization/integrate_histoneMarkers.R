rm(list = ls())

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(patchwork)
set.seed(1234)

## Read annotation
ensdb = EnsDb.Mmusculus.v79
gene.coords <- ensembldb::genes(ensdb, filter = ~ gene_biotype == "protein_coding")
lncRNA.coords <- ensembldb::genes(ensdb, filter = ~ gene_biotype == "lincRNA")
gene.coords <- c(gene.coords,lncRNA.coords)
seqlevelsStyle(gene.coords) <- 'UCSC'
genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')

# Flatten the overlapping genes and extend by 2kb upstream of promoters
genebody.coords.flat <- GenomicRanges::reduce(x = genebody.coords)
genebodyandpromoter.coords.flat <- Signac::Extend(genebody.coords.flat,upstream = 2000)

# Retrieve gene names from the original annotation (lost because of flatenning)
genebodyandpromoter.coords.flat$name<- gene.coords[nearest(genebodyandpromoter.coords.flat,genebody.coords)]$gene_name


## H3K4me3
# Create peak counts matrix
peaks_file = paste0('./H3K4me3_peaks.narrowPeak')
peaks <- rtracklayer::import(peaks_file)

metadata_K4me3 <- read.table(
  file = "./H3K4me3_metadata.txt",
  header = TRUE,
  row.names = 1
)
row.names(metadata_K4me3) <- paste0(row.names(metadata_K4me3), '-1')

fragments_K4m3 <- CreateFragmentObject('./fragments.tsv.gz', cells = row.names(metadata_K4me3))

counts.matrix.peaks <- FeatureMatrix(fragments = fragments_K4m3,
                                     features = peaks,
                                     cells = row.names(metadata_K4me3),
                                     process_n = 20)

counts.matrix.bins <- GenomeBinMatrix(fragments = fragments_K4m3, 
                                      genome = seqlengths(BSgenome.Mmusculus.UCSC.mm10),
                                      binsize = 5000,
                                      cells = row.names(metadata_K4me3), 
                                      process_n = 5000)

K4me3 <- CreateSeuratObject(
  counts = counts.matrix.bins,
  assay = 'bins',
  project = 'H3K4me3',
  meta.data = metadata_K4me3
)

K4me3[['peaks']] <- CreateAssayObject(counts = counts.matrix.peaks[,colnames(counts.matrix.peaks) %in% colnames(K4me3)])
K4me3

# Add sample id to cell names
K4me3 <- RenameCells(object = K4me3, add.cell.id = 'H3K4me3')

# Create Gene activity matrix
gene.matrix <- FeatureMatrix(fragments = fragments_K4m3,
                             features = genebodyandpromoter.coords.flat,
                             cells = gsub(paste0('H3K4me3',"_"),"",colnames(K4me3)))

genes.key             <- genebodyandpromoter.coords.flat$name
names(genes.key)      <- GRangesToString(genebodyandpromoter.coords.flat)
rownames(gene.matrix) <- genes.key[rownames(gene.matrix)]

gene.matrix           <- gene.matrix[rownames(gene.matrix) != "",]
colnames(gene.matrix) <- paste0('H3K4me3',"_",colnames(gene.matrix))

K4me3[['GA']] <- CreateAssayObject(counts = gene.matrix)

DefaultAssay(K4me3) <- 'bins'
K4me3 <- RunTFIDF(K4me3)
K4me3 <- FindTopFeatures(K4me3, min.cutoff = 'q20')
K4me3 <- RunSVD(object = K4me3)
DepthCor(K4me3)

K4me3 <- NormalizeData(
  object = K4me3,
  assay = 'GA',
  normalization.method = 'LogNormalize',
  scale.factor = median(K4me3$nCount_GA)
)


## H3K27ac
peaks_file = paste0('./H3K27ac_peaks.narrowPeak')
peaks <- rtracklayer::import(peaks_file)

metadata_K4me3 <- read.table(
  file = "./H3K27ac_metadata.txt",
  header = TRUE,
  row.names = 1
)
row.names(metadata_K27ac) <- paste0(row.names(metadata_K27ac), '-1')

fragments_K27ac <- CreateFragmentObject('./K27ac_ME11_20um/fragments.tsv.gz', cells = row.names(metadata_K27ac))

counts.matrix.bins <- GenomeBinMatrix(fragments = fragments_K27ac, 
                                      genome = seqlengths(BSgenome.Mmusculus.UCSC.mm10),
                                      binsize = 5000,
                                      cells = row.names(metadata_K27ac), 
                                      process_n = 5000)

K27ac <- CreateSeuratObject(
  counts = counts.matrix.bins,
  assay = 'bins',
  project = 'H3K27ac',
  meta.data = metadata_K27ac
)

K27ac[['peaks']] <- CreateAssayObject(counts = counts.matrix.peaks[,colnames(counts.matrix.peaks) %in% colnames(K27ac)])

# Add sample id to cell names
K27ac <- RenameCells(object = K27ac,add.cell.id = 'H3K27ac')

# Create Gene activity matrix
gene.matrix <- FeatureMatrix(fragments = fragments_K27ac,
                             features = genebodyandpromoter.coords.flat,
                             cells = gsub(paste0('H3K27ac',"_"),"",colnames(K27ac)))

genes.key             <- genebodyandpromoter.coords.flat$name
names(genes.key)      <- GRangesToString(genebodyandpromoter.coords.flat)
rownames(gene.matrix) <- genes.key[rownames(gene.matrix)]

gene.matrix           <- gene.matrix[rownames(gene.matrix) != "",]
colnames(gene.matrix) <- paste0('H3K27ac',"_",colnames(gene.matrix))

K27ac[['GA']] <- CreateAssayObject(counts = gene.matrix)

DefaultAssay(K27ac) <- 'bins'
K27ac <- RunTFIDF(K27ac)
K27ac <- FindTopFeatures(K27ac, min.cutoff = 'q20')
K27ac <- RunSVD(object = K27ac)

DepthCor(K27ac)

K27ac <- NormalizeData(
  object = K27ac,
  assay = 'GA',
  normalization.method = 'LogNormalize',
  scale.factor = median(K27ac$nCount_GA)
)
##


## Integrate data across different histone marks
assay = 'GA'

DefaultAssay(K4me3)  <- assay
DefaultAssay(K27ac)  <- assay

min_reads = 5

features.common.table <- table(c(rownames(K4me3)[Matrix::rowSums(K4me3[[assay]]@counts) > min_reads],
                                 rownames(K27ac)[Matrix::rowSums(K27ac[[assay]]@counts) > min_reads]))
features.common.table

peaks.use <- names(features.common.table[features.common.table == 2])

anchors <- FindIntegrationAnchors(
  object.list = list(K4me3, K27ac),
  anchor.features = peaks.use,
  assay = rep(assay,2),
  k.filter = NA,
  reference = 1
)

integrated <- IntegrateData(
  anchorset = anchors,
  preserve.order = TRUE
)

integrated <- RunSVD(
  object = integrated,
  n = 50,
  reduction.name = 'integratedLSI'
)

integrated <- RunUMAP(
  object = integrated,
  dims = 2:40,
  reduction = 'integratedLSI')

integrated <- FindNeighbors(
  object = integrated,
  reduction = 'integratedLSI',
  dims = 2:40
)

integrated <- FindClusters(
  object = integrated,
  resolution = 1,
  verbose = FALSE
)

cols <- c("#90D5E4", "#208A42", "#C06CAB", "#D51F26", "#89288F", "#FEE500")
p <- DimPlot(integrated, split.by = 'orig.ident', cols = cols)
p