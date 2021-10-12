rm(list=ls())
library(ArchR)
library(Seurat)
library(grid)

threads = 8
addArchRThreads(threads = threads)

addArchRGenome("mm10")

inputFiles <- './fragments.tsv.gz'
sampleNames <- 'ME11'

## Create ArchRProject
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = sampleNames,
  filterTSS = 0,
  filterFrags = 0,
  minFrags = 0,
  maxFrags = 1e+07,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  offsetPlus = 0,
  offsetMinus = 0,
  TileMatParams = list(tileSize = 5000)
)
ArrowFiles

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = sampleNames,
  copyArrows = TRUE
)
proj


## Select pixels in tissue
meta.data <- as.data.frame(getCellColData(ArchRProj = proj))
meta.data['cellID_archr'] <- row.names(meta.data)
data.dir <- getwd()
assay = "Spatial"
filter.matrix = TRUE
slice = "slice1"
image <- Read10X_Image(image.dir = file.path(data.dir, "spatial"), 
                       filter.matrix = filter.matrix)
meta.data.spatial <- meta.data[row.names(image@coordinates), ]
proj_in_tissue <- proj[meta.data.spatial$cellID_archr, ]
proj_in_tissue


## Data normalization and dimensionality reduction 
proj_in_tissue <- addIterativeLSI(
  ArchRProj = proj_in_tissue,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list(
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,
  force = TRUE
)

proj_in_tissue <- addClusters(
  input = proj_in_tissue,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8,
  force = TRUE
)

proj_in_tissue <- addUMAP(
  ArchRProj = proj_in_tissue, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)

plotEmbedding(ArchRProj = proj_in_tissue, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", size = 1.5)

proj_in_tissue <- addImputeWeights(proj_in_tissue)


## Identify the marker genes for each cluster 
markersGS <- getMarkerFeatures(
  ArchRProj = proj_in_tissue, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  testMethod = "wilcoxon"
)

markerList_pos <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1")
markerList_neg <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC <= -1")

# Export marker genes
for (i in seq_len(length(markerList_pos))) {
  write.table(markerList_pos[[i]], file=paste0('./markers_list/', sampleNames, '_C', i, '_markers.txt'),
              quote=FALSE, sep='\t', col.names = TRUE, row.names = TRUE)
}

markerGenes <- list()
for (i in seq_len(length(markerList_pos))) {
  markerGenes <- c(markerGenes, markerList_pos[[i]]$name)
}
markerGenes <- unlist(markerGenes)


## ChromVAR Deviatons Enrichment
proj_in_tissue <- addBgdPeaks(proj_in_tissue, force = TRUE)

proj_in_tissue <- addDeviationsMatrix(
  ArchRProj = proj_in_tissue, 
  peakAnnotation = "Motif",
  force = TRUE
)

markersMotifs <- getMarkerFeatures(
  ArchRProj = proj_in_tissue, 
  useMatrix = "MotifMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useSeqnames = 'z'
)


## bulk sample (ENCODE) projection
library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)

register(MulticoreParam(8))
addArchRGenome("mm10")
getGenome()

seq_info <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)
excludeChr <- grep("chrM|chrY", seq_info@seqnames)
seq_info_sub <- seq_info[seq_info@seqnames[-excludeChr]]

tiles <- tileGenome(seq_info_sub, tilewidth=5000, cut.last.tile.in.chrom=TRUE)

bamfile_dir <- "/ref_data/bulk_chipseq_ENCODE/H3K27me3_E11.5/"
bamfile <- list.files(bamfile_dir, pattern = 'bam')
bamfile_names <- tools::file_path_sans_ext(bamfile)
meta_data <- read.table(paste0(bamfile_dir, 'metadata.tsv'), sep = '\t', header = TRUE, fill = TRUE, stringsAsFactors = FALSE)

meta_data_sub <- meta_data[which(meta_data$File.accession %in% bamfile_names), ]
row.names(meta_data_sub) <- meta_data_sub$File.accession
meta_data_sub <- meta_data_sub[bamfile_names,]
bamfile_dir <- paste0(bamfile_dir, bamfile)

fragment_counts <- getCounts(bamfile_dir, 
                             tiles, 
                             paired = FALSE,
                             by_rg = FALSE, 
                             format = "bam", 
                             colData = DataFrame(celltype = meta_data_sub$Biosample.term.name))

projBulk <- projectBulkATAC(proj_in_tissue, seATAC = fragment_counts, n = 250)

encode_cell_types <- function(x) {
  x <- tools::file_path_sans_ext(x)
  meta_data_sub[x,]$Biosample.term.name
}
projBulk[[1]][,3] <- sapply(projBulk[[1]]$Type, encode_cell_types)

plotProj <- rbind(projBulk[[2]], projBulk[[1]])
pal <- paletteDiscrete(unique(as.vector(plotProj[,3])))
pal["spatial_CUT_Tag"] <- "lightgrey"

p <- ggPoint(plotProj[,1], plotProj[,2], as.vector(plotProj[,3]), rastr = TRUE, pal = pal)
p


## Pseudo-time analysis
trajectory <- c('Radial glia', "Postmitotic premature neurons", 'Excitatory neurons')

proj_in_tissue <- addTrajectory(
  ArchRProj = proj_in_tissue, 
  name = "Neuron_U", 
  groupBy = "cell_type",
  trajectory = trajectory, 
  embedding = "UMAP", 
  force = TRUE
)