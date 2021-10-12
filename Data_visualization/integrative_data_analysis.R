library(ArchR)
library(Seurat)
library(grid)

source('SpatialDimPlot_new.R')

### Integrate with scRNA-seq reference data 
## https://oncoscape.v3.sttrcancer.org/atlas.gs.washington.edu.mouse.rna/downloads
## Prepare MOCA data
MOCA_dir <- "./ref_data/MOCA/"

meta.data.RNA <- read.csv(file = paste0(MOCA_dir, 'cell_annotate.csv'), header = TRUE, row.names = 1, stringsAsFactors = FALSE)
gene.ANN.RNA <- read.csv(file = paste0(MOCA_dir, 'gene_annotate.csv'), header = TRUE, row.names = 1, stringsAsFactors = FALSE)
gene.ANN.RNA <- gene.ANN.RNA[, 'gene_short_name', drop = FALSE]

cds <- readRDS(paste0(MOCA_dir, 'gene_count_cleaned_sampled_100k.RDS'))

MOCA <- CreateSeuratObject(counts = cds, project = 'MOCA')
meta.data.RNA <- meta.data.RNA[colnames(MOCA), ]
meta.data.RNA <- meta.data.RNA[, c('Main_cell_type', 'development_stage')]

MOCA <- AddMetaData(object = MOCA, metadata = meta.data.RNA)
MOCA_E11 <- subset(MOCA, development_stage == 11.5)
MOCA_E11.raw.data <- as.matrix(GetAssayData(MOCA_E11, slot = 'counts'))
MOCA_E11.raw.data <- as.data.frame(MOCA_E11.raw.data)
MOCA_E11.raw.data <- merge(gene.ANN.RNA, MOCA_E11.raw.data, by=0, all=TRUE)
which(is.na(MOCA_E11.raw.data$gene_short_name))

tt <- table(MOCA_E11.raw.data$gene_short_name)
name_rep <- names(which(tt > 1))
row_del_fun <- function(x){
  rows <- which(MOCA_E11.raw.data$gene_short_name == x)
  return(rows[2:length(rows)] )
}
row_del <- unlist(lapply(name_rep, row_del_fun))
MOCA_E11.raw.data <- MOCA_E11.raw.data[-row_del, ]

row.names(MOCA_E11.raw.data) <- MOCA_E11.raw.data$gene_short_name
MOCA_E11.raw.data <- MOCA_E11.raw.data[, -c(1:2), drop=FALSE]
MOCA_E11 <- CreateSeuratObject(counts = MOCA_E11.raw.data, project = 'MOCA_E11', meta.data = MOCA_E11@meta.data)


## Integration with ArchR oject
proj_in_tissue <- addGeneIntegrationMatrix(
  ArchRProj = proj_in_tissue, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = MOCA_E11,
  addToArrow = TRUE,
  groupRNA = "Main_cell_type",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  force = TRUE
)


## Plot results
meta.data.integration <- as.data.frame(getCellColData(ArchRProj = proj_in_tissue))[, c('predictedCell', 'predictedGroup', 'predictedScore')]
new_row_names <- row.names(meta.data.integration)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub(".*#","", x)))
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(meta.data.integration) <- new_row_names

spatial.obj <- AddMetaData(object = spatial.obj, metadata = meta.data.integration)

Idents(spatial.obj) <- 'predictedGroup'

ids.highlight <- names(table(spatial.obj$predictedGroup))[1]
ids.highlight

p <- SpatialDimPlot_new(spatial.obj, cells.highlight = CellsByIdentities(object = spatial.obj, idents = ids.highlight), 
                        facet.highlight = TRUE, pt.size.factor = 2.5, alpha = c(1,0), stroke = 0)
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22)
p



### Integrate with scCUT&Tag reference data (H3K4me3)
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163532
K4.spatial      <- readRDS(file='/data/proj/GCB_MB/spatial_cut-tag/results/CT/H3K4me3/clustering/Seurat_object.Rds')
K4.spatial      <- RenameAssays(object = K4.spatial,bin_5000='bins_5000')

p1 <- DimPlot(K4.spatial)
p2 <- DimPlot(K4.spatial,group.by = 'YD.clusters')
p1+p2

K4.single_cell <- readRDS(file='/data/proj/GCB_MB/single-cell-CUT-Tag/nbiotech_paper/analysis/results/H3K4me3/clustering/01.clustering.Rds')

assay = 'bins_5000'

DefaultAssay(K4.spatial)      <- assay
DefaultAssay(K4.single_cell)  <- assay

min_reads = 5
peaks.common <- table(c(rownames(K4.spatial),rownames(K4.single_cell))) == 2
peaks.common <- peaks.common[peaks.common]

features.common.table <- table(c(rownames(K4.spatial)[rowSums(K4.spatial[[assay]]@counts) > min_reads],
                                 rownames(K4.single_cell)[rowSums(K4.single_cell[[assay]]@counts) > min_reads]))

peaks.use <- names(features.common.table[features.common.table == 2])

anchors <- FindIntegrationAnchors(
  object.list = list(K4.single_cell,K4.spatial),
  anchor.features = peaks.use,
  assay = rep(assay,2),
  k.filter = NA,reference = 2
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

p1 <- DimPlot(integrated[,integrated$orig.ident != 'spatial'],pt.size=0.1,label=TRUE)
p2 <- DimPlot(integrated[,integrated$orig.ident == 'spatial'],pt.size=0.1,label=TRUE,group.by='YD.clusters')
p1+p2



### Integrate with scCUT&Tag reference data (H3K27me3)
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163532
bin=5000
assay = paste0('bins_',bin)

K27.spatial <- readRDS(file='/data/proj/GCB_MB/spatial_cut-tag/results/CT/H3K27me3/clustering/H3K27me3_clustering.Rds')

p1 <- DimPlot(K27.spatial)
p2 <- DimPlot(K27.spatial,group.by = 'YD.clusters')
p1+p2

K27.single_cell <- readRDS(file='/data/proj/GCB_MB/single-cell-CUT-Tag/nbiotech_paper/analysis/results/H3K27me3/clustering/01.clustering.Rds')

DefaultAssay(K27.spatial)      <- assay
DefaultAssay(K27.single_cell)  <- assay

min_reads = 5
features.common.table <- table(c(rownames(K27.spatial)[rowSums(K27.spatial[[assay]]@counts) > min_reads],
                                 rownames(K27.single_cell)[rowSums(K27.single_cell[[assay]]@counts) > min_reads]))

peaks.use <- names(features.common.table[features.common.table == 2])

anchors <- FindIntegrationAnchors(
  object.list = list(K27.single_cell,K27.spatial),
  anchor.features = peaks.use,
  assay = rep(assay,2),
  k.filter = NA,reference = 1
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

p1 <- DimPlot(integrated[,integrated$orig.ident != 'spatial'],pt.size=0.1,label=TRUE)
p2 <- DimPlot(integrated[,integrated$orig.ident == 'spatial'],pt.size=0.1,label=TRUE,group.by='YD.clusters')
p1+p2

SpatialDimPlot(K27.spatial,group.by = 'YD.clusters',pt.size.factor = 6)





