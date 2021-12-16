
#' ---
#' title: "Cluster analysis of Neonatal Spinal Cord"
#' author: "James Choi"
#' date: "`r Sys.Date()`"
#' output: pdf_document
#' ---

#+ setup, include=TRUE
knitr::opts_chunk$set(warning=FALSE, message=FALSE, fig.align='center', 
                      tidy.opts=list(width.cutoff=60), tidy=TRUE)

# Load libraries
require('Seurat')
require('ggplot2')
require('dplyr')

# Directory setup
# setwd('scripts')
results_out <- '../results/NeonatalClusterAnalysis/'
dir.create(path = results_out)
data_path <- '../data/GSE150871_RAW/'
sample_names <- list.dirs(path = data_path, recursive = FALSE, full.names = FALSE)



# Cell barcode calling ----------------------------------------------------

#' Deposited data at GSE150871 contains the full, unfiltered count matrix i.e. 
#' all of the counts for every droplet sequenced via 10X (~7x10^6). First, need
#' to determine which to call cells via EmptyDrops() function of the 
#' DropletUtils R package (by Aaron Lun). 

data_path <- '../data/GSE150871_RAW/'
dir.create(path = results_out)

# Run EmptyDroplets algorithm
drops <- vector(mode = 'list', length = length(sample_names))
names(drops) <- sample_names
for (id in sample_names) {
  # Load data
  message(paste('Loading data:', id))
  counts <- Read10X(data.dir = paste0(data_path, id))
  
  # Calculate barcode ranks
  message('Ranking...')
  ranks <- DropletUtils::barcodeRanks(
    m = counts,
    lower = 200,
    fit.bounds = c(200, 1e6)
  )
  knee <- ranks@metadata$knee
  inflection <- ranks@metadata$inflection
  
  # Empty vs non-empty droplet
  message('Running EmptyDrops...')
  drops[[id]] <- DropletUtils::emptyDrops(
    m = counts,
    retain = knee,
    lower = inflection
  )
  drops[[id]]$Rank <- rank(-drops[[id]]$Total)
  drops[[id]]$Retain <- (drops[[id]]$FDR < 0.001)
  drops[[id]]$Col <- 'black'
  drops[[id]]$Col[drops[[id]]$Retain] <- 'Red'
}

#' Cell-barcode calling is usually inspected with a barcode-rank plot, which 
#' orders each cell barcode by the number of UMIs associated with it. Barcodes
#' with more UMIs are assumed to have been from droplet containing true cells
#' than those with very few UMIs.

#+ barcode_rank, fig.height=5, fig.width=8, fig.caption='Barcode-rank plot for each neonatal microglia sample. Red circles denote barcodes considered to be actual cells, whereas black circles denote empty droplets.'
# Plot EmptyDrops output
do_plot <- function() {
  par(mfrow = c(2,3))
  for (id in names(drops)) {
    tmp_uniq <- !duplicated(drops[[id]]$Rank)
    # plot points
    plot(x = drops[[id]][['Rank']][tmp_uniq],
         y = drops[[id]][['Total']][tmp_uniq],
         col = drops[[id]]$Col[tmp_uniq],
         log = 'xy',
         main = id, xlab = 'Rank', ylab = 'Total UMI',
         cex.lab = 1.4, cex.main = 1.4, cex.axis = 1.2)
  }
  par(mfrow = c(1,1))
}
# tiff(file = paste0(results_out, 'barcodeRankPlots.tiff'),
#     height = 5, width = 8, units = 'in')
# do_plot()
# dev.off()
do_plot()

# Get cell barcodes that pass
retain_bc <- vector(mode = 'list', length = length(drops))
names(retain_bc) <- names(drops)
for (id in names(drops)) {
  retain_bc[[id]] <- rownames(drops[[id]])[which(drops[[id]]$Retain)]
}

#+ retention_count
knitr::kable(x = sapply(retain_bc, length))

# Load data, filter retained cells
counts <- vector(mode = 'list', length = length(retain_bc))
names(counts) <- names(retain_bc)
for (id in names(counts)) {
  tmp <- Read10X(data.dir = paste0(data_path, id))
  counts[[id]] <- tmp[,colnames(tmp) %in% retain_bc[[id]]]
}

# Save filtered count matrix
dir.create(path = '../data/counts/')
for (id in names(counts)) {
  tmp_name <- paste0('../data/counts/', id, '_filtered_feature_bc_matrix.rds')
  saveRDS(object = counts[[id]], file = tmp_name)
} 

rm(data_path, drops, counts, ranks, knee, 
   inflection, tmp_uniq, retain_bc, tmp_name)
gc(verbose=FALSE)



# Quality Control ---------------------------------------------------------

# Load counts matrices 
counts_path <- list.files('../data/counts/', full.names = TRUE)
counts <- vector(mode = 'list', length = length(counts_path))
names(counts) <- sample_names
for (id in names(counts)) {
  counts[[id]] <- readRDS(file = counts_path[grepl(id, counts_path)])
}

# Calculate percent.mt and percent.rp
for (id in names(counts)) {
  counts[[id]] <- CreateSeuratObject(
    counts = counts[[id]],
    project = id
  )
  counts[[id]] <- PercentageFeatureSet(
    object = counts[[id]],
    pattern = '^mt-',
    col.name = 'percent.mt'
  )
  counts[[id]] <- PercentageFeatureSet(
    object = counts[[id]],
    pattern = '^Rp[ls]',
    col.name = 'percent.rp'
  )
}


# Extract and plot data
#+ qc_metrics, fig.height=4, fig.width=5.5, fig.cap='Violin plot of various quality control metrics. nCount_RNA = number of UMIs per cell. nFeature_RNA = number of unique genes detected per cell. percent.mt = % of UMIs mapping to mitochondrial genes. percent.rp = % of UMIs mapping to ribosomal genes.
meta_feats <- c('orig.ident','nCount_RNA','nFeature_RNA','percent.mt','percent.rp')
metadata <- lapply(
  X = counts,
  FUN = function(x) {
    x@meta.data[meta_feats]
  }
)
metadata <- Reduce(rbind, metadata)
p1 <- metadata %>%
  reshape2::melt(id.vars = c('orig.ident')) %>%
  ggplot(mapping = aes(x = orig.ident, y = value)) + 
  geom_violin(mapping = aes(fill = orig.ident), scale = 'width') +
  facet_wrap(. ~ variable, scales = 'free_y') +
  theme_bw() +
  xlab(label = 'Sample ID') +
  ylab(label = 'Value') +
  theme(axis.text.x = element_blank())
# tiff(filename = paste0(results_out, 'qualityControlMetrics.tiff'),
#     height = 4, width = 6, units = 'in')
# p1
# dev.off()
p1


#' Calculate median-absolute-deviation based thresholds for nCount_RNA, which 
#' allows for variability in mean sequencing depth by sample which can be 
#' affected by number of cells. Set flat threshold at 15% for percent.mt based
#' on similarity in distribution between samples. Filter cells.  
get_mad_max <- function(x, dev = 3) {
  return(median(x) + dev * mad(x, center = 1))
}

high_qc <- vector(mode = 'list', length = length(counts))
names(high_qc) <- names(counts)
for (id in names(high_qc)) {
  umi_mad_max <- get_mad_max(counts[[id]]$nCount_RNA, dev = 3)
  high_qc[[id]] <- counts[[id]]@meta.data$nCount_RNA <= umi_mad_max &
    counts[[id]]@meta.data$percent.mt <= 15
}


#' Next, we remove putative doublets using the Scrublet python package. 
require('reticulate')
doublet_rate <- read.table(file = '../ref/DoubletRates_10x.tsv', header = TRUE)

# Import package and std_out config details.
scrublet_out <- paste0(results_out, 'scrublet_outs/')
dir.create(path = scrublet_out)
scrub <- import(module = 'scrublet', convert = FALSE)
writeLines(text = str(py_config()), con = paste0(scrublet_out, 'py_config.txt'))

# Run Scrublet (https://doi.org/10.1016/j.cels.2018.11.005)
doublet_results <- vector(mode = 'list', length = length(counts))
names(doublet_results) <- names(counts)
svg(filename = paste0(scrublet_out, 'score_hist.svg'),
    height = 5, width = 9)
par(mfrow = c(2, 3))
for (ii in 1:length(counts)) {
  sample_id <- names(counts)[ii]
  ncells_est <- round(ncol(counts[[sample_id]]), digits = -3)
  rate <- doublet_rate$Multiplet_rate[doublet_rate$nCells_Recovered == ncells_est]
  scrublet_obj <- scrub$Scrublet(
    counts_matrix = r_to_py(
      x = Matrix::t(counts[[sample_id]]))$tocsc(),
    expected_doublet_rate = rate
  )
  scrublet_result <- py_capture_output(
    scrublet_obj$scrub_doublets(
      min_counts = 2,
      min_cells = 3,
      min_gene_variability_pctl = 85,
      verbose = TRUE
    )
  )
  writeLines(
    text = c(sample_id, scrublet_result),
    con = paste0(results_out, 'scrublet_outs/', sample_id, '.txt')
  )
  scores <- py_to_r(scrublet_obj$doublet_scores_obs_)
  threshold <- py_to_r(scrublet_obj$threshold_)
  
  # Plot for inspection
  hist(x = scores, breaks = 100, freq = FALSE,
       main = sample_id, xlab = 'doublet score', ylab = 'prob. density')
  abline(v = threshold, lty = 'dashed')
  
  doublet_results[[sample_id]][['Doublet_score']] <- py_to_r(scrublet_obj$doublet_scores_obs_)
  doublet_results[[sample_id]][['is_doublet']] <- 
    doublet_results[[sample_id]][['Doublet_score']] > py_to_r(scrublet_obj$threshold_)
  names(doublet_results[[sample_id]]) <- c('doublet_score', 'is_doublet')
  doublet_results[[sample_id]][['nCells_estimated']] <- ncells_est
  doublet_results[[sample_id]][['Multiplet_rate_10X']] <- rate
  doublet_results[[sample_id]][['Threshold_score']] <- py_to_r(scrublet_obj$threshold_)
}
dev.off()

# Plot again for pdf output
#+ doublet_scoring, fig.height=5, fig.width=9, fig.cap='Scrublet algorithm summary. Vertical lines denote thresholds used to call doublets.'
par(mfrow = c(2,3))
for (i in 1:length(doublet_results)) {
  hist(x = doublet_results[[i]]$doublet_score, breaks = 100, freq = FALSE,
       main = names(doublet_results)[i], 
       xlab = 'doublet score', 
       ylab = 'prob. density',
       xaxt = 'n')
  axis(side = 1, at = seq(0, 1, 0.05))
  abline(v = doublet_results[[i]]$Threshold_score, lty = 'dashed')
}
par(mfrow = c(1,1))

#' As is apparent in the distributions, using Scrublet with default parameters 
#' fails to identify an acceptable threshold score. Instead, we opt to set the
#' threshold manually. 
doublet_results$GSM4559967_Ctl_1$Threshold_score <- 0.20
doublet_results$GSM4559968_3dpi_1$Threshold_score <- 0.125
doublet_results$GSM4559969_3dpi_2$Threshold_score <- 0.125
doublet_results$GSM4559970_5dpi_1$Threshold_score <- 0.1
doublet_results$GSM4559971_5dpi_2$Threshold_score <- 0.1

for (i in 1:length(doublet_results)) {
  doublet_results[[i]]$is_doublet <- 
    doublet_results[[i]]$doublet_score > doublet_results[[i]]$Threshold_score
}

# Summary of doublet call results
knitr::kable(x = t(sapply(doublet_results, 
                          FUN = function(x) table(x[['is_doublet']]))),
             caption = 'Scrublet summary: is doublet?')


#' We retain all cells that pass QC metrics and doublet detection.  
good_cells <- vector(mode = 'list', length = length(counts))
names(good_cells) <- names(counts)
for (i in 1:length(counts)) {
  good_cells[[i]] <- which(high_qc[[i]] & !doublet_results[[i]]$is_doublet)
}

# Summary table
filter_results <- cbind(
  'Before' = sapply(counts, ncol),
  'Filtered' = sapply(good_cells, length)
)
knitr::kable(x = filter_results, caption = 'Number of cells retained after QC and doublet detection.')

# Filter cells
for (id in names(counts)) {
  counts[[id]] <- counts[[id]][, good_cells[[id]]]
}

dir.create('../data/qc_filtered_counts/')
for (id in names(counts)) {
  saveRDS(counts[[id]][['RNA']]@counts, 
          file = paste0('../data/qc_filtered_counts/', id, '.rds'))
}

rm(id, p1, get_mad_max, high_qc, filter_results, counts, doublet_results, 
   filter_results, good_cells, scrub, scrublet_obj, i, ii, sample_id, 
   ncells_est, scrublet_result, meta_feats, rate, scores, threshold)
gc(verbose=FALSE)




# Normalization testing -------------------------------------------------------


#' ### Log-normalization, column-merge
counts_path <- list.files('../data/qc_filtered_counts', full.names = TRUE)
counts <- vector(mode = 'list', length = length(counts_path))
names(counts) <- sample_names
for (id in sample_names) {
  counts[[id]] <- readRDS(file = counts_path[grepl(id, counts_path)])
}

# Give unique cell barcodes across samples
for (id in names(counts)) {
  colnames(counts[[id]]) <- paste(
    substr(colnames(counts[[id]]), start = 1, stop = 16),
    id,
    sep = '_'
  )
}

#' Initially, we test a column merge before data normalization and downstream. 
#' This is to test whether significant batch effects are present between 
#' samples.  
neo <- Reduce(f = cbind, x = counts)
neo <- CreateSeuratObject(counts = neo)
neo <- NormalizeData(neo, verbose = FALSE)
firstup <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}
neo <- CellCycleScoring(
  object = neo, 
  s.features = firstup(cc.genes.updated.2019$s.genes),
  g2m.features = firstup(cc.genes.updated.2019$g2m.genes)
)
neo$CC.difference <- neo$S.Score - neo$G2M.Score
neo <- FindVariableFeatures(neo, verbose = FALSE)
neo <- ScaleData(neo, vars.to.regress = 'CC.difference', verbose = FALSE)
neo <- RunPCA(neo, npcs = 30, verbose = FALSE)
# ElbowPlot(neo, ndims = 30)
neo <- FindNeighbors(neo, dims = 1:15, verbose = FALSE)
neo <- RunUMAP(neo, dims = 1:15, verbose = FALSE)
neo <- FindClusters(neo, resolution = 0.8, verbose = FALSE)

# Pull sample information from barcodes
tmp_bc <- strsplit(x = colnames(neo), split = '_')
tmp_bc <- sapply(X = tmp_bc, FUN = function(x) paste(x[2:4], collapse = '_'))
neo$sample_id <- tmp_bc
tmp_group <- sapply(strsplit(x = neo$sample_id, split = '_'), `[`, 2)
neo$group <- tmp_group


#+ columnMerge_umap, fig.height=4, fig.width=14.5, fig.cap='UMAP of cells by cluster (left), injury time-point (middle), and sample (right). No batch correction on samples.'
p1 <- DimPlot(neo, group.by = 'seurat_clusters', label = TRUE, label.size = 5,
              shuffle = TRUE) + theme_bw() + NoLegend()
p2 <- DimPlot(neo, group.by = 'group', shuffle = TRUE) + theme_bw()
p3 <- DimPlot(neo, group.by = 'sample_id', shuffle = TRUE) + theme_bw()
summary_umap <- p1 + p2 + p3
summary_umap
# ggsave(filename = paste0(results_out, 'columnMergeLogNormalization_umap.tiff'),
#        plot = summary_umap, height = 4, width = 14.5)



#' ### Log-normalization, batch-corrected
counts_path <- list.files('../data/qc_filtered_counts', full.names = TRUE)
counts <- vector(mode = 'list', length = length(counts_path))
names(counts) <- sample_names
for (id in sample_names) {
  counts[[id]] <- readRDS(file = counts_path[grepl(id, counts_path)])
}

# Give unique cell barcodes across samples
for (id in names(counts)) {
  colnames(counts[[id]]) <- paste(
    substr(colnames(counts[[id]]), start = 1, stop = 16),
    id,
    sep = '_'
  )
}

neo <- list(
  'Ctl' = counts$GSM4559967_Ctl_1,
  'D3' = cbind(counts$GSM4559968_3dpi_1, counts$GSM4559969_3dpi_2),
  'D5' = cbind(counts$GSM4559970_5dpi_1, counts$GSM4559971_5dpi_2)
)
for (i in 1:length(neo)) {
  neo[[i]] <- CreateSeuratObject(
    counts = neo[[i]],
    project = names(neo)[i]
  )
  neo[[i]] <- NormalizeData(neo[[i]], verbose = FALSE)
  neo[[i]] <- CellCycleScoring(
    object = neo[[i]],
    s.features = firstup(cc.genes.updated.2019$s.genes),
    g2m.features = firstup(cc.genes.updated.2019$g2m.genes)
  )
  neo[[i]]$CC.difference <-  neo[[i]]$S.Score - neo[[i]]$G2M.Score
  neo[[i]] <- FindVariableFeatures(neo[[i]], verbose = FALSE)
}

neo_anchors <- FindIntegrationAnchors(object.list = neo, verbose = FALSE)
neo <- IntegrateData(anchorset = neo_anchors, verbose = FALSE)
rm(neo_anchors); gc(verbose = FALSE)
DefaultAssay(neo) <- 'integrated'
neo <- ScaleData(neo, vars.to.regress = 'CC.difference', verbose = FALSE)
neo <- RunPCA(neo, npcs = 40, verbose = FALSE)
ElbowPlot(neo, ndims = 40)
pcs <- 1:15
neo <- FindNeighbors(neo, dims = pcs, verbose = FALSE)
neo <- RunUMAP(neo, dims = pcs, verbose = FALSE)
neo <- FindClusters(neo, resolution = 0.8, verbose = FALSE)

# Pull sample information from barcodes
tmp_bc <- strsplit(x = colnames(neo), split = '_')
tmp_bc <- sapply(X = tmp_bc, FUN = function(x) paste(x[2:4], collapse = '_'))
neo$sample_id <- tmp_bc
tmp_group <- sapply(strsplit(x = neo$sample_id, split = '_'), `[`, 2)
neo$group <- tmp_group
neo$group <- factor(x = neo$group,
                    levels = c('Ctl', '3dpi', '5dpi'))

#+ columnMerge_umap, fig.height=4, fig.width=14.5, fig.cap='UMAP of cells by cluster (left), injury time-point (middle), and sample (right). Batch-corrected across injury group, replicates column-merged.'
p1 <- DimPlot(neo, group.by = 'seurat_clusters', label = TRUE, label.size = 5,
              shuffle = TRUE) + theme_bw() + NoLegend()
p2 <- DimPlot(neo, group.by = 'group', shuffle = TRUE) + theme_bw()
p3 <- DimPlot(neo, group.by = 'sample_id', shuffle = TRUE) + theme_bw()
summary_umap <- p1 + p2 + p3
summary_umap
# ggsave(filename = paste0(results_out, 'integratedLogNormalization_umap.tiff'),
#        plot = summary_umap, height = 4, width = 14.5)





# Cell-type annotation and marker genes -----------------------------------

#' ## Identification of cell-types in neonatal spinal cord

# Add some lost metadata
DefaultAssay(neo) <- 'RNA'
neo <- PercentageFeatureSet(
  object = neo,
  pattern = '^mt-',
  col.name = 'percent.mt'
)
neo <- PercentageFeatureSet(
  object = neo,
  pattern = '^Rp[ls]',
  col.name = 'percent.rp'
)

# ' First, identify differential expressed genes per cluster.
Idents(neo) <- 'seurat_clusters'
markers <- FindAllMarkers(
  object = neo,
  assay = 'RNA',
  only.pos = TRUE,
  logfc.threshold = 0.25
)
top_markers <- markers %>%
  group_by(cluster) %>%
  filter(p_val_adj == 0) %>%
  top_n(n = 3, wt = -avg_logFC)
knitr::kable(x = top_markers)


#' Plot marker genes as used Extended Data Fig 6e.'
prior_genes <- list(
  'microglia1' = c('P2ry12','Tmem119'),
  'microglia2' = c('Cd9','Serpine2'),
  'microglia3' = c('Spp1','Fn1'),
  'microglia4' = c('Ms4a7','Ms4a6c'),
  'div_microglia1' = c('Stmn1','Ube2c'),
  'div_microglia2' = c('Mcm3','Mcm6'),
  'macrophage' = c('Mrc1','Lyve1'),
  'monocyte' = c('Ccr2','Ly6c2'),
  'neutrophil' = c('S100a8','S100a9'),
  'b_cell' = c('Cd74','Ly6d'),
  't_cell' = c('Trbc1','Trbc2'),
  'astrocyte' = c('Sparcl1','Gfap'),
  'oligodendrocyte' = c('Olig1','Olig2')
)
DefaultAssay(neo) <- 'RNA'
gene_umap <- FeaturePlot(
  object = neo,
  features = unlist(prior_genes, use.names = FALSE),
  # order = TRUE, 
  combine = FALSE
)
tmp_names <- rep(names(prior_genes), each = 2)
for (i in 1:length(gene_umap)) {
  gene_umap[[i]] <- gene_umap[[i]] + 
    scale_color_viridis_c(option = 'A') + 
    labs(title = unlist(prior_genes)[i],
         subtitle = tmp_names[i]) +
    theme_bw() +
    NoLegend()
}
gene_umap <- cowplot::plot_grid(
  plotlist = gene_umap,
  ncol = floor(length(prior_genes)/2),
  byrow = TRUE
)
# ggsave(filename = paste0(results_out, 'gene_umap.tiff'),
#        plot = gene_umap, device = 'tiff', height = 14, width = 14)

#+ neonatal_gene_umap, fig.height=14, fig.width=14, fig.cap='Expression of marker genes tabled in Extended Data Figure 6E.'
gene_umap


gene_vln <- FetchData(
  object = neo, vars = c(unlist(prior_genes), 'seurat_clusters'),
  slot = 'data'
) %>%
  reshape2::melt(id.vars = 'seurat_clusters') %>%
  ggplot(mapping = aes(x = seurat_clusters, y = value)) + 
  geom_violin(mapping = aes(fill = seurat_clusters),
              scale = 'width') +
  facet_wrap(. ~ variable, scales = 'free_y', ncol = 4) +
  theme(legend.position = 'none') +
  ylab(label = 'Normalized expression') +
  xlab(label = 'Seurat cluster') +
  theme_bw() +
  theme(strip.text = element_text(size = 12)) +
  NoLegend()
# ggsave(filename = paste0(results_out, 'cluster_marker_gene_vln.tiff'),
#        plot = gene_vln, height = 7.5, width = 18, device = 'tiff')

#+ neonatal_gene_vln, fig.height=7.5, fig.width=18, fig.cap='Violin plots of marker gene expression. Genes taken from Extended Data Figure 6E.'
gene_vln


#' Based on expression, we assign each cluster to the following cell-types:
celltype <- c(
  '0' = 'Microglia',
  '1' = 'Microglia',
  '2' = 'Microglia',
  '3' = 'Microglia',
  '4' = 'Monocyte',
  '5' = 'Microglia',
  '6' = 'Macrophage',
  '7' = 'Microglia',
  '8' = 'Macrophage', # Maybe microglia
  '9' = 'Microglia',
  '10' = 'Neutrophil',
  '11' = 'Microglia',
  '12' = 'Microglia',
  '13' = 'Macrophage', # Maybe microglia
  '14' = 'T cell',
  '15' = 'Microglia',
  '16' = 'Monocyte',
  '17' = 'Microglia',
  '18' = 'Astrocyte',
  '19' = 'B cell',
  '20' = 'Monocyte',
  '21' = 'Oligodendrocyte'
)

neo$celltype <- plyr::mapvalues(
  x = neo$seurat_clusters,
  from = names(celltype),
  to = celltype
)

#+ celltype_umap, fig.height=4, fig.width=5.75, fig.cap='UMAP of cells by annotated cell-type (using marker definitions provided by Li et al, 2020)
p1 <- DimPlot(neo, group.by = 'celltype',
              label = TRUE, label.size = 4, repel = TRUE) +
  theme_bw()
p1
# ggsave(filename = paste0(results_out, 'celltype_umap.tiff'),
#        plot = p1, height = 4, width = 5.75, device = 'tiff')

saveRDS(object = neo, file = '../data/neonatal.rds')



# Reproducting results by Li et al ----------------------------------------

neo <- readRDS(file = '../data/neonatal.rds')

#' ## Reproducing results by Li et al., 2020.  

#' In the article's methods section, they performed quality control by using 
#' "cells with mitochondrial expression level lower than 5%, ribosomal 
#' expression level lower than 30%, and a number of features larger than 2000.  

#+ qc_reproduce, fig.height=3, fig.width=9, fig.cap='Violin plot of QC metrics used by Li et al. Red dashed lines denote thresholds used (min, max, max, respectively).'
p1 <- VlnPlot(
  object = neo, 
  features = c('nFeature_RNA','percent.mt','percent.rp'),
  pt.size = 0, 
  group.by = 'sample_id',
  combine = FALSE)
p1 <- lapply(p1, function(x) {
  x + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank()) +
    NoLegend()
  }
)
p1[[1]] <- p1[[1]] + 
  geom_hline(yintercept = 2000, linetype = 'dashed', color  = 'red', size = 1)
p1[[2]] <- p1[[2]] + 
  geom_hline(yintercept = 5, linetype = 'dashed', color = 'red', size = 1)
p1[[3]] <- p1[[3]] + 
  geom_hline(yintercept = 30, linetype = 'dashed', color = 'red', size = 1)
p1 <- cowplot::plot_grid(plotlist = p1, ncol = 3)
p1

# Produce table of # of cells that qualify.  
knitr::kable(
  x = table('Above unique gene threshold?' = neo$nFeature_RNA > 2000, 
      'Injury time-point' =  neo$orig.ident)
)
knitr::kable(
  x = table('Below mito % threshold?' = neo$percent.mt <= 5,
            'Injury time-point' = neo$orig.ident)
)
knitr::kable(
  x = table('Below ribosomal % threshold?' = neo$percent.rp <= 30,
            'Injury time-point' = neo$orig.ident)
)

#' We see that using the thresholds provided cannot reproduce the number of 
#' that are reported in their Extended Data Fig. 6B. They report 8015 cells
#' from Ctl, 9778 cells from 3dpi, and 7754 cells from 5dpi. If we were to set
#' a cutoff of 5% mitochondrial expression, we would be left with no more than
#' 4000 cells per sample. The drastic difference probably cannot be explained
#' by differences in cell-barcode calling. It is also unclear what mitochondrial
#' genes they used for calculating %.  

knitr::kable(
  x = prop.table(table('Is microglia:' = neo$celltype == 'Microglia',
                       'Injury time-point' = neo$orig.ident), margin = 2),
  caption = 'Proportion of microglia per injury time-point'
)

#' Comparing the table above to Extended Data Fig. 6c, we see general agreement.
#' Microglia proportion drops at 3dpi and increases modestly by 5dpi.  


