
#' ---
#' title: "Integrating adult and neonatal SCI"
#' author: "James Choi"
#' date: "`r Sys.Date()`"
#' output: pdf_document
#' ---


#+ setup, include=TRUE
knitr::opts_chunk$set(warning=FALSE, message=FALSE, fig.align='center', 
                      tidy.opts=list(width.cutoff=60), tidy=TRUE)

#' ## Batch correction between adult SCI myeloid datasets
#' I previously found that linear-regression based methods of batch correction
#' (using the batchelor R package) are able to correct for differences between
#' replicates of our SCI time-points. The injury time-point replicates of the 
#' neonatal SCI dataset produced by Li et al. appear to require no correction
#' based on UMAP visualization.  
#' 
#' In order to create a combined SCI immune cell atlas, I modify the typical 
#' integration workflow in the following ways:  
#' 
#' 1. Apply batchelor::rescaleBatches() between replicates adult SCI time-
#' points. This will generate a single corrected matrix for each adult SCI time-
#' point containing log-normalized expression values.
#' 2. Find variable features in log-expression space(?). FindVariableFeatures()
#' function in Seurat selects features using log(mean) vs log(variance) of raw
#' count values. scran tutorial suggests methods to select features by computing
#' variance of log-normalized expression values. 

# Load libraries
require('Seurat')
require('ggplot2')
require('dplyr')
require('SingleCellExperiment')
require('batchelor')
require('scran')


# Util functions
firstup <- function(x) {
  x <- tolower(x); substr(x,1,1) <- toupper(substr(x,1,1)); return(x)
}

# Directories
if (!grepl('scripts', getwd())) {
  setwd('scripts/')
}
results_out <- '../results/IntegratingAdultNeonatal/'
dir.create(path = results_out)

# Import data
adult <- readRDS(file = '../../sci_scRNAseq/data/myeloid.rds')
neo <- readRDS(file = '../data/neonatal.rds')

# Rename myeloid cells
adult@meta.data[['myeloid_functional']] <- plyr::mapvalues(
  x = adult@meta.data[['integrated_snn_res.0.35']],
  from = levels(adult@meta.data[['integrated_snn_res.0.35']]),
  to = c('Homeostatic Microglia',
         'Inflammatory Microglia',
         'Chemotaxis-Inducing Mac',
         'Dividing Microglia',
         'Inflammatory Mac',
         'Monocyte',
         'Neutrophil',
         'Migrating Microglia',
         'Dendritic',
         'Interferon Myeloid',
         'Dividing Myeloid',
         'Border-Associated Mac')
)
adult@meta.data[['myeloid_functional']] <- factor(
  x = adult@meta.data[['myeloid_functional']],
  levels = c('Neutrophil',
             'Monocyte',
             'Chemotaxis-Inducing Mac',
             'Inflammatory Mac',
             'Border-Associated Mac',
             'Dendritic',
             'Dividing Myeloid',
             'Homeostatic Microglia',
             'Inflammatory Microglia',
             'Dividing Microglia',
             'Migrating Microglia',
             'Interferon Myeloid')
)

# Remove other assay slots
DefaultAssay(adult) <- 'RNA'
DefaultAssay(neo) <- 'RNA'
adult[['RNAcorrected']] <- NULL
adult[['integrated']] <- NULL
neo[['integrated']] <- NULL
neo <- neo[,neo$celltype %in% c('Microglia','Monocyte','Macrophage','Neutrophil')]

# Define common set of genes
universe <- intersect(rownames(adult), rownames(neo))
adult <- adult[universe,]
neo <- neo[universe,]


#' To correct for differences in sequencing depth between replicates of the 
#' adult SCI time-points (Choi and Milich study), I use the rescaleBathces() 
#' function from the batchelor R package. In brief, this method will 1) revert 
#' log-expression values and remove pseudocount, 2) per-gene, scale down counts
#' such that the average of each replicate is equal to the lowest average, 3) 
#' then re-transforming. According to Lun et al, this preserves sparsity of the
#' matrix and mitigates differences in variance due to varying means. The 
#' function below was used for our SCI scSeq manuscript and adapted for this 
#' analysis. Only difference is that the 'corrected_log' Seurat object is built
#' by setting 'data =' argument instead of 'counts ='.

# Inputs:
#   obj_list: list of Seurat objects. These should be technical replicates of a 
#     given condition/group.
correct_replicates <- function(obj_list) {
  
  # Identify shared genes across datasets and take common subset
  gene_universe <- lapply(X = obj_list, 
                          FUN = function(x) {
                            return(rownames(slot(x[['RNA']], 'data')))
                          }
  )
  gene_universe <- Reduce(intersect, gene_universe)
  obj_list_sce <- lapply(X = obj_list, 
                         FUN = function(x) {
                           x <- as.SingleCellExperiment(x)
                           x <- x[gene_universe,]
                           return(x)
                         }
  )
  
  # Extract raw count and cell-level meta data. Merge all into
  # single matrix. These values derived directly from Seurat objects.
  raw_counts <- do.call(cbind, lapply(X = obj_list_sce, FUN = counts))
  log_raw_counts <- do.call(cbind, lapply(X = obj_list_sce, FUN = logcounts))
  cell_metadata <- do.call(rbind, lapply(X = obj_list_sce, FUN = colData))
  cell_metadata <- data.frame(cell_metadata)
  
  # Perform linear regression-based batch correction (ie scale down counts).
  # Automatically generates single, merged matrix. Recalculate the corrected 
  # "raw counts".
  corrected_log <- batchelor::rescaleBatches(
    obj_list_sce, 
    log.base = exp(1),
    pseudo.count = 1
  )
  corrected_log <- Matrix::Matrix(
    data = assays(corrected_log)[['corrected']],
    sparse = TRUE
  )
  
  # Set corrected values as new Seurat Assay slot
  corrected_log <- CreateAssayObject(data = corrected_log)
  out_seurat <- CreateSeuratObject(counts = raw_counts,
                                   assay = 'RNA')
  out_seurat[['RNAcorrected']] <- corrected_log
  
  # Import log counts and metadata. Return assembled Seurat object.
  slot(out_seurat[['RNA']], 'data') <- log_raw_counts
  out_seurat@meta.data <- cell_metadata
  return(out_seurat)
}


# Set Seurat objects and result lists by time-point
neo <- SplitObject(neo, split.by = 'orig.ident')
Idents(adult) <- 'time'
inj_groups <- levels(adult$time)
adult_corrected <- vector(mode = 'list', length = length(inj_groups))
names(adult_corrected) <- inj_groups

# Perform the correction
for(i in 1:length(inj_groups)) {
  inj_time <- inj_groups[i]
  inj_cells <- unlist(x = CellsByIdentities(adult, idents = inj_time),
                      use.names = FALSE)
  inj_cells <- subset(adult, cells = inj_cells)
  inj_cells <- SplitObject(inj_cells, split.by = 'sample_id')
  inj_cells <- correct_replicates(obj_list = inj_cells)
  adult_corrected[[inj_time]] <- inj_cells
  DefaultAssay(adult_corrected[[inj_time]]) <- 'RNAcorrected'
  # message(paste('Done with:', inj_time))
}
rm(inj_time, inj_cells)


#' In normal scSeq workflows, I perform feature selection to improve biological
#' signal:technical noise ratio and computation speed. By default in Seurat,
#' this is done by modeling the log(variance of counts) ~ log(mean of counts)
#' per gene, fitting a local polynomial regression (loess), standardizing gene
#' count matrix by gene means and variances, then computing standardized 
#' variance for each gene. Genes with greatest variance are used for downstream.
#'   
#' Alternatively, I can model the mean(log-expression) ~ var(log-expression)
#' as is demonstrated in [OSCA guidebook](http://bioconductor.org/books/release/OSCA/feature-selection.html#quantifying-per-gene-variation). The procedure is approximately the same, with the exception
#' that Seurat models gene-variance relationship of the counts while OSCA models
#' relationship for log-expression values. I opt for this method because it 
#' allows us to use the batchelor-normalized expression values from before (
#' which does not recompute a raw "counts" matrix).  


#+ feature selection, fig.height=5.5, fig.width=11, fig.cap='Modeling mean-variance relationships for each dataset."
neo_var_models <- vector(mode = 'list', length = length(neo))
names(neo_var_models) <- names(neo)
adult_var_models <- vector(mode = 'list', length = length(adult_corrected))
names(adult_var_models) <- names(adult_corrected)

tiff(filename = '../results/IntegratingAdultNeonatal/feature_selection.tiff',
     height = 6, width = 12, res = 440, units = 'in')
par(mfrow = c(2,4))
for (i in 1:length(adult_corrected)) {
  tmp_sce <- as.SingleCellExperiment(
    x = adult_corrected[[i]],
    assay = 'RNAcorrected'
  )
  adult_var_models[[i]] <- modelGeneVar(logcounts(tmp_sce))
  tmp_fit <- metadata(adult_var_models[[i]])
  plot(x = tmp_fit$mean, y = tmp_fit$var,
       xlab = 'Mean of log-expression',
       ylab = 'Variance of log-expression',
       main = names(adult_corrected)[i])
  curve(tmp_fit$trend(x), col="dodgerblue", add=TRUE, lwd=2)
}
for (i in 1:length(neo)) {
  tmp_sce <- as.SingleCellExperiment(
    x = neo[[i]],
    assay = 'RNA'
  )
  neo_var_models[[i]] <- modelGeneVar(logcounts(tmp_sce))
  tmp_fit <- metadata(neo_var_models[[i]])
  plot(x = tmp_fit$mean, y = tmp_fit$var,
       xlab = 'Mean of log-expression',
       ylab = 'Variance of log-expression',
       main = names(neo)[i])
  curve(tmp_fit$trend(x), col="dodgerblue", add=TRUE, lwd=2)
}
dev.off()
par(mfrow = c(1,1)); rm(list = ls()[grepl('tmp', ls())])


#' Based on the trend-fitting, it looks like the Li et al neonatal samples have
#' smoother curves than our adult samples. Whether this is a reliable proxy for
#' sample quality is uncertain. Nevertheless, I now combine the model results
#' into a single DataFrame and pull the top 3000 variable genes.  

combined_models <- combineVar(adult_var_models[['Uninjured']], 
                              adult_var_models[['1dpi']],
                              adult_var_models[['3dpi']],
                              adult_var_models[['7dpi']],
                              neo_var_models[['Ctl']],
                              neo_var_models[['D3']],
                              neo_var_models[['D5']])
hvg <- getTopHVGs(combined_models, n = 3000)


#' Inspecting the top 20 of these, I see some usual suspects: Spp1, Lyz2, C1qa, 
#' Lgals3, etc. Of note, I also see Xist and Fos. Xist is a sex-specific gene 
#' (see [here](https://www.genecards.org/cgi-bin/carddisp.pl?gene=XIST)), while 
#' Fos has been shown to be upregulated in microglia due to single
#'  cell dissociation (see [work from Beth Stevens & Evan Macosko labs](https://www.biorxiv.org/content/10.1101/2020.12.03.408542v1)). In downstream, I will have to check whether cells separate by these
#' genes.  

#' Finally, I integrate the datasets using these genes. Unlike the previous 
#' SCI myeloid analysis, here I use reciprocal PCA ("RPCA") as implemented in 
#' Seurata. My primary reason was memory and speed efficiency. Seurat devs 
#' recommend RPCA when cell-types between datasets strongly differ, data come 
#' from same platform, or there are many cells. With drastic composition 
#' difference between injured and uninjured adult samples due to peripheral 
#' myeloid infiltration, all three points are met and RPCA seems appropriate.  
sci <- c(unlist(neo), unlist(adult_corrected))
sci <- lapply(
  X = sci,
  FUN = function(x) {
    x <- ScaleData(x, features = hvg_names)
    x <- RunPCA(x, features = hvg_names)
  }
)
sci <- FindIntegrationAnchors(
  object.list = sci,
  anchor.features = hvg_names,
  reduction = 'rpca'
)
sci <- IntegrateData(anchorset = sci)


# Metadata cleanup
ccdiff <- sci$CC.difference
ccdiff[is.na(ccdiff)] <- sci$CC.Difference[!is.na(sci$CC.Difference)]
any(is.na(ccdiff)) | any(names(ccdiff) != rownames(sci@meta.data)) # if fine, F
sci$CC_difference <- ccdiff
percent_mt <- sci$percent_mt * 100
percent_mt[is.na(percent_mt)] <- sci$percent.mt[!is.na(sci$percent.mt)]
any(is.na(percent_mt)) | any(names(percent_mt) != rownames(sci@meta.data))
sci$percent_mt <- percent_mt
percent_rp <- sci$percent_rp * 100
percent_rp[is.na(percent_rp)] <- sci$percent.rp[!is.na(sci$percent.rp)]
any(is.na(percent_rp)) | any(names(percent_rp) != rownames(sci@meta.data))
sci$percent_rp <- percent_rp
remove_cols <- c('integrated_snn_res.0.8','seurat_clusters','group',
                 'nCount_SCT','nFeature_SCT','library_size','pass_umi',
                 'n_genes','pass_n_genes','pass_percent_mt','pass_percent_rp',
                 'percent_hbb','pass_percent_hbb','default_cluster',
                 'integrated_snn_res.0.35','CC.Difference','CC.difference',
                 'percent.mt','percent.rp','is_doublet',
                 'default_myeloid_subcluster','integrated_snn_res.0.4')
sci@meta.data[remove_cols] <- NULL
sci$time[!is.na(sci$time)] <- paste('adult', sci$time[!is.na(sci$time)], sep = '_')
sci$time[is.na(sci$time)] <- sci$orig.ident[is.na(sci$time)]
sci$time <- plyr::mapvalues(
  x = sci$time,
  from = c('Ctl','D3','D5'),
  to = c('neonatal_Uninjured','neonatal_3dpi','neonatal_5dpi')
)
sci$time <- factor(
  x = sci$time,
  levels = c('neonatal_Uninjured','neonatal_3dpi','neonatal_5dpi',
             'adult_Uninjured','adult_1dpi','adult_3dpi','adult_7dpi')
)


# Dimensional reduction
DefaultAssay(sci) <- 'integrated'
sci <- ScaleData(sci, vars.to.regress = 'CC_difference', verbose = FALSE)
sci <- RunPCA(sci, npcs = 50, verbose = FALSE)
ElbowPlot(sci, ndims = 50)
sci <- RunUMAP(sci, dims = 1:15, verbose = FALSE)
sci <- FindNeighbors(sci, dims = 1:15, verbose = FALSE)
sci <- FindClusters(sci, resolution = 0.8, verbose = FALSE)


#' To confirm successful integration, I plot UMAPs and label cells by study, 
#' injury time-point, and previously annotated cell-types.  
study <- c(rep('neonatal', 5), rep('adult', 10))
names(study) <- unique(sci$sample_id)
sci$study <- plyr::mapvalues(
  x = sci$sample_id,
  from = names(study),
  to = study
)
inj_group <- strsplit(
  x = sci$sample_id,
  split = '_'
)
inj_group <- sapply(
  X = inj_group,
  FUN = function(x) {rev(x)[2]}
)
inj_group <- paste(sci$study, inj_group, sep = '_')
inj_group <- gsub(pattern = 'Ctl', replacement = 'uninj', x = inj_group)
sci$inj_group <- factor(
  x = inj_group,
  levels = c('adult_uninj','adult_1dpi','adult_3dpi','adult_7dpi','neonatal_uninj','neonatal_3dpi','neonatal_5dpi')
)

p1 <- DimPlot(sci, group.by = 'study', shuffle = TRUE) + theme_bw()
p2 <- DimPlot(sci, group.by = 'myeloid_subcluster', shuffle = TRUE) + theme_bw()
p3 <- DimPlot(sci, group.by = 'inj_group', split.by = 'inj_group', ncol = 4) + 
  theme_bw() + NoLegend()
p4 <- (p1 + p2) / p3
p4 <- cowplot::plot_grid(cowplot::plot_grid(p1, p2, ncol = 2), p3, ncol = 1)
ggsave(filename = paste0(results_out, 'UMAP_summary.tiff'), plot = p4,
       height = 8, width = 10, device = 'tiff')


saveRDS(sci, file = '../data/myeloid_combined.rds')




# Alternative approaches --------------------------------------------------


# Convert to sce to use scran/batchelor fxns
adult <- as.SingleCellExperiment(adult, assay = 'RNA')
neo <- as.SingleCellExperiment(neo, assay = 'RNA')

# Define common set of genes
universe <- intersect(rownames(adult), rownames(neo))
adult <- adult[universe,]
neo <- neo[universe,]

tmp <- multiBatchNorm(adult, batch = adult$orig.ident)
tmp <- multiBatchNorm(adult)
adult <- SplitObject(adult, split.by = 'sample_id')



# Feature selection via loess(log10(mean) ~ log10(var)) -------------------

#+ mean_var, fig.height=3, fig.width=3.5, fig.cap='log10(variance(counts)) vs log10(mean(counts))'
tmp <- adult_corrected[['3dpi']]
DefaultAssay(tmp) <- 'RNA'
tmp <- FindVariableFeatures(tmp, verbose = FALSE)
head(tmp@assays$RNA@meta.features)
gene_means <- Matrix::rowMeans(tmp[['RNA']]@counts)
rowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}
# microbenchmark::microbenchmark(rowVar(tmp[['RNA']]@counts), 
#                                apply(tmp[['RNA']]@counts,1,var),
#                                times = 3)
# Unit: seconds
# expr       min        lq     mean    median        uq      max
# rowVar(tmp[["RNA"]]@counts)  6.307895  6.536965  7.07475  6.766034  7.458177  8.15032
# apply(tmp[["RNA"]]@counts, 1, var) 34.850113 35.324619 35.50165 35.799124 35.827418 35.85571
gene_vars <- rowVar(tmp[['RNA']]@counts)
plot(log10(gene_means), log10(gene_vars),
     xlab = 'Log10(count mean)', ylab = 'Log10(count variance)')

hvg_info <- data.frame(mean = gene_means, variance = gene_vars)
hvg_info$variance_expected <- 0
hvg_info$variance_standardized <- 0
not_const_var <- which(hvg_info$variance > 0)
loess_model <- loess(
  formula = log10(variance) ~ log10(mean),
  data = hvg_info[not_const_var,],
  span = 0.3
)
hvg_info$variance_expected[not_const_var] <- 10^(loess_model$fitted)
standardizeVar <- function(x, mu, sd, vmax) {
  std <- sqrt(sd)
  for (i in 1:nrow(x)) {
    x[i,] <- min((x[i,] - mu[i])/std[i], vmax)
  }
  outs <- apply(X = x, FUN = var, MARGIN = 1)
  return(outs)
}
hvg_info$variance_standardized[not_const_var] <- standardizeVar(
  x = tmp[['RNA']]@counts,
  mu = gene_means,
  sd = gene_vars,
  vmax = 10
)
