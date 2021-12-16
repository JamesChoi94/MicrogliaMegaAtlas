
#' ---
#' title: "Cluster Analysis of SCI Immune Cells"
#' author: "James Choi"
#' date: "Last run: `r Sys.Date()`"
#' output: pdf_document
#' urlcolor: blue
#' ---


#+ setup, include=TRUE
knitr::opts_chunk$set(warning=FALSE, message=FALSE, fig.align='center', 
                      tidy.opts=list(width.cutoff=60), tidy=TRUE)

#' ## Cluster Analysis of SCI Myeloid cells
#' The goal of this analysis is to identify subsets of immune cells after SCI 
#' and to determine which are shared or unique to the injury responses between
#' neonatal and adult mouse models.  

if (!grepl('scripts', getwd())) {
  setwd('./scripts/')
}
require('Seurat')
require('ggplot2')
require('dplyr')
require('cowplot')

results_out <- '../results/ClusterAnalysisAllMyeloid/'
dir.create(results_out)

sci <- readRDS(file = '../data/myeloid_combined.rds')

#' First, we cluster the cells. We use the graph-based louvain clustering
#' approach that is implemented in Seurat. We use default settings and elect to
#' merge non-distinct or consider finer clustering after intial identification
#' of cell-types and subsets.  

#+ myeloid_umap, fig.height=6, fig.width=9, fig.cap='UMAP of SCI immune cells. Top left: cells colored by cluster. Top right: cells colored by study of origin. Bottom: cells colored by cluster and split by injury time-point across both studies.'
DefaultAssay(sci) <- 'integrated'
sci <- FindClusters(sci, resolution = 0.8, verbose = FALSE)
p1 <- DimPlot(sci, label = TRUE, label.size = 5) + 
  theme_bw() + NoLegend()
p2 <- DimPlot(sci, group.by = 'study', shuffle = TRUE) +
  theme_bw() +
  theme(title = element_blank())
p3 <- DimPlot(sci, split.by = 'time') + 
  theme_bw() + NoLegend()
p4 <- plot_grid(plot_grid(p1, p2, ncol = 2, rel_widths = c(0.8,1)), p3,
                ncol = 1, rel_heights = c(2,1))
ggsave(filename = paste0(results_out, 'allMyeloid_cluster-summary_umap.tiff'),
       plot = p4, height = 6, width = 9, device = 'tiff')
p4


#' Inspecting the UMAPs, we see that there are some age-specific immune cells as
#' well as subsets that are enriched in neonatal vs. adult. For example, cluster
#' 13 and cluster 17 are specific to adult and neonatal, respectively, while 
#' many more cells in cluster 12 come from the neonatal dataset. Cluster 10
#' appears to be comprised almost entirely from neonatal dataset. To better 
#' characterize these clusters, we next identify the globally differentially
#' expressed genes per cluster and examine the top few.  

#+ eval=FALSE
markers <- FindAllMarkers(
  object = sci, 
  assay = 'RNA',
  slot = 'data', 
  test.use = 'wilcox',
  logfc.threshold = log(2),
  # only.pos = TRUE
)
write.table(
  x = markers, 
  file = paste0(results_out, 'defaultClusterMarkers_minFClog2.txt'), 
  sep = '\t',
  quote = FALSE,
  row.names = FALSE
)

#+ markers,
markers <- read.table(
  file = paste0(results_out, 'defaultClusterMarkers_minFClog2.txt'),
  header = TRUE
)
top_markers <- markers %>%
  group_by(cluster) %>%
  filter(avg_logFC > 0) %>%
  top_n(n = 3, wt = -p_val_adj) %>%
  top_n(n = 3, wt = avg_logFC)
knitr::kable(
  x = top_markers, 
  caption = 'Top 3 differentially expressed genes per cluster'
)

#' Because this test is a global comparison, it would benefit to generate violin
#' plots to visualize distribution of gene expression.  

#+ top_gene_vln, fig.height=20, fig.width=6, fig.cap='Violin plot of top 2 marker genes per cluster.' 
top_gene <- markers %>%
  group_by(cluster) %>%
  filter(avg_logFC > 0) %>%
  top_n(n = 2, wt = -p_val_adj) %>%
  top_n(n = 2, wt = avg_logFC)
top_gene <- unique(top_gene$gene)

sci$seurat_clusters <- sci$integrated_snn_res.0.8
Idents(sci) <- 'seurat_clusters'
DefaultAssay(sci) <- 'RNA'
expr_dat <- FetchData(
  object = sci,
  vars = c('seurat_clusters', top_gene),
  slot = 'data'
)
top_gene_vln <- expr_dat %>%
  reshape2::melt(id.vars = 'seurat_clusters') %>%
  ggplot(mapping = aes(x = seurat_clusters, y = value)) +
  geom_violin(mapping = aes(fill = seurat_clusters), scale = 'width') +
  facet_grid(variable ~ ., scales = 'free_y') +
  xlab(label = 'Cluster') + 
  ylab(label = 'log-normalized expression') + 
  theme_bw() + 
  theme(legend.position = 'none')
ggsave(filename = paste0(results_out, 'allMyeloid_default-cluster_top2genes_vln.tiff'),
       plot = top_gene_vln, device = 'tiff', height = 20, width = 6)
top_gene_vln


#' Based on the previous marker gene violin plots and UMAP, clusters 3, 6, 8, 9,
#' 13, 15, 16, and 17 correspond to peripherally-derived myeloid cells. Cluster 
#' 13 and 17 are subsets of neutrophils, but what distinguishes them needs to be
#' further investigated. Cluster 8 is Ccr2+/Ly6c2+ monocyte. Cluster 9 is 
#' macrophages early in differentiation from monocytes and express Arg1 and Fn1. 
#' Accordinly, cluster 3 is likely a more "mature" macrophage in transition to 
#' becoming a foamy macrophage as seen by expression of Fabp4, Fabp5, and Plin2.
#' Clusters 9 and 3 together likely comprise what we' previously described as 
#' macrophage-A. Not surprisingly, cluster 6 likely corresponds to 
#' macrophage-B. It is enriched in Gpnmb, Apoe, Ms4a7, and Lyz2. Cluster 15 
#' highly expresses Cd74 and other MHC-class related genes. These are likely 
#' what we previously described as Dendritic cells. Cluster 16 expresses cell
#' division genes.  
#' 
#' Clusters 0, 1, 2, 4, 5, 7, 11, 12, and 18 are likely different subsets of 
#' microglia. Upregulation of Tnf and Cd83 by cluster 0 suggests an inflammatory
#' and/or "activated" signature. P2ry12 is enriched in cluster 1, indicating a
#' homeostatic subset. Cluster 2 is very similar to cluster 0 with few DE genes
#' - possibly separated due to library size differences. Cluster 4 shows 
#' upregulation of Igf1, Lpl, and Spp1, which we previously identified to be
#' a potentially [regeneration permissive subset of microglia](https://pubmed.ncbi.nlm.nih.gov/30705270/)
#' peaking later at 7dpi in adult SCI. Clusters 5 and 12 highly express genes
#' related to proliferation such as Mki67 and Cenpa. Cluster 11 is enriched for
#' genes such as Nme2, Ran, and ribosomal genes and is low in expression of
#' Cx3cr1. Malat1 and mithondrial genes are upregulated in cluster 18, which 
#' means it is probably a low-quality cluster that should be removed. Cluster 7
#' expresses immediate-early genes such as Jun, Klf2, and Fos. This signature is
#' similar to the ex-DAMs found in a [report describing dissociation-induced effects on glial cells](https://www.biorxiv.org/content/10.1101/2020.12.03.408542v1).  
#' 
#' Cluster 19 expresses Igkc, which is a B cell marker. Cluster 14 is an immune 
#' subset enriched for interferon-pathway genes, which has been described in 
#' [other studies](https://pubmed.ncbi.nlm.nih.gov/29020624/) but their role in 
#' SCI is unclear . Cluster 10 are likely border-associated macrophages as 
#' described by  [Van Hove et al 2018](https://pubmed.ncbi.nlm.nih.gov/31061494/).
#'   

#' After identifying marker genes per cluster, we can also identify how much of
#' each cluster comprise each of the samples. We can further cluster the samples
#' by similarity of cell-type composition to understanding similarities and 
#' differences in total immune cell response to SCI in adult vs neonatal.  
#+ cell_proportion_heatmap, fig.height=6, fig.width=6, fig.cap='Heatmap of proportion of cells from each cluster by sample.'
clust_prop <- proportions(
  x = table(sci$integrated_snn_res.0.8, sci$sample_id), 
  margin = 2
) * 100
attributes(clust_prop)$class <- 'matrix'
clust_prop_heatmap <- ComplexHeatmap::Heatmap(
  matrix = clust_prop,
  col = circlize::colorRamp2(
    breaks = c(0, 5, 10, 20, 30, 60),
    colors = viridis::inferno(n = 6)
  ),
  clustering_method_rows = 'ward.D2',
  clustering_method_columns = 'ward.D2',
  row_title = 'Cluster Identity',
  column_title = 'Sample ID',
  heatmap_legend_param = list(
    title = 'Percent of sample'
  ),
  column_dend_height = unit(15, units = 'mm')
)
tiff(filename = paste0(results_out, 'allMyeloid_cluster-proportion-by-sample_heatmap.tiff'),
     res = 440, height = 6, width = 6, units = 'in')
clust_prop_heatmap
invisible(dev.off())
clust_prop_heatmap


#' The cluster dendrogram of the samples shows three larger groupings:  
#' 
#' * uninj_sample1 + uninj_sample2
#' * uninj_sample3 + neonatal datasets
#' * all other adult SCI datasets
#' 
#' The groupings likely reflect that the adult SCI datasets contain peripherally
#' derived immune cells while the neonatal data had them removed via FACS. 
#' Interesting to note that uninj_sample3 has low proportion of cells from 
#' cluster 1, which we described as the homeostatic microglia cluster. This 
#' might speak to sample quality. Regardless, because of the neonatal FACS, a 
#' better comparison of proportions would be among the microglia clusters.  

#+ microglia_proportion_heatmap, fig.height=5, fig.width=6, fig.cap='Heatmap of proportion of cells from each microglia cluster by sample.'
micro_sub <- c(0, 1, 2, 4, 5, 7, 11, 12, 18)
micro_sub <- sci@meta.data[c('integrated_snn_res.0.8', 'sample_id')] %>%
  filter(integrated_snn_res.0.8 %in% c(micro_sub))
micro_sub$integrated_snn_res.0.8 <- as.character(micro_sub$integrated_snn_res.0.8)
micro_prop <- proportions(
  x = table(micro_sub$integrated_snn_res.0.8, micro_sub$sample_id),
  margin = 2
) * 100
attributes(micro_prop)$class <- 'matrix'
micro_prop_heatmap <- ComplexHeatmap::Heatmap(
  matrix = micro_prop,
  col = circlize::colorRamp2(
    breaks = c(0, 10, 20, 30, 60),
    colors = viridis::inferno(n = 5)
  ),
  clustering_method_rows = 'ward.D2',
  clustering_method_columns = 'ward.D2',
  row_title = 'Cluster Identity',
  column_title = 'Sample ID',
  heatmap_legend_param = list(
    title = 'Percent of sample'
  ),
  column_dend_height = unit(15, units = 'mm')
)
tiff(filename = paste0(results_out, 'microgliaClusters_cluster-proportion-by-sample_heatmap.tiff'),
     res = 440, height = 5, width = 6, units = 'in')
micro_prop_heatmap
invisible(dev.off())
micro_prop_heatmap


#' 

# saveRDS(sci, file = '../data/myeloid_combined.rds')

sessionInfo()
