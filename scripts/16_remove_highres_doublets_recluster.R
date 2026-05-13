#!/usr/bin/env Rscript

source(file.path("config", "analysis_config.R"))
source(file.path("R", "project_helpers.R"))

suppressPackageStartupMessages({
  library(ArchR)
})

set.seed(random_seed)
setup_archr_session()

review_resolution <- 1.4
review_cluster_col <- "Clusters_Combined_highres_r1p4"
source_reduced_dims <- "LSI_Combined"
doublet_clusters_to_remove <- c("C7", "C8", "C9")

ensure_dirs(file.path(figure_dir, "umap"), file.path(figure_dir, "combined_umap"), log_dir)

proj <- loadArchRProject(path = output_dir)
cell_metadata_before <- as.data.frame(getCellColData(proj))

if (!source_reduced_dims %in% names(proj@reducedDims)) {
  stop("Missing reducedDims `", source_reduced_dims, "`. Run scripts 04 and 12 first.")
}
if (!"GeneExpressionMatrix" %in% getAvailableMatrices(proj)) {
  stop("Missing GeneExpressionMatrix. Run script 04 before reclustering.")
}

message(
  "Recreating high-resolution combined clusters at resolution ", review_resolution,
  " to identify doublet clusters: ", paste(doublet_clusters_to_remove, collapse = ", ")
)
proj <- addClusters(
  input = proj,
  reducedDims = source_reduced_dims,
  method = "Seurat",
  name = review_cluster_col,
  resolution = review_resolution,
  force = TRUE
)

cell_metadata_with_review <- as.data.frame(getCellColData(proj))
exclude_cells <- rownames(cell_metadata_with_review)[
  cell_metadata_with_review[[review_cluster_col]] %in% doublet_clusters_to_remove
]

if (length(exclude_cells) == 0) {
  stop(
    "No cells matched high-resolution clusters to remove: ",
    paste(doublet_clusters_to_remove, collapse = ", ")
  )
}

message("Removing ", length(exclude_cells), " cells before reclustering.")

removed_metadata <- cell_metadata_with_review[
  rownames(cell_metadata_with_review) %in% exclude_cells,
  ,
  drop = FALSE
]
write.csv(
  removed_metadata,
  file = file.path(output_dir, "removed_highres_doublet_cells_r1p4_C7_C8_C9.csv"),
  row.names = TRUE
)

removed_counts <- as.data.frame(table(
  removed_highres_cluster = removed_metadata[[review_cluster_col]],
  sample_id = removed_metadata$Sample
), stringsAsFactors = FALSE)
colnames(removed_counts) <- c("removed_highres_cluster", "sample_id", "n_cells_removed")
write.csv(
  removed_counts,
  file = file.path(output_dir, "removed_highres_doublet_counts_r1p4_C7_C8_C9.csv"),
  row.names = FALSE
)

keep_cells <- setdiff(rownames(cell_metadata_with_review), exclude_cells)
proj <- proj[keep_cells]

proj <- addIterativeLSI(
  ArchRProj = proj,
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

proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8,
  force = TRUE
)

proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine",
  force = TRUE
)

lsi_rna_args <- list(
  ArchRProj = proj,
  useMatrix = "GeneExpressionMatrix",
  depthCol = "Gex_nUMI",
  name = "LSI_RNA",
  iterations = 2,
  clusterParams = list(
    resolution = c(0.2),
    sampleCells = 10000,
    n.start = 10
  ),
  varFeatures = rna_lsi_var_features,
  firstSelection = "variable",
  filterQuantile = rna_lsi_filter_quantile,
  binarize = FALSE,
  dimsToUse = 1:30,
  force = TRUE
)

if (!is.null(rna_lsi_total_features)) {
  lsi_rna_args$totalFeatures <- as.integer(rna_lsi_total_features)
}

proj <- do.call(addIterativeLSI, lsi_rna_args)

proj <- addCombinedDims(
  ArchRProj = proj,
  reducedDims = c("IterativeLSI", "LSI_RNA"),
  name = "LSI_Combined"
)

proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "LSI_RNA",
  name = "UMAP_RNA",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine",
  force = TRUE
)

proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "LSI_Combined",
  name = "UMAP_Combined",
  nNeighbors = 30,
  minDist = 0.4,
  metric = "cosine",
  force = TRUE
)

proj <- addClusters(
  input = proj,
  reducedDims = "LSI_RNA",
  method = "Seurat",
  name = "Clusters_RNA",
  resolution = 0.8,
  force = TRUE
)

proj <- addClusters(
  input = proj,
  reducedDims = "LSI_Combined",
  method = "Seurat",
  name = "Clusters_Combined",
  resolution = 1.0,
  force = TRUE
)

p_umap_sample <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "UMAP"
)
p_umap_cluster <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "Clusters",
  embedding = "UMAP"
)
plotPDF(
  p_umap_sample,
  p_umap_cluster,
  name = "umap_sample_clusters.pdf",
  ArchRProj = proj,
  addDOC = FALSE,
  width = 6,
  height = 6
)
file.copy(file.path(output_dir, "Plots", "umap_sample_clusters.pdf"), file.path(figure_dir, "umap"), overwrite = TRUE)

p_rna_sample <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "UMAP_RNA"
)
p_combined_sample <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "UMAP_Combined"
)
p_combined_cluster <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "Clusters_Combined",
  embedding = "UMAP_Combined"
)
plotPDF(
  p_rna_sample,
  p_combined_sample,
  p_combined_cluster,
  name = "rna_combined_umap.pdf",
  ArchRProj = proj,
  addDOC = FALSE,
  width = 6,
  height = 6
)
file.copy(file.path(output_dir, "Plots", "rna_combined_umap.pdf"), file.path(figure_dir, "combined_umap"), overwrite = TRUE)

cell_metadata_after <- as.data.frame(getCellColData(proj))
write.csv(cell_metadata_after, file = file.path(output_dir, "cell_metadata_after_rna_combined.csv"), row.names = TRUE)

cell_counts_summary <- rbind(
  count_cells_by_sample(cell_metadata_before, "before_highres_doublet_exclusion"),
  count_cells_by_sample(removed_metadata, "removed_highres_r1p4_C7_C8_C9"),
  count_cells_by_sample(cell_metadata_after, "after_highres_doublet_exclusion_reclustered")
)
cell_counts_summary <- rbind(
  cell_counts_summary,
  aggregate(n_cells ~ stage, data = cell_counts_summary, sum) |>
    transform(sample_id = "Total") |>
    subset(select = c("stage", "sample_id", "n_cells"))
)
write.csv(
  cell_counts_summary[order(cell_counts_summary$stage, cell_counts_summary$sample_id), ],
  file = file.path(output_dir, "cell_counts_after_highres_doublet_exclusion_reclustered.csv"),
  row.names = FALSE
)

write_filter_settings(
  path = file.path(output_dir, "highres_doublet_exclusion_reclustering_settings.csv"),
  settings = list(
    review_cluster_col = review_cluster_col,
    review_resolution = review_resolution,
    removed_clusters = paste(doublet_clusters_to_remove, collapse = ","),
    final_combined_cluster_resolution = 1.0,
    atac_lsi_recomputed = TRUE,
    rna_lsi_recomputed = TRUE,
    combined_lsi_recomputed = TRUE,
    note = "Cells assigned to high-resolution combined clusters C7, C8, and C9 at resolution 1.4 were removed as likely doublets, then ATAC/RNA/combined embeddings and clusters were recomputed at final combined resolution 1.0."
  )
)

saveArchRProject(ArchRProj = proj, outputDirectory = output_dir, load = FALSE)

message("Reclustered project after high-resolution doublet exclusion.")
message("Removed cells: ", length(exclude_cells))
message("Remaining cells: ", nrow(cell_metadata_after))
