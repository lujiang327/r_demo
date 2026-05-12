#!/usr/bin/env Rscript

source(file.path("config", "analysis_config.R"))
source(file.path("R", "project_helpers.R"))

suppressPackageStartupMessages({
  library(ArchR)
})

set.seed(random_seed)
setup_archr_session()

clusters_to_remove <- c("C1", "C10")
cluster_col <- "Clusters_Combined"

ensure_dirs(file.path(figure_dir, "umap"), file.path(figure_dir, "combined_umap"), log_dir)

proj <- loadArchRProject(path = output_dir)
cell_metadata_before <- as.data.frame(getCellColData(proj))

if (!cluster_col %in% colnames(cell_metadata_before)) {
  stop("Missing `", cluster_col, "`. Run script 04 before removing clusters.")
}
if (!"GeneExpressionMatrix" %in% getAvailableMatrices(proj)) {
  stop("Missing GeneExpressionMatrix. Run script 04 before removing clusters.")
}

exclude_cells <- rownames(cell_metadata_before)[cell_metadata_before[[cluster_col]] %in% clusters_to_remove]
if (length(exclude_cells) == 0) {
  warning("No cells matched clusters_to_remove: ", paste(clusters_to_remove, collapse = ", "))
} else {
  message(
    "Removing ", length(exclude_cells), " cells from ",
    paste(clusters_to_remove, collapse = ", "), " before reclustering."
  )
}

exclusion_summary <- as.data.frame(table(
  removed_cluster = cell_metadata_before[[cluster_col]][rownames(cell_metadata_before) %in% exclude_cells],
  sample_id = cell_metadata_before$Sample[rownames(cell_metadata_before) %in% exclude_cells]
))
colnames(exclusion_summary) <- c("removed_cluster", "sample_id", "n_cells_removed")
write.csv(
  exclusion_summary,
  file = file.path(output_dir, "removed_clusters_before_reclustering.csv"),
  row.names = FALSE
)

keep_cells <- setdiff(rownames(cell_metadata_before), exclude_cells)
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
  count_cells_by_sample(cell_metadata_before, "before_cluster_exclusion"),
  count_cells_by_sample(
    cell_metadata_before[rownames(cell_metadata_before) %in% exclude_cells, , drop = FALSE],
    "removed_clusters_C1_C10"
  ),
  count_cells_by_sample(cell_metadata_after, "after_cluster_exclusion_reclustered")
)
cell_counts_summary <- rbind(
  cell_counts_summary,
  aggregate(n_cells ~ stage, data = cell_counts_summary, sum) |>
    transform(sample_id = "Total") |>
    subset(select = c("stage", "sample_id", "n_cells"))
)
write.csv(
  cell_counts_summary[order(cell_counts_summary$stage, cell_counts_summary$sample_id), ],
  file = file.path(output_dir, "cell_counts_after_cluster_exclusion_reclustered.csv"),
  row.names = FALSE
)

write_filter_settings(
  path = file.path(output_dir, "cluster_exclusion_reclustering_settings.csv"),
  settings = list(
    source_cluster_col = cluster_col,
    removed_clusters = paste(clusters_to_remove, collapse = ","),
    atac_lsi_recomputed = TRUE,
    rna_lsi_recomputed = TRUE,
    combined_lsi_recomputed = TRUE,
    clusters_recomputed = "Clusters,Clusters_RNA,Clusters_Combined",
    note = "Cells in the source Clusters_Combined labels C1 and C10 were removed after inspection, then ATAC/RNA/combined embeddings and clusters were recomputed. New cluster labels are renumbered by ArchR and may reuse C1 or C10 for different cells."
  )
)

saveArchRProject(ArchRProj = proj, outputDirectory = output_dir, load = FALSE)

message("Reclustered project after removing ", paste(clusters_to_remove, collapse = ", "), ".")
message("Remaining cells: ", nrow(cell_metadata_after))
