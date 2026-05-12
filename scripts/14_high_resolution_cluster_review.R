#!/usr/bin/env Rscript

source(file.path("config", "analysis_config.R"))
source(file.path("R", "project_helpers.R"))

suppressPackageStartupMessages({
  library(ArchR)
  library(ggplot2)
})

set.seed(random_seed)
setup_archr_session()

review_dir <- file.path("results", "high_resolution_cluster_review")
ensure_dirs(review_dir)

target_cluster <- "C20"
source_cluster_col <- "Clusters_Combined"
reduced_dims <- "LSI_Combined"
embedding_name <- "UMAP_Combined"

global_resolutions <- c(1.2, 1.4, 1.6)
c20_resolutions <- c(0.8, 1.2, 1.6, 2.0, 2.5)

resolution_label <- function(x) {
  gsub("\\.", "p", as.character(x))
}

sort_cluster_levels <- function(x) {
  x <- as.character(x)
  numeric_part <- suppressWarnings(as.integer(sub("^C", "", x)))
  if (all(!is.na(numeric_part))) {
    return(x[order(numeric_part)])
  }
  sort(unique(x))
}

write_count_table <- function(metadata, cluster_col, path) {
  counts <- as.data.frame(table(metadata[[cluster_col]], metadata$Sample), stringsAsFactors = FALSE)
  colnames(counts) <- c("cluster", "sample_id", "n_cells")
  counts$cluster <- factor(counts$cluster, levels = sort_cluster_levels(unique(counts$cluster)))
  counts <- counts[order(counts$cluster, counts$sample_id), , drop = FALSE]
  counts$cluster <- as.character(counts$cluster)
  write.csv(counts, path, row.names = FALSE)
  invisible(counts)
}

plot_embedding_png <- function(proj, cluster_col, file_name, title = NULL) {
  p <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "cellColData",
    name = cluster_col,
    embedding = embedding_name
  )
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  ggsave(
    filename = file.path(review_dir, file_name),
    plot = p,
    width = 6,
    height = 5.5,
    dpi = 180
  )
}

proj <- loadArchRProject(path = output_dir)
metadata <- as.data.frame(getCellColData(proj))

if (!source_cluster_col %in% colnames(metadata)) {
  stop("Missing `", source_cluster_col, "`. Run scripts 04 and 12 first.")
}
if (!target_cluster %in% metadata[[source_cluster_col]]) {
  stop("Target cluster `", target_cluster, "` was not found in `", source_cluster_col, "`.")
}

highres_summary <- list()

for (resolution in global_resolutions) {
  label <- resolution_label(resolution)
  cluster_name <- paste0("Clusters_Combined_highres_r", label)
  message("Adding global high-resolution clusters: ", cluster_name)

  proj <- addClusters(
    input = proj,
    reducedDims = reduced_dims,
    method = "Seurat",
    name = cluster_name,
    resolution = resolution,
    force = TRUE
  )

  metadata <- as.data.frame(getCellColData(proj))
  write_count_table(
    metadata,
    cluster_col = cluster_name,
    path = file.path(review_dir, paste0(cluster_name, "_counts_by_sample.csv"))
  )

  crosstab <- as.data.frame(table(
    source_cluster = metadata[[source_cluster_col]],
    highres_cluster = metadata[[cluster_name]]
  ), stringsAsFactors = FALSE)
  colnames(crosstab) <- c("source_cluster", "highres_cluster", "n_cells")
  crosstab <- crosstab[crosstab$n_cells > 0, , drop = FALSE]
  write.csv(
    crosstab,
    file = file.path(review_dir, paste0(cluster_name, "_source_crosstab.csv")),
    row.names = FALSE
  )

  c20_split <- crosstab[crosstab$source_cluster == target_cluster, , drop = FALSE]
  c20_split$resolution <- resolution
  c20_split$cluster_mode <- "global_highres"
  highres_summary[[cluster_name]] <- c20_split

  plot_embedding_png(
    proj,
    cluster_col = cluster_name,
    file_name = paste0(cluster_name, "_umap.png"),
    title = paste0("Combined UMAP - global high-res r", resolution)
  )
}

c20_cells <- rownames(metadata)[metadata[[source_cluster_col]] == target_cluster]
c20_proj <- proj[c20_cells]
c20_metadata <- as.data.frame(getCellColData(c20_proj))

message("C20 contains ", nrow(c20_metadata), " cells before C20-only subclustering.")

c20_assignments <- data.frame(
  cell_id = rownames(c20_metadata),
  sample_id = c20_metadata$Sample,
  source_cluster = c20_metadata[[source_cluster_col]],
  stringsAsFactors = FALSE
)

for (resolution in c20_resolutions) {
  label <- resolution_label(resolution)
  cluster_name <- paste0("C20_subclusters_r", label)
  message("Adding C20-only subclusters: ", cluster_name)

  c20_proj <- addClusters(
    input = c20_proj,
    reducedDims = reduced_dims,
    method = "Seurat",
    name = cluster_name,
    resolution = resolution,
    force = TRUE
  )

  c20_metadata <- as.data.frame(getCellColData(c20_proj))
  c20_assignments[[cluster_name]] <- c20_metadata[[cluster_name]][match(c20_assignments$cell_id, rownames(c20_metadata))]

  write_count_table(
    c20_metadata,
    cluster_col = cluster_name,
    path = file.path(review_dir, paste0(cluster_name, "_counts_by_sample.csv"))
  )

  c20_counts <- as.data.frame(table(c20_metadata[[cluster_name]]), stringsAsFactors = FALSE)
  colnames(c20_counts) <- c("highres_cluster", "n_cells")
  c20_counts$resolution <- resolution
  c20_counts$source_cluster <- target_cluster
  c20_counts$cluster_mode <- "c20_only"
  highres_summary[[cluster_name]] <- c20_counts[, c("source_cluster", "highres_cluster", "n_cells", "resolution", "cluster_mode")]

  plot_embedding_png(
    c20_proj,
    cluster_col = cluster_name,
    file_name = paste0(cluster_name, "_umap.png"),
    title = paste0("C20-only UMAP - subclusters r", resolution)
  )
}

write.csv(
  c20_assignments,
  file = file.path(review_dir, "c20_high_resolution_cluster_assignments.csv"),
  row.names = FALSE
)

summary_df <- do.call(rbind, highres_summary)
write.csv(
  summary_df,
  file = file.path(review_dir, "c20_high_resolution_split_summary.csv"),
  row.names = FALSE
)

write_filter_settings(
  path = file.path(review_dir, "high_resolution_cluster_review_settings.csv"),
  settings = list(
    source_cluster_col = source_cluster_col,
    target_cluster = target_cluster,
    reduced_dims = reduced_dims,
    embedding = embedding_name,
    global_resolutions = paste(global_resolutions, collapse = ","),
    c20_resolutions = paste(c20_resolutions, collapse = ","),
    note = "Review-only high-resolution clustering. This script does not save the ArchRProject or remove cells."
  )
)

message("High-resolution cluster review outputs written to: ", review_dir)
