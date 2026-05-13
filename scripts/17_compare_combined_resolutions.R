#!/usr/bin/env Rscript

source(file.path("config", "analysis_config.R"))
source(file.path("R", "project_helpers.R"))

suppressPackageStartupMessages({
  library(ArchR)
  library(ggplot2)
})

set.seed(random_seed)
setup_archr_session()

review_dir <- file.path("results", "combined_resolution_review")
ensure_dirs(review_dir)

resolutions_to_try <- c(1.2, 1.4)
current_cluster_col <- "Clusters_Combined"
reduced_dims <- "LSI_Combined"
embedding_name <- "UMAP_Combined"
clusters_of_interest <- c("C18", "C19", "C20")

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

if (!current_cluster_col %in% colnames(metadata)) {
  stop("Missing `", current_cluster_col, "`. Run reclustering before this script.")
}
if (!reduced_dims %in% names(proj@reducedDims)) {
  stop("Missing reducedDims `", reduced_dims, "`.")
}

summary_rows <- list()

for (resolution in resolutions_to_try) {
  label <- resolution_label(resolution)
  cluster_col <- paste0("Clusters_Combined_r", label)

  message("Adding combined clusters at resolution ", resolution, ": ", cluster_col)
  proj <- addClusters(
    input = proj,
    reducedDims = reduced_dims,
    method = "Seurat",
    name = cluster_col,
    resolution = resolution,
    force = TRUE
  )

  metadata <- as.data.frame(getCellColData(proj))

  write_count_table(
    metadata,
    cluster_col = cluster_col,
    path = file.path(review_dir, paste0(cluster_col, "_counts_by_sample.csv"))
  )

  crosstab <- as.data.frame(table(
    current_cluster = metadata[[current_cluster_col]],
    review_cluster = metadata[[cluster_col]]
  ), stringsAsFactors = FALSE)
  colnames(crosstab) <- c("current_cluster", "review_cluster", "n_cells")
  crosstab <- crosstab[crosstab$n_cells > 0, , drop = FALSE]
  crosstab$resolution <- resolution
  write.csv(
    crosstab,
    file = file.path(review_dir, paste0(cluster_col, "_current_cluster_crosstab.csv")),
    row.names = FALSE
  )

  interest_split <- crosstab[crosstab$current_cluster %in% clusters_of_interest, , drop = FALSE]
  if (nrow(interest_split) > 0) {
    summary_rows[[cluster_col]] <- interest_split
  }

  plot_embedding_png(
    proj,
    cluster_col = cluster_col,
    file_name = paste0(cluster_col, "_umap.png"),
    title = paste0("Combined UMAP - resolution ", resolution)
  )
}

if (length(summary_rows) > 0) {
  write.csv(
    do.call(rbind, summary_rows),
    file = file.path(review_dir, "clusters_C18_C19_C20_split_summary.csv"),
    row.names = FALSE
  )
}

write_filter_settings(
  path = file.path(review_dir, "combined_resolution_review_settings.csv"),
  settings = list(
    current_cluster_col = current_cluster_col,
    reduced_dims = reduced_dims,
    embedding = embedding_name,
    resolutions_to_try = paste(resolutions_to_try, collapse = ","),
    clusters_of_interest = paste(clusters_of_interest, collapse = ","),
    note = "Review-only combined clustering. This script does not save the ArchRProject or overwrite final clusters."
  )
)

message("Combined resolution review outputs written to: ", review_dir)
