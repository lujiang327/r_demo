#!/usr/bin/env Rscript

source(file.path("config", "analysis_config.R"))
source(file.path("R", "project_helpers.R"))

suppressPackageStartupMessages({
  library(ArchR)
  library(ggplot2)
})

set.seed(random_seed)
setup_archr_session()

review_dir <- file.path("results", "combined_umap_parameter_review")
ensure_dirs(review_dir)

cluster_col <- "Clusters_Combined"
reduced_dims <- "LSI_Combined"
clusters_of_interest <- c("C18", "C19", "C20")

umap_configs <- list(
  list(id = "dims01_30_nn30_md04", dims = 1:30, n_neighbors = 30, min_dist = 0.4),
  list(id = "dims02_30_nn30_md04", dims = 2:30, n_neighbors = 30, min_dist = 0.4),
  list(id = "dims01_20_nn30_md04", dims = 1:20, n_neighbors = 30, min_dist = 0.4),
  list(id = "dims01_40_nn30_md04", dims = 1:40, n_neighbors = 30, min_dist = 0.4),
  list(id = "dims01_30_nn15_md04", dims = 1:30, n_neighbors = 15, min_dist = 0.4),
  list(id = "dims01_30_nn50_md04", dims = 1:30, n_neighbors = 50, min_dist = 0.4),
  list(id = "dims01_30_nn30_md02", dims = 1:30, n_neighbors = 30, min_dist = 0.2),
  list(id = "dims01_30_nn30_md08", dims = 1:30, n_neighbors = 30, min_dist = 0.8)
)

plot_embedding_png <- function(proj, embedding_name, color_col, file_name, title = NULL) {
  p <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "cellColData",
    name = color_col,
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

summarize_embedding_centroids <- function(proj, embedding_name, config_id) {
  embedding <- getEmbedding(proj, embedding = embedding_name, returnDF = TRUE)
  if (ncol(embedding) != 2) {
    stop("Embedding `", embedding_name, "` does not have exactly 2 dimensions.")
  }
  colnames(embedding) <- c("UMAP_1", "UMAP_2")

  metadata <- as.data.frame(getCellColData(proj))
  metadata <- metadata[rownames(embedding), , drop = FALSE]

  df <- data.frame(
    cell_id = rownames(embedding),
    cluster = metadata[[cluster_col]],
    sample_id = metadata$Sample,
    UMAP_1 = embedding$UMAP_1,
    UMAP_2 = embedding$UMAP_2,
    stringsAsFactors = FALSE
  )

  centroids <- aggregate(
    cbind(UMAP_1, UMAP_2) ~ cluster,
    data = df,
    FUN = mean
  )
  counts <- as.data.frame(table(df$cluster), stringsAsFactors = FALSE)
  colnames(counts) <- c("cluster", "n_cells")
  centroids <- merge(centroids, counts, by = "cluster", all.x = TRUE)
  centroids$config_id <- config_id
  centroids[, c("config_id", "cluster", "n_cells", "UMAP_1", "UMAP_2")]
}

proj <- loadArchRProject(path = output_dir)
metadata <- as.data.frame(getCellColData(proj))

if (!cluster_col %in% colnames(metadata)) {
  stop("Missing `", cluster_col, "`. Run final reclustering first.")
}
if (!reduced_dims %in% names(proj@reducedDims)) {
  stop("Missing reducedDims `", reduced_dims, "`.")
}

centroid_rows <- list()
settings_rows <- list()

for (config in umap_configs) {
  embedding_name <- paste0("UMAP_Combined_", config$id)
  message("Adding ", embedding_name)

  proj <- addUMAP(
    ArchRProj = proj,
    reducedDims = reduced_dims,
    name = embedding_name,
    nNeighbors = config$n_neighbors,
    minDist = config$min_dist,
    metric = "cosine",
    dimsToUse = config$dims,
    force = TRUE
  )

  plot_embedding_png(
    proj,
    embedding_name = embedding_name,
    color_col = cluster_col,
    file_name = paste0(embedding_name, "_clusters.png"),
    title = paste0("Combined UMAP ", config$id, " - clusters")
  )
  plot_embedding_png(
    proj,
    embedding_name = embedding_name,
    color_col = "Sample",
    file_name = paste0(embedding_name, "_sample.png"),
    title = paste0("Combined UMAP ", config$id, " - sample")
  )

  centroids <- summarize_embedding_centroids(proj, embedding_name, config$id)
  centroid_rows[[config$id]] <- centroids

  settings_rows[[config$id]] <- data.frame(
    config_id = config$id,
    reduced_dims = reduced_dims,
    dims_to_use = paste(range(config$dims), collapse = ":"),
    n_neighbors = config$n_neighbors,
    min_dist = config$min_dist,
    metric = "cosine",
    stringsAsFactors = FALSE
  )
}

centroid_df <- do.call(rbind, centroid_rows)
write.csv(
  centroid_df,
  file = file.path(review_dir, "combined_umap_parameter_cluster_centroids.csv"),
  row.names = FALSE
)

interest_centroids <- centroid_df[centroid_df$cluster %in% clusters_of_interest, , drop = FALSE]
write.csv(
  interest_centroids,
  file = file.path(review_dir, "clusters_C18_C19_C20_umap_centroids.csv"),
  row.names = FALSE
)

write.csv(
  do.call(rbind, settings_rows),
  file = file.path(review_dir, "combined_umap_parameter_review_settings.csv"),
  row.names = FALSE
)

message("Combined UMAP parameter review outputs written to: ", review_dir)
