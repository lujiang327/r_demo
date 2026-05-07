#!/usr/bin/env Rscript

source(file.path("config", "analysis_config.R"))
source(file.path("R", "project_helpers.R"))

suppressPackageStartupMessages(library(ArchR))

set.seed(random_seed)
setup_archr_session()

sweep_dir <- file.path(figure_dir, "parameter_sweep")
ensure_dirs(sweep_dir, log_dir)

proj <- loadArchRProject(path = output_dir)

available_reduced_dims <- names(proj@reducedDims)
candidate_reduced_dims <- c("IterativeLSI", "LSI_RNA", "LSI_Combined")
reduced_dims_to_use <- intersect(candidate_reduced_dims, available_reduced_dims)

if (length(reduced_dims_to_use) == 0) {
  stop(
    "No supported reducedDims found. Expected one of: ",
    paste(candidate_reduced_dims, collapse = ", ")
  )
}

umap_grid <- expand.grid(
  reducedDims = reduced_dims_to_use,
  dimsToUse = c(20, 30),
  nNeighbors = c(15, 30, 50),
  minDist = c(0.2, 0.5, 0.8),
  stringsAsFactors = FALSE
)

cluster_grid <- expand.grid(
  reducedDims = reduced_dims_to_use,
  dimsToUse = c(20, 30),
  resolution = c(0.2, 0.4, 0.8, 1.2),
  stringsAsFactors = FALSE
)

message("Running UMAP sweep with ", nrow(umap_grid), " settings.")
for (i in seq_len(nrow(umap_grid))) {
  params <- umap_grid[i, ]
  embedding_name <- paste(
    "UMAP",
    params$reducedDims,
    paste0("d", params$dimsToUse),
    paste0("n", params$nNeighbors),
    paste0("md", gsub("[.]", "", params$minDist)),
    sep = "_"
  )

  proj <- addUMAP(
    ArchRProj = proj,
    reducedDims = params$reducedDims,
    dimsToUse = seq_len(params$dimsToUse),
    name = embedding_name,
    nNeighbors = params$nNeighbors,
    minDist = params$minDist,
    metric = "cosine",
    force = TRUE
  )

  p_sample <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "cellColData",
    name = "Sample",
    embedding = embedding_name
  )

  plotPDF(
    p_sample,
    name = paste0(embedding_name, "_sample.pdf"),
    ArchRProj = proj,
    addDOC = FALSE,
    width = 6,
    height = 6
  )
  file.copy(
    file.path(output_dir, "Plots", paste0(embedding_name, "_sample.pdf")),
    sweep_dir,
    overwrite = TRUE
  )
}

message("Running clustering sweep with ", nrow(cluster_grid), " settings.")
cluster_summaries <- list()
for (i in seq_len(nrow(cluster_grid))) {
  params <- cluster_grid[i, ]
  cluster_name <- paste(
    "Clusters",
    params$reducedDims,
    paste0("d", params$dimsToUse),
    paste0("res", gsub("[.]", "", params$resolution)),
    sep = "_"
  )

  proj <- addClusters(
    input = proj,
    reducedDims = params$reducedDims,
    dimsToUse = seq_len(params$dimsToUse),
    method = "Seurat",
    name = cluster_name,
    resolution = params$resolution,
    force = TRUE
  )

  cell_metadata <- as.data.frame(getCellColData(proj))
  composition <- as.data.frame(table(cell_metadata[[cluster_name]], cell_metadata$Sample))
  colnames(composition) <- c("cluster", "sample_id", "n_cells")
  composition$cluster_setting <- cluster_name
  composition$reducedDims <- params$reducedDims
  composition$dimsToUse <- params$dimsToUse
  composition$resolution <- params$resolution
  cluster_summaries[[cluster_name]] <- composition
}

cluster_summary <- do.call(rbind, cluster_summaries)
write.csv(
  cluster_summary,
  file = file.path(output_dir, "clustering_parameter_sweep_counts.csv"),
  row.names = FALSE
)

cell_metadata <- as.data.frame(getCellColData(proj))
write.csv(
  cell_metadata,
  file = file.path(output_dir, "cell_metadata_after_parameter_sweep.csv"),
  row.names = TRUE
)

write.csv(
  umap_grid,
  file = file.path(output_dir, "umap_parameter_sweep_settings.csv"),
  row.names = FALSE
)
write.csv(
  cluster_grid,
  file = file.path(output_dir, "clustering_parameter_sweep_settings.csv"),
  row.names = FALSE
)

saveArchRProject(ArchRProj = proj, outputDirectory = output_dir, load = FALSE)

message("Parameter sweep plots copied to: ", sweep_dir)
message("Cluster composition summary: ", file.path(output_dir, "clustering_parameter_sweep_counts.csv"))
