#!/usr/bin/env Rscript

source(file.path("config", "analysis_config.R"))
source(file.path("R", "project_helpers.R"))

suppressPackageStartupMessages({
  library(ArchR)
})

set.seed(random_seed)
setup_archr_session()
ensure_dirs(file.path(figure_dir, "combined_umap"))

proj <- loadArchRProject(path = output_dir)

metadata <- as.data.frame(getCellColData(proj))
if (!"Clusters_Combined" %in% colnames(metadata)) {
  stop("Missing `Clusters_Combined`. Run combined reclustering before this script.")
}

proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "LSI_Combined",
  name = "UMAP_Combined",
  dimsToUse = combined_umap_dims_to_use,
  nNeighbors = combined_umap_n_neighbors,
  minDist = combined_umap_min_dist,
  metric = combined_umap_metric,
  force = TRUE
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
  p_combined_sample,
  p_combined_cluster,
  name = "primary_combined_umap.pdf",
  ArchRProj = proj,
  addDOC = FALSE,
  width = 6,
  height = 6
)
file.copy(
  file.path(output_dir, "Plots", "primary_combined_umap.pdf"),
  file.path(figure_dir, "combined_umap"),
  overwrite = TRUE
)

write_filter_settings(
  path = file.path(output_dir, "primary_combined_umap_settings.csv"),
  settings = list(
    reduced_dims = "LSI_Combined",
    embedding = "UMAP_Combined",
    dims_to_use = paste(range(combined_umap_dims_to_use), collapse = ":"),
    n_neighbors = combined_umap_n_neighbors,
    min_dist = combined_umap_min_dist,
    metric = combined_umap_metric,
    note = "Primary combined UMAP selected after parameter review."
  )
)

saveArchRProject(ArchRProj = proj, outputDirectory = output_dir, load = FALSE)

message("Primary combined UMAP updated and saved in: ", output_dir)
message("Settings written to: ", file.path(output_dir, "primary_combined_umap_settings.csv"))
