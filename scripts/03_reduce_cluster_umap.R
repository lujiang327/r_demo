#!/usr/bin/env Rscript

source(file.path("config", "analysis_config.R"))
source(file.path("R", "project_helpers.R"))

suppressPackageStartupMessages(library(ArchR))

set.seed(random_seed)
setup_archr_session()

ensure_dirs(file.path(figure_dir, "umap"), log_dir)

proj <- loadArchRProject(path = output_dir)

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
  dimsToUse = 1:30
)

proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)

proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine"
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

plotPDF(p_umap_sample, p_umap_cluster, name = "umap_sample_clusters.pdf", ArchRProj = proj, addDOC = FALSE, width = 6, height = 6)
file.copy(file.path(output_dir, "Plots", "umap_sample_clusters.pdf"), file.path(figure_dir, "umap"), overwrite = TRUE)

cell_col_data <- as.data.frame(getCellColData(proj))
write.csv(cell_col_data, file = file.path(output_dir, "cell_metadata_after_clustering.csv"), row.names = TRUE)

saveArchRProject(ArchRProj = proj, outputDirectory = output_dir, load = FALSE)

message("Reduced dimensions, clusters, and UMAP saved in: ", output_dir)
