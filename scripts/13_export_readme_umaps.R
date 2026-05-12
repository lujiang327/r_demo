#!/usr/bin/env Rscript

source(file.path("config", "analysis_config.R"))
source(file.path("R", "project_helpers.R"))

suppressPackageStartupMessages({
  library(ArchR)
  library(ggplot2)
})

set.seed(random_seed)
setup_archr_session()

readme_figure_dir <- file.path("results", "readme_figures")
ensure_dirs(readme_figure_dir)

proj <- loadArchRProject(path = output_dir)

save_embedding_png <- function(embedding, color_by, name, file_name) {
  p <- plotEmbedding(
    ArchRProj = proj,
    colorBy = color_by,
    name = name,
    embedding = embedding
  )

  ggsave(
    filename = file.path(readme_figure_dir, file_name),
    plot = p,
    width = 5.5,
    height = 5,
    dpi = 180
  )
}

save_embedding_png("UMAP", "cellColData", "Sample", "atac_umap_sample.png")
save_embedding_png("UMAP", "cellColData", "Clusters", "atac_umap_clusters.png")

save_embedding_png("UMAP_RNA", "cellColData", "Sample", "rna_umap_sample.png")
save_embedding_png("UMAP_Combined", "cellColData", "Sample", "combined_umap_sample.png")
save_embedding_png("UMAP_Combined", "cellColData", "Clusters_Combined", "combined_umap_clusters.png")

message("Wrote README UMAP PNGs to: ", readme_figure_dir)
