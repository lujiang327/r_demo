#!/usr/bin/env Rscript

source(file.path("config", "analysis_config.R"))
source(file.path("R", "project_helpers.R"))

suppressPackageStartupMessages({
  library(ArchR)
  library(ggplot2)
})

set.seed(random_seed)
setup_archr_session()

label_free_dir <- file.path(figure_dir, "label_free_umaps")
readme_dir <- file.path("results", "readme_figures")
ensure_dirs(label_free_dir, readme_dir)

proj <- loadArchRProject(path = output_dir)
metadata <- as.data.frame(getCellColData(proj))

sort_cluster_levels <- function(x) {
  x <- as.character(unique(x))
  numeric_part <- suppressWarnings(as.integer(sub("^C", "", x)))
  x[order(is.na(numeric_part), numeric_part, x)]
}

manual_celltype_colors <- c(
  AC = "#D62728",
  Cone = "#2B2C7E",
  HC = "#1B8E3E",
  nPRC = "#8A2C91",
  `Pre-AC` = "#F07F2F",
  `Pre-Rod` = "#FFD900",
  pRPC = "#8FA8E3",
  RGC = "#C95BAE",
  Rod = "#D8A45F",
  `Starburst AC` = "#81D4E5"
)

get_umap_df <- function(proj, embedding_name, metadata) {
  embedding <- getEmbedding(proj, embedding = embedding_name, returnDF = TRUE)
  if (ncol(embedding) != 2) {
    stop("Embedding `", embedding_name, "` does not have exactly 2 columns.")
  }
  colnames(embedding) <- c("UMAP_1", "UMAP_2")
  common_cells <- intersect(rownames(embedding), rownames(metadata))
  cbind(
    data.frame(cell_id = common_cells, embedding[common_cells, , drop = FALSE], check.names = FALSE),
    metadata[common_cells, , drop = FALSE]
  )
}

plot_label_free_umap <- function(df, color_col, title, path_prefix, cluster_like = FALSE, color_values = NULL) {
  plot_df <- df[!is.na(df[[color_col]]), , drop = FALSE]
  plot_df[[color_col]] <- as.character(plot_df[[color_col]])
  if (cluster_like) {
    plot_df[[color_col]] <- factor(plot_df[[color_col]], levels = sort_cluster_levels(plot_df[[color_col]]))
  } else if (!is.null(color_values)) {
    plot_df[[color_col]] <- factor(plot_df[[color_col]], levels = names(color_values)[names(color_values) %in% unique(plot_df[[color_col]])])
  }

  p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = .data[[color_col]])) +
    geom_point(size = 0.12, alpha = 0.9) +
    guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1), ncol = 2)) +
    labs(x = "UMAP_1", y = "UMAP_2", color = NULL, title = title) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(linewidth = 0.25, color = "black"),
      legend.key.height = unit(0.42, "lines"),
      legend.text = element_text(size = 8)
    )
  if (!is.null(color_values)) {
    p <- p + scale_color_manual(values = color_values, drop = FALSE)
  }

  ggsave(paste0(path_prefix, ".pdf"), p, width = 6.2, height = 5.6)
  ggsave(paste0(path_prefix, ".png"), p, width = 6.2, height = 5.6, dpi = 220)
  invisible(p)
}

combined_df <- get_umap_df(proj, "UMAP_Combined", metadata)
plot_label_free_umap(
  combined_df,
  color_col = "Clusters_Combined",
  title = "Combined UMAP by cluster",
  path_prefix = file.path(label_free_dir, "combined_umap_clusters_no_labels"),
  cluster_like = TRUE
)
plot_label_free_umap(
  combined_df,
  color_col = "Sample",
  title = "Combined UMAP by sample",
  path_prefix = file.path(label_free_dir, "combined_umap_sample_no_labels")
)

file.copy(
  file.path(label_free_dir, "combined_umap_clusters_no_labels.png"),
  file.path(readme_dir, "combined_umap_clusters_no_labels.png"),
  overwrite = TRUE
)
file.copy(
  file.path(label_free_dir, "combined_umap_sample_no_labels.png"),
  file.path(readme_dir, "combined_umap_sample_no_labels.png"),
  overwrite = TRUE
)

if ("ManualCelltype" %in% colnames(combined_df)) {
  plot_label_free_umap(
    combined_df,
    color_col = "ManualCelltype",
    title = "Combined UMAP by manual cell type",
    path_prefix = file.path(label_free_dir, "combined_umap_manual_celltype_no_labels"),
    color_values = manual_celltype_colors
  )
  file.copy(
    file.path(label_free_dir, "combined_umap_manual_celltype_no_labels.png"),
    file.path(readme_dir, "combined_umap_manual_celltype_no_labels.png"),
    overwrite = TRUE
  )
}

message("Label-free UMAPs written to: ", label_free_dir)
