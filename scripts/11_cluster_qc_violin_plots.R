#!/usr/bin/env Rscript

source(file.path("config", "analysis_config.R"))
source(file.path("R", "project_helpers.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

qc_cluster_dir <- file.path(figure_dir, "cluster_qc")
ensure_dirs(output_dir, qc_cluster_dir)

metadata_candidates <- c(
  file.path(output_dir, "cell_metadata_after_rna_combined.csv"),
  file.path(output_dir, "cell_metadata_with_tentative_celltypes.csv"),
  file.path(output_dir, "cell_metadata_after_clustering.csv")
)
metadata_path <- metadata_candidates[file.exists(metadata_candidates)][1]
if (is.na(metadata_path)) {
  stop("No cell metadata CSV found. Run scripts 03/04 before this script.")
}

cell_col_data <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
cell_col_data$cell_id <- cell_col_data[[1]]

if ("nFrags" %in% colnames(cell_col_data)) {
  cell_col_data$log10_nFrags <- log10(cell_col_data$nFrags + 1)
}
if ("Gex_nUMI" %in% colnames(cell_col_data)) {
  cell_col_data$Gex_Log10_nUMI <- log10(cell_col_data$Gex_nUMI + 1)
}
if ("Gex_nGenes" %in% colnames(cell_col_data)) {
  cell_col_data$Gex_Log10_nGenes <- log10(cell_col_data$Gex_nGenes + 1)
}

if ("Gex_MitoRatio" %in% colnames(cell_col_data)) {
  mito <- cell_col_data$Gex_MitoRatio
  if (max(mito, na.rm = TRUE) <= 1) {
    mito <- mito * 100
  }
  cell_col_data$Gex_MitoPercent <- mito
}
if ("Gex_RiboRatio" %in% colnames(cell_col_data)) {
  ribo <- cell_col_data$Gex_RiboRatio
  if (max(ribo, na.rm = TRUE) <= 1) {
    ribo <- ribo * 100
  }
  cell_col_data$Gex_RiboPercent <- ribo
}

qc_metrics_atac <- c(
  "TSSEnrichment",
  "nFrags",
  "log10_nFrags",
  "NucleosomeRatio",
  "ReadsInTSS",
  "ReadsInPromoter",
  "PromoterRatio"
)
qc_metrics_atac <- qc_metrics_atac[qc_metrics_atac %in% colnames(cell_col_data)]

qc_metrics_rna <- c(
  "Gex_nUMI",
  "Gex_nGenes",
  "Gex_Log10_nUMI",
  "Gex_Log10_nGenes",
  "Gex_MitoPercent",
  "Gex_RiboPercent"
)
qc_metrics_rna <- qc_metrics_rna[qc_metrics_rna %in% colnames(cell_col_data)]

cluster_sets <- list(
  Combined = list(cluster_col = "Clusters_Combined", metrics = c(qc_metrics_atac, qc_metrics_rna)),
  RNA = list(cluster_col = "Clusters_RNA", metrics = qc_metrics_rna),
  ATAC = list(cluster_col = "Clusters", metrics = qc_metrics_atac)
)

sort_cluster_levels <- function(x) {
  x <- as.character(x)
  numeric_part <- suppressWarnings(as.integer(sub("^C", "", x)))
  if (all(!is.na(numeric_part))) {
    return(x[order(numeric_part)])
  }
  sort(unique(x))
}

make_qc_violin <- function(metadata, cluster_col, qc_metric, metric_type) {
  df <- data.frame(
    Cluster = metadata[[cluster_col]],
    Value = metadata[[qc_metric]],
    stringsAsFactors = FALSE
  )
  df <- df[complete.cases(df), , drop = FALSE]
  if (nrow(df) == 0) {
    return(NULL)
  }

  df$Cluster <- factor(df$Cluster, levels = sort_cluster_levels(unique(df$Cluster)))

  ggplot(df, aes(x = Cluster, y = Value, fill = Cluster)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, linewidth = 0.25) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 11)
    ) +
    labs(
      title = paste0(qc_metric, " (", metric_type, ")"),
      y = qc_metric,
      x = "Cluster"
    )
}

write_cluster_qc_summary <- function(metadata, cluster_col, metrics, out_path) {
  summary_rows <- lapply(metrics, function(metric) {
    df <- metadata[, c(cluster_col, metric), drop = FALSE]
    colnames(df) <- c("cluster", "value")
    df <- df[complete.cases(df), , drop = FALSE]
    if (nrow(df) == 0) {
      return(data.frame())
    }
    by_cluster <- split(df$value, df$cluster)
    do.call(rbind, lapply(names(by_cluster), function(cluster) {
      values <- by_cluster[[cluster]]
      data.frame(
        cluster = cluster,
        metric = metric,
        n_cells = length(values),
        mean = mean(values),
        median = stats::median(values),
        sd = stats::sd(values),
        q05 = as.numeric(stats::quantile(values, 0.05)),
        q25 = as.numeric(stats::quantile(values, 0.25)),
        q75 = as.numeric(stats::quantile(values, 0.75)),
        q95 = as.numeric(stats::quantile(values, 0.95)),
        stringsAsFactors = FALSE
      )
    }))
  })
  summary_df <- do.call(rbind, summary_rows)
  write.csv(summary_df, out_path, row.names = FALSE)
  invisible(summary_df)
}

for (cluster_set_name in names(cluster_sets)) {
  cluster_col <- cluster_sets[[cluster_set_name]]$cluster_col
  metrics <- unique(cluster_sets[[cluster_set_name]]$metrics)

  if (!cluster_col %in% colnames(cell_col_data) || length(metrics) == 0) {
    next
  }

  qc_plot_list <- list()
  for (qc_metric in metrics) {
    metric_type <- if (qc_metric %in% qc_metrics_atac) "ATAC" else "RNA"
    p <- make_qc_violin(cell_col_data, cluster_col, qc_metric, metric_type)
    if (!is.null(p)) {
      qc_plot_list[[qc_metric]] <- p
    }
  }

  if (length(qc_plot_list) == 0) {
    next
  }

  n_cols <- min(2, length(qc_plot_list))
  n_rows <- ceiling(length(qc_plot_list) / n_cols)
  combined_qc_plot <- wrap_plots(qc_plot_list, ncol = n_cols, nrow = n_rows)

  base_name <- file.path(qc_cluster_dir, paste0("Clusters_", cluster_set_name, "_QC"))
  ggsave(
    filename = paste0(base_name, ".pdf"),
    plot = combined_qc_plot,
    width = 8 * n_cols,
    height = 5.5 * n_rows,
    limitsize = FALSE
  )
  ggsave(
    filename = paste0(base_name, ".png"),
    plot = combined_qc_plot,
    width = 8 * n_cols,
    height = 5.5 * n_rows,
    dpi = 160,
    limitsize = FALSE
  )

  write_cluster_qc_summary(
    cell_col_data,
    cluster_col = cluster_col,
    metrics = metrics,
    out_path = file.path(output_dir, paste0("cluster_qc_summary_", cluster_set_name, ".csv"))
  )

  message("Wrote ", cluster_set_name, " cluster QC plots: ", base_name, ".pdf/.png")
}

message("Cluster QC outputs written to: ", qc_cluster_dir)
