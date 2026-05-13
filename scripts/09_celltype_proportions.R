#!/usr/bin/env Rscript

source(file.path("config", "analysis_config.R"))
source(file.path("R", "project_helpers.R"))

suppressPackageStartupMessages({
  library(ggplot2)
})

proportion_dir <- file.path(figure_dir, "proportions")
ensure_dirs(output_dir, proportion_dir)

metadata_path <- file.path(output_dir, "cell_metadata_with_tentative_celltypes.csv")
if (!file.exists(metadata_path)) {
  stop("Missing cell type metadata: ", metadata_path, "\nRun script 07 first.")
}

cell_metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
cell_metadata$cell_id <- cell_metadata[[1]]

required_cols <- c("Sample", "Clusters_Combined")
missing_cols <- setdiff(required_cols, colnames(cell_metadata))
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
}

if ("tentative_celltype" %in% colnames(cell_metadata)) {
  cell_metadata$celltype_for_proportion <- cell_metadata$tentative_celltype
} else if ("TentativeCelltype" %in% colnames(cell_metadata)) {
  cell_metadata$celltype_for_proportion <- cell_metadata$TentativeCelltype
} else {
  stop("Missing tentative cell type column. Run script 07 first.")
}

if ("MGCandidate" %in% colnames(cell_metadata)) {
  cell_metadata$mg_candidate_for_proportion <- cell_metadata$MGCandidate
} else if ("mg_candidate" %in% colnames(cell_metadata)) {
  cell_metadata$mg_candidate_for_proportion <- ifelse(cell_metadata$mg_candidate, "MG_candidate", "Not_MG_candidate")
} else {
  cell_metadata$mg_candidate_for_proportion <- NA_character_
}

sample_order <- unique(cell_metadata$Sample)
if (all(c("TH1", "TH2") %in% sample_order)) {
  sample_order <- c("TH1", "TH2")
}
cell_metadata$Sample <- factor(cell_metadata$Sample, levels = sample_order)

count_proportions <- function(metadata, group_col, stage) {
  counts <- as.data.frame(
    table(metadata$Sample, metadata[[group_col]]),
    stringsAsFactors = FALSE
  )
  colnames(counts) <- c("sample_id", "group", "n_cells")
  counts <- counts[counts$n_cells > 0, , drop = FALSE]

  totals <- aggregate(n_cells ~ sample_id, data = counts, sum)
  colnames(totals)[2] <- "sample_total"
  counts <- merge(counts, totals, by = "sample_id", all.x = TRUE)
  counts$proportion <- counts$n_cells / counts$sample_total
  counts$percent <- counts$proportion * 100
  counts$stage <- stage
  counts[, c("stage", "sample_id", "group", "n_cells", "sample_total", "proportion", "percent")]
}

compare_two_samples <- function(prop_df, sample_a = "TH1", sample_b = "TH2") {
  groups <- sort(unique(prop_df$group))
  wide <- data.frame(group = groups, stringsAsFactors = FALSE)

  for (sample_id in c(sample_a, sample_b)) {
    sample_df <- prop_df[prop_df$sample_id == sample_id, c("group", "n_cells", "percent"), drop = FALSE]
    colnames(sample_df) <- c("group", paste0(sample_id, "_n_cells"), paste0(sample_id, "_percent"))
    wide <- merge(wide, sample_df, by = "group", all.x = TRUE)
  }

  numeric_cols <- setdiff(colnames(wide), "group")
  for (col_name in numeric_cols) {
    wide[[col_name]][is.na(wide[[col_name]])] <- 0
  }

  wide$delta_percent_TH2_minus_TH1 <- wide[[paste0(sample_b, "_percent")]] - wide[[paste0(sample_a, "_percent")]]
  wide$ratio_percent_TH2_over_TH1 <- ifelse(
    wide[[paste0(sample_a, "_percent")]] == 0,
    NA_real_,
    wide[[paste0(sample_b, "_percent")]] / wide[[paste0(sample_a, "_percent")]]
  )
  wide[order(abs(wide$delta_percent_TH2_minus_TH1), decreasing = TRUE), ]
}

plot_stacked_proportions <- function(prop_df, title, path_prefix) {
  prop_df$sample_id <- factor(prop_df$sample_id, levels = sample_order)
  prop_df$group <- factor(prop_df$group, levels = sort(unique(prop_df$group)))

  p <- ggplot(prop_df, aes(x = sample_id, y = percent, fill = group)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.2) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(x = NULL, y = "Cells (%)", fill = NULL, title = title) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )

  ggsave(paste0(path_prefix, ".pdf"), p, width = 7, height = 4.5)
  ggsave(paste0(path_prefix, ".png"), p, width = 7, height = 4.5, dpi = 180)
  invisible(p)
}

plot_delta <- function(comparison_df, title, path_prefix, max_groups = Inf) {
  plot_df <- comparison_df[order(abs(comparison_df$delta_percent_TH2_minus_TH1), decreasing = TRUE), , drop = FALSE]
  if (is.finite(max_groups) && nrow(plot_df) > max_groups) {
    plot_df <- plot_df[seq_len(max_groups), , drop = FALSE]
  }
  plot_df$group <- factor(plot_df$group, levels = rev(plot_df$group))

  p <- ggplot(plot_df, aes(x = group, y = delta_percent_TH2_minus_TH1, fill = delta_percent_TH2_minus_TH1 > 0)) +
    geom_col(width = 0.72) +
    coord_flip() +
    scale_fill_manual(values = c(`TRUE` = "#b2182b", `FALSE` = "#2166ac"), guide = "none") +
    labs(x = NULL, y = "Delta cells (% TH2 - % TH1)", title = title) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  height <- max(4, min(10, 0.28 * nrow(plot_df) + 1.8))
  ggsave(paste0(path_prefix, ".pdf"), p, width = 7, height = height)
  ggsave(paste0(path_prefix, ".png"), p, width = 7, height = height, dpi = 180)
  invisible(p)
}

plot_grouped_cluster_percent <- function(prop_df, title, path_prefix) {
  plot_df <- prop_df
  plot_df$sample_id <- factor(plot_df$sample_id, levels = sample_order)
  plot_df$group <- factor(plot_df$group, levels = sort(unique(plot_df$group)))

  p <- ggplot(plot_df, aes(x = group, y = percent, fill = sample_id)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.7) +
    scale_fill_manual(values = c(TH1 = "#4C78A8", TH2 = "#E45756"), drop = FALSE) +
    labs(x = "Combined cluster", y = "Cells (%)", fill = NULL, title = title) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )

  ggsave(paste0(path_prefix, ".pdf"), p, width = 9, height = 4.8)
  ggsave(paste0(path_prefix, ".png"), p, width = 9, height = 4.8, dpi = 180)
  invisible(p)
}

plot_cluster_ratio <- function(comparison_df, title, path_prefix) {
  plot_df <- comparison_df[order(comparison_df$group), , drop = FALSE]
  plot_df$group <- factor(plot_df$group, levels = plot_df$group)

  p <- ggplot(plot_df, aes(x = group, y = ratio_percent_TH2_over_TH1)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.35) +
    geom_col(width = 0.72, fill = "#7B3294") +
    labs(x = "Combined cluster", y = "Ratio of cell percentage (TH2 / TH1)", title = title) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(paste0(path_prefix, ".pdf"), p, width = 9, height = 4.8)
  ggsave(paste0(path_prefix, ".png"), p, width = 9, height = 4.8, dpi = 180)
  invisible(p)
}

celltype_props <- count_proportions(cell_metadata, "celltype_for_proportion", "tentative_celltype")
cluster_props <- count_proportions(cell_metadata, "Clusters_Combined", "clusters_combined")
mg_candidate_props <- count_proportions(cell_metadata, "mg_candidate_for_proportion", "mg_candidate")

celltype_comparison <- compare_two_samples(celltype_props)
cluster_comparison <- compare_two_samples(cluster_props)
mg_candidate_comparison <- compare_two_samples(mg_candidate_props)

write.csv(celltype_props, file.path(output_dir, "celltype_proportions_by_sample.csv"), row.names = FALSE)
write.csv(cluster_props, file.path(output_dir, "cluster_proportions_by_sample.csv"), row.names = FALSE)
write.csv(mg_candidate_props, file.path(output_dir, "mg_candidate_proportions_by_sample.csv"), row.names = FALSE)

write.csv(celltype_comparison, file.path(output_dir, "celltype_proportion_comparison_TH2_vs_TH1.csv"), row.names = FALSE)
write.csv(cluster_comparison, file.path(output_dir, "cluster_proportion_comparison_TH2_vs_TH1.csv"), row.names = FALSE)
write.csv(mg_candidate_comparison, file.path(output_dir, "mg_candidate_proportion_comparison_TH2_vs_TH1.csv"), row.names = FALSE)

plot_stacked_proportions(
  celltype_props,
  title = "Tentative cell type proportions",
  path_prefix = file.path(proportion_dir, "celltype_proportions_stacked")
)
plot_delta(
  celltype_comparison,
  title = "Tentative cell type proportion difference",
  path_prefix = file.path(proportion_dir, "celltype_proportion_delta")
)
plot_stacked_proportions(
  cluster_props,
  title = "Combined cluster proportions",
  path_prefix = file.path(proportion_dir, "cluster_proportions_stacked")
)
plot_delta(
  cluster_comparison,
  title = "Top combined cluster proportion differences",
  path_prefix = file.path(proportion_dir, "cluster_proportion_delta_top20"),
  max_groups = 20
)
plot_grouped_cluster_percent(
  cluster_props,
  title = "Combined cluster percentages by sample",
  path_prefix = file.path(proportion_dir, "cluster_percent_by_sample_grouped")
)
plot_cluster_ratio(
  cluster_comparison,
  title = "Combined cluster percentage ratio",
  path_prefix = file.path(proportion_dir, "cluster_percent_ratio_TH2_over_TH1")
)
plot_stacked_proportions(
  mg_candidate_props,
  title = "MG-candidate proportions",
  path_prefix = file.path(proportion_dir, "mg_candidate_proportions_stacked")
)
plot_delta(
  mg_candidate_comparison,
  title = "MG-candidate proportion difference",
  path_prefix = file.path(proportion_dir, "mg_candidate_proportion_delta")
)

message("Wrote proportion tables to: ", output_dir)
message("Wrote proportion plots to: ", proportion_dir)
