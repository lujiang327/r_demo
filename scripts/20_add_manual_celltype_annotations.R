#!/usr/bin/env Rscript

source(file.path("config", "analysis_config.R"))
source(file.path("R", "project_helpers.R"))

suppressPackageStartupMessages({
  library(ArchR)
  library(ggplot2)
})

set.seed(random_seed)
setup_archr_session()

manual_dir <- file.path(figure_dir, "manual_celltype_annotation")
proportion_dir <- file.path(figure_dir, "proportions")
ensure_dirs(output_dir, manual_dir, proportion_dir)

cluster_col <- "Clusters_Combined"
manual_col <- "ManualCelltype"

manual_annotation <- data.frame(
  cluster = paste0("C", c(12, 13, 14, 15, 16, 18, 19, 20, 2, 3, 1, 21, 10, 11, 8, 17, 4, 5, 6, 7, 9)),
  manual_celltype = c(
    rep("pRPC", 5),
    rep("nPRC", 3),
    "Cone",
    "Rod",
    rep("Pre-Rod", 2),
    rep("RGC", 2),
    "HC",
    "Pre-AC",
    rep("AC", 4),
    "Starburst AC"
  ),
  manual_celltype_full = c(
    rep("pRPC (primary progenitors)", 5),
    rep("nPRC (neurogenic progenitors)", 3),
    "Cone",
    "Rod",
    rep("Pre-Rod (rod precursors)", 2),
    rep("RGC", 2),
    "HC",
    "Pre-AC (amacrine precursors)",
    rep("AC", 4),
    "Starburst AC"
  ),
  stringsAsFactors = FALSE
)

manual_celltype_order <- c(
  "pRPC",
  "nPRC",
  "Pre-Rod",
  "Rod",
  "Cone",
  "RGC",
  "HC",
  "Pre-AC",
  "AC",
  "Starburst AC",
  "Unannotated"
)

sample_colors <- c(TH1 = "#4C78A8", TH2 = "#E45756")

sort_manual_levels <- function(x) {
  c(manual_celltype_order[manual_celltype_order %in% x], setdiff(sort(unique(x)), manual_celltype_order))
}

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
  groups <- sort_manual_levels(unique(prop_df$group))
  wide <- data.frame(group = groups, stringsAsFactors = FALSE)

  for (sample_id in c(sample_a, sample_b)) {
    sample_df <- prop_df[prop_df$sample_id == sample_id, c("group", "n_cells", "percent"), drop = FALSE]
    colnames(sample_df) <- c("group", paste0(sample_id, "_n_cells"), paste0(sample_id, "_percent"))
    wide <- merge(wide, sample_df, by = "group", all.x = TRUE, sort = FALSE)
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
  wide
}

plot_stacked_proportions <- function(prop_df, title, path_prefix) {
  plot_df <- prop_df
  plot_df$sample_id <- factor(plot_df$sample_id, levels = sample_order)
  plot_df$group <- factor(plot_df$group, levels = sort_manual_levels(plot_df$group))

  p <- ggplot(plot_df, aes(x = sample_id, y = percent, fill = group)) +
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

plot_grouped_percent <- function(prop_df, title, path_prefix) {
  plot_df <- prop_df
  plot_df$sample_id <- factor(plot_df$sample_id, levels = sample_order)
  plot_df$group <- factor(plot_df$group, levels = sort_manual_levels(plot_df$group))

  p <- ggplot(plot_df, aes(x = group, y = percent, fill = sample_id)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.7) +
    scale_fill_manual(values = sample_colors, drop = FALSE) +
    labs(x = "Manual cell type", y = "Cells (%)", fill = NULL, title = title) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "top"
    )

  ggsave(paste0(path_prefix, ".pdf"), p, width = 8.5, height = 4.8)
  ggsave(paste0(path_prefix, ".png"), p, width = 8.5, height = 4.8, dpi = 180)
  invisible(p)
}

plot_ratio <- function(comparison_df, title, path_prefix) {
  plot_df <- comparison_df
  plot_df$group <- factor(plot_df$group, levels = sort_manual_levels(plot_df$group))

  p <- ggplot(plot_df, aes(x = group, y = ratio_percent_TH2_over_TH1)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.35) +
    geom_col(width = 0.72, fill = "#7B3294") +
    labs(x = "Manual cell type", y = "Ratio of cell percentage (TH2 / TH1)", title = title) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 35, hjust = 1)
    )

  ggsave(paste0(path_prefix, ".pdf"), p, width = 8.5, height = 4.8)
  ggsave(paste0(path_prefix, ".png"), p, width = 8.5, height = 4.8, dpi = 180)
  invisible(p)
}

plot_delta <- function(comparison_df, title, path_prefix) {
  plot_df <- comparison_df[order(abs(comparison_df$delta_percent_TH2_minus_TH1), decreasing = TRUE), , drop = FALSE]
  plot_df$group <- factor(plot_df$group, levels = rev(plot_df$group))

  p <- ggplot(plot_df, aes(x = group, y = delta_percent_TH2_minus_TH1, fill = delta_percent_TH2_minus_TH1 > 0)) +
    geom_col(width = 0.72) +
    coord_flip() +
    scale_fill_manual(values = c(`TRUE` = "#b2182b", `FALSE` = "#2166ac"), guide = "none") +
    labs(x = NULL, y = "Delta cells (% TH2 - % TH1)", title = title) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(paste0(path_prefix, ".pdf"), p, width = 7, height = 4.8)
  ggsave(paste0(path_prefix, ".png"), p, width = 7, height = 4.8, dpi = 180)
  invisible(p)
}

proj <- loadArchRProject(path = output_dir)
cell_metadata <- as.data.frame(getCellColData(proj))

if (!cluster_col %in% colnames(cell_metadata)) {
  stop("Missing `", cluster_col, "`. Run combined clustering before manual annotation.")
}

cell_annotation <- manual_annotation[match(cell_metadata[[cluster_col]], manual_annotation$cluster), ]
cell_metadata$manual_celltype <- cell_annotation$manual_celltype
cell_metadata$manual_celltype_full <- cell_annotation$manual_celltype_full
cell_metadata$manual_celltype[is.na(cell_metadata$manual_celltype)] <- "Unannotated"
cell_metadata$manual_celltype_full[is.na(cell_metadata$manual_celltype_full)] <- "Unannotated"

proj <- addCellColData(
  ArchRProj = proj,
  data = cell_metadata$manual_celltype,
  cells = rownames(cell_metadata),
  name = manual_col,
  force = TRUE
)
proj <- addCellColData(
  ArchRProj = proj,
  data = cell_metadata$manual_celltype_full,
  cells = rownames(cell_metadata),
  name = "ManualCelltypeFull",
  force = TRUE
)

p_manual <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = manual_col,
  embedding = "UMAP_Combined"
)
p_cluster <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = cluster_col,
  embedding = "UMAP_Combined"
)
p_sample <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "UMAP_Combined"
)

plotPDF(
  p_manual,
  p_cluster,
  p_sample,
  name = "manual_celltypes_clusters_samples.pdf",
  ArchRProj = proj,
  addDOC = FALSE,
  width = 6,
  height = 6
)
file.copy(
  file.path(output_dir, "Plots", "manual_celltypes_clusters_samples.pdf"),
  manual_dir,
  overwrite = TRUE
)
ggsave(file.path(manual_dir, "manual_celltypes_umap.png"), p_manual, width = 6, height = 6, dpi = 180)

sample_order <- unique(cell_metadata$Sample)
if (all(c("TH1", "TH2") %in% sample_order)) {
  sample_order <- c("TH1", "TH2")
}
cell_metadata$Sample <- factor(cell_metadata$Sample, levels = sample_order)

manual_props <- count_proportions(cell_metadata, "manual_celltype", "manual_celltype")
manual_comparison <- compare_two_samples(manual_props)
manual_counts <- as.data.frame(table(cell_metadata$manual_celltype, cell_metadata$Sample), stringsAsFactors = FALSE)
colnames(manual_counts) <- c("manual_celltype", "sample_id", "n_cells")
manual_counts <- manual_counts[manual_counts$n_cells > 0, , drop = FALSE]

write.csv(manual_annotation, file.path(output_dir, "manual_cluster_celltype_annotations.csv"), row.names = FALSE)
write.csv(cell_metadata, file.path(output_dir, "cell_metadata_with_manual_celltypes.csv"), row.names = TRUE)
write.csv(manual_counts, file.path(output_dir, "manual_celltype_counts_by_sample.csv"), row.names = FALSE)
write.csv(manual_props, file.path(output_dir, "manual_celltype_proportions_by_sample.csv"), row.names = FALSE)
write.csv(manual_comparison, file.path(output_dir, "manual_celltype_proportion_comparison_TH2_vs_TH1.csv"), row.names = FALSE)

plot_stacked_proportions(
  manual_props,
  title = "Manual cell type proportions",
  path_prefix = file.path(proportion_dir, "manual_celltype_proportions_stacked")
)
plot_grouped_percent(
  manual_props,
  title = "Manual cell type percentages by sample",
  path_prefix = file.path(proportion_dir, "manual_celltype_percent_by_sample_grouped")
)
plot_ratio(
  manual_comparison,
  title = "Manual cell type percentage ratio",
  path_prefix = file.path(proportion_dir, "manual_celltype_percent_ratio_TH2_over_TH1")
)
plot_delta(
  manual_comparison,
  title = "Manual cell type proportion difference",
  path_prefix = file.path(proportion_dir, "manual_celltype_proportion_delta")
)

saveArchRProject(ArchRProj = proj, outputDirectory = output_dir, load = FALSE)

message("Manual cell type annotation written to: ", file.path(output_dir, "manual_cluster_celltype_annotations.csv"))
message("Manual proportion plots written to: ", proportion_dir)
