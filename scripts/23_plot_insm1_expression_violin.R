#!/usr/bin/env Rscript

source(file.path("config", "analysis_config.R"))
source(file.path("R", "project_helpers.R"))

suppressPackageStartupMessages({
  library(ArchR)
  library(SummarizedExperiment)
  library(Matrix)
  library(ggplot2)
})

set.seed(random_seed)
setup_archr_session()

expression_dir <- file.path(figure_dir, "gene_expression_compare")
ensure_dirs(output_dir, expression_dir)

gene_symbol <- "Insm1"
sample_colors <- c(TH1 = "#4C78A8", TH2 = "#E45756")
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

get_feature_names <- function(se) {
  rd <- as.data.frame(rowData(se))
  if ("name" %in% colnames(rd)) {
    return(as.character(rd$name))
  }
  if ("symbol" %in% colnames(rd)) {
    return(as.character(rd$symbol))
  }
  rownames(se)
}

get_first_assay <- function(se) {
  assay_names <- names(assays(se))
  preferred <- c("GeneExpressionMatrix", "GeneScoreMatrix", "z", "mat", "logcounts", "counts")
  assay_name <- preferred[preferred %in% assay_names][1]
  if (is.na(assay_name)) {
    assay_name <- assay_names[1]
  }
  assay(se, assay_name)
}

wilcox_two_sample <- function(df) {
  samples <- unique(as.character(df$Sample))
  if (!all(c("TH1", "TH2") %in% samples)) {
    return(NA_real_)
  }
  stats::wilcox.test(Expression ~ Sample, data = df[df$Sample %in% c("TH1", "TH2"), , drop = FALSE])$p.value
}

plot_violin <- function(df, title, y_label, path_prefix, width = 4.8, height = 4.5) {
  df$Sample <- factor(df$Sample, levels = sample_order)

  p <- ggplot(df, aes(x = Sample, y = Expression, fill = Sample)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.25, color = "grey25") +
    geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.25, color = "grey20", fill = "white") +
    scale_fill_manual(values = sample_colors, drop = FALSE) +
    labs(x = NULL, y = y_label, title = title) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none"
    )

  ggsave(paste0(path_prefix, ".pdf"), p, width = width, height = height)
  ggsave(paste0(path_prefix, ".png"), p, width = width, height = height, dpi = 200)
  invisible(p)
}

plot_faceted_violin <- function(df, title, y_label, path_prefix) {
  df$Sample <- factor(df$Sample, levels = sample_order)
  df$manual_celltype <- factor(
    df$manual_celltype,
    levels = manual_celltype_order[manual_celltype_order %in% unique(df$manual_celltype)]
  )

  p <- ggplot(df, aes(x = Sample, y = Expression, fill = Sample)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.2, color = "grey30") +
    geom_boxplot(width = 0.14, outlier.shape = NA, linewidth = 0.2, color = "grey20", fill = "white") +
    facet_wrap(~manual_celltype, scales = "free_y", ncol = 5) +
    scale_fill_manual(values = sample_colors, drop = FALSE) +
    labs(x = NULL, y = y_label, title = title) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      strip.background = element_rect(fill = "white", color = "grey30", linewidth = 0.25),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "none"
    )

  ggsave(paste0(path_prefix, ".pdf"), p, width = 10, height = 6)
  ggsave(paste0(path_prefix, ".png"), p, width = 10, height = 6, dpi = 200)
  invisible(p)
}

plot_umap_expression_split_by_sample <- function(proj, df, gene_symbol, path_prefix) {
  embedding <- getEmbedding(proj, embedding = "UMAP_Combined", returnDF = TRUE)
  if (ncol(embedding) != 2) {
    stop("UMAP_Combined does not have exactly 2 columns.")
  }
  colnames(embedding) <- c("UMAP_1", "UMAP_2")

  common_cells <- intersect(rownames(embedding), df$cell_id)
  plot_df <- merge(
    data.frame(cell_id = common_cells, embedding[common_cells, , drop = FALSE], check.names = FALSE),
    df,
    by = "cell_id"
  )
  plot_df$Sample <- factor(plot_df$Sample, levels = sample_order)

  color_limits <- stats::quantile(plot_df$Expression, probs = c(0, 0.99), na.rm = TRUE)
  if (!all(is.finite(color_limits)) || color_limits[1] == color_limits[2]) {
    color_limits <- range(plot_df$Expression, na.rm = TRUE)
  }

  plot_df <- plot_df[order(plot_df$Expression), , drop = FALSE]

  p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = Expression)) +
    geom_point(size = 0.12, alpha = 0.9) +
    facet_wrap(~Sample, nrow = 1) +
    scale_color_gradientn(
      colors = c("#D9D9D9", "#E7B8A8", "#C15A44", "#7F0000"),
      limits = color_limits,
      oob = scales::squish,
      name = "Expression"
    ) +
    labs(x = "UMAP_1", y = "UMAP_2", title = paste0(gene_symbol, " RNA expression on combined UMAP")) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "white", color = "grey30", linewidth = 0.25),
      strip.text = element_text(face = "bold"),
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(linewidth = 0.25, color = "black")
    )

  ggsave(paste0(path_prefix, ".pdf"), p, width = 9.2, height = 4.6)
  ggsave(paste0(path_prefix, ".png"), p, width = 9.2, height = 4.6, dpi = 220)
  invisible(p)
}

plot_umap_expression_split_by_sample_with_background <- function(proj, df, gene_symbol, path_prefix) {
  embedding <- getEmbedding(proj, embedding = "UMAP_Combined", returnDF = TRUE)
  if (ncol(embedding) != 2) {
    stop("UMAP_Combined does not have exactly 2 columns.")
  }
  colnames(embedding) <- c("UMAP_1", "UMAP_2")

  all_df <- data.frame(cell_id = rownames(embedding), embedding, check.names = FALSE)
  common_cells <- intersect(rownames(embedding), df$cell_id)
  plot_df <- merge(
    data.frame(cell_id = common_cells, embedding[common_cells, , drop = FALSE], check.names = FALSE),
    df,
    by = "cell_id"
  )
  plot_df$Sample <- factor(plot_df$Sample, levels = sample_order)
  background_df <- do.call(rbind, lapply(sample_order, function(sample_id) {
    cbind(all_df, Sample = sample_id, stringsAsFactors = FALSE)
  }))
  background_df$Sample <- factor(background_df$Sample, levels = sample_order)

  color_limits <- stats::quantile(plot_df$Expression, probs = c(0, 0.99), na.rm = TRUE)
  if (!all(is.finite(color_limits)) || color_limits[1] == color_limits[2]) {
    color_limits <- range(plot_df$Expression, na.rm = TRUE)
  }

  plot_df <- plot_df[order(plot_df$Expression), , drop = FALSE]

  p <- ggplot() +
    geom_point(data = background_df, aes(x = UMAP_1, y = UMAP_2), color = "#E6E6E6", size = 0.08, alpha = 0.28) +
    geom_point(data = plot_df, aes(x = UMAP_1, y = UMAP_2, color = Expression), size = 0.12, alpha = 0.9) +
    facet_wrap(~Sample, nrow = 1) +
    scale_color_gradientn(
      colors = c("#D9D9D9", "#E7B8A8", "#C15A44", "#7F0000"),
      limits = color_limits,
      oob = scales::squish,
      name = "Expression"
    ) +
    labs(x = "UMAP_1", y = "UMAP_2", title = paste0(gene_symbol, " RNA expression split by sample")) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "white", color = "grey30", linewidth = 0.25),
      strip.text = element_text(face = "bold"),
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(linewidth = 0.25, color = "black")
    )

  ggsave(paste0(path_prefix, ".pdf"), p, width = 9.2, height = 4.6)
  ggsave(paste0(path_prefix, ".png"), p, width = 9.2, height = 4.6, dpi = 220)
  invisible(p)
}

proj <- loadArchRProject(path = output_dir)
cell_metadata <- as.data.frame(getCellColData(proj))
if (!"Sample" %in% colnames(cell_metadata)) {
  stop("Missing `Sample` in cell metadata.")
}
if (!"ManualCelltype" %in% colnames(cell_metadata)) {
  cell_metadata$ManualCelltype <- "Unannotated"
}

expr_matrix <- getMatrixFromProject(ArchRProj = proj, useMatrix = "GeneExpressionMatrix")
expr <- get_first_assay(expr_matrix)
gene_names <- get_feature_names(expr_matrix)
rownames(expr) <- make.unique(gene_names)

gene_lookup <- stats::setNames(rownames(expr), tolower(rownames(expr)))
gene_key <- tolower(gene_symbol)
if (!gene_key %in% names(gene_lookup)) {
  stop("Gene not found in GeneExpressionMatrix: ", gene_symbol)
}

matrix_gene_name <- gene_lookup[[gene_key]]
expression_values <- as.numeric(expr[matrix_gene_name, , drop = TRUE])
names(expression_values) <- colnames(expr)

common_cells <- intersect(rownames(cell_metadata), names(expression_values))
plot_df <- data.frame(
  cell_id = common_cells,
  Sample = cell_metadata[common_cells, "Sample"],
  manual_celltype = cell_metadata[common_cells, "ManualCelltype"],
  Expression = expression_values[common_cells],
  stringsAsFactors = FALSE
)
plot_df$manual_celltype[is.na(plot_df$manual_celltype)] <- "Unannotated"

sample_order <- unique(plot_df$Sample)
if (all(c("TH1", "TH2") %in% sample_order)) {
  sample_order <- c("TH1", "TH2")
}

summary_by_sample <- aggregate(
  Expression ~ Sample,
  data = plot_df,
  FUN = function(x) {
    c(
      n_cells = length(x),
      mean = mean(x),
      median = stats::median(x),
      pct_detected = mean(x > 0) * 100
    )
  }
)
summary_by_sample <- do.call(data.frame, summary_by_sample)
colnames(summary_by_sample) <- c("sample_id", "n_cells", "mean_expression", "median_expression", "percent_detected")
summary_by_sample$wilcox_p_value_TH2_vs_TH1 <- wilcox_two_sample(plot_df)

summary_by_celltype <- aggregate(
  Expression ~ manual_celltype + Sample,
  data = plot_df,
  FUN = function(x) {
    c(
      n_cells = length(x),
      mean = mean(x),
      median = stats::median(x),
      pct_detected = mean(x > 0) * 100
    )
  }
)
summary_by_celltype <- do.call(data.frame, summary_by_celltype)
colnames(summary_by_celltype) <- c("manual_celltype", "sample_id", "n_cells", "mean_expression", "median_expression", "percent_detected")

write.csv(plot_df, file.path(output_dir, "insm1_expression_by_cell.csv"), row.names = FALSE)
write.csv(summary_by_sample, file.path(output_dir, "insm1_expression_summary_by_sample.csv"), row.names = FALSE)
write.csv(summary_by_celltype, file.path(output_dir, "insm1_expression_summary_by_manual_celltype_sample.csv"), row.names = FALSE)

plot_violin(
  plot_df,
  title = "Insm1 RNA expression by sample",
  y_label = "Insm1 expression",
  path_prefix = file.path(expression_dir, "Insm1_violin_by_sample")
)
plot_faceted_violin(
  plot_df,
  title = "Insm1 RNA expression by manual cell type and sample",
  y_label = "Insm1 expression",
  path_prefix = file.path(expression_dir, "Insm1_violin_by_manual_celltype_sample")
)
plot_umap_expression_split_by_sample(
  proj = proj,
  df = plot_df,
  gene_symbol = gene_symbol,
  path_prefix = file.path(expression_dir, "Insm1_umap_expression_split_by_sample")
)
plot_umap_expression_split_by_sample_with_background(
  proj = proj,
  df = plot_df,
  gene_symbol = gene_symbol,
  path_prefix = file.path(expression_dir, "Insm1_umap_expression_split_by_sample_with_background")
)

message("Insm1 expression summaries written to: ", output_dir)
message("Insm1 expression plots written to: ", expression_dir)
