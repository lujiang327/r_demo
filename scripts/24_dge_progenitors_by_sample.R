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

dge_dir <- file.path("results", "dge_progenitors_by_sample")
figure_dge_dir <- file.path(figure_dir, "dge_progenitors_by_sample")
ensure_dirs(output_dir, dge_dir, figure_dge_dir)

manual_col <- "ManualCelltype"
sample_col <- "Sample"
sample_a <- "TH1"
sample_b <- "TH2"

# Output label nRPC follows the user's wording; project metadata label is nPRC.
target_celltypes <- data.frame(
  output_label = c("pRPC", "nRPC"),
  manual_celltype = c("pRPC", "nPRC"),
  stringsAsFactors = FALSE
)

min_pct_detected <- 0.05
padj_cutoff <- 0.05
logfc_cutoff <- 0.25

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

run_presto_dge <- function(expr_subset, samples) {
  if (!requireNamespace("presto", quietly = TRUE)) {
    stop("Package `presto` is required for this script.")
  }

  # presto expects features x cells for matrix input.
  dge <- presto::wilcoxauc(
    X = expr_subset,
    y = samples
  )
  dge <- dge[dge$group == sample_b, , drop = FALSE]

  col_rename <- c(
    feature = "gene",
    logFC = "log2FC_TH2_vs_TH1",
    pval = "p_value",
    padj = "padj",
    pct_in = "pct_detected_TH2",
    pct_out = "pct_detected_TH1",
    avgExpr = "mean_expression_TH2"
  )
  for (old_name in names(col_rename)) {
    if (old_name %in% colnames(dge)) {
      colnames(dge)[colnames(dge) == old_name] <- col_rename[[old_name]]
    }
  }

  dge
}

add_group_summaries <- function(dge, expr_subset, samples) {
  th1_cells <- names(samples)[samples == sample_a]
  th2_cells <- names(samples)[samples == sample_b]

  mean_th1 <- Matrix::rowMeans(expr_subset[, th1_cells, drop = FALSE])
  mean_th2 <- Matrix::rowMeans(expr_subset[, th2_cells, drop = FALSE])
  pct_th1 <- Matrix::rowMeans(expr_subset[, th1_cells, drop = FALSE] > 0)
  pct_th2 <- Matrix::rowMeans(expr_subset[, th2_cells, drop = FALSE] > 0)

  summary_df <- data.frame(
    gene = rownames(expr_subset),
    mean_expression_TH1 = as.numeric(mean_th1),
    mean_expression_TH2 = as.numeric(mean_th2),
    pct_detected_TH1 = as.numeric(pct_th1),
    pct_detected_TH2 = as.numeric(pct_th2),
    stringsAsFactors = FALSE
  )

  summary_cols <- c(
    "mean_expression_TH1",
    "mean_expression_TH2",
    "pct_detected_TH1",
    "pct_detected_TH2",
    "max_pct_detected",
    "delta_pct_detected_TH2_minus_TH1"
  )
  dge <- dge[, setdiff(colnames(dge), summary_cols), drop = FALSE]

  out <- merge(dge, summary_df, by = "gene", all.x = TRUE)
  if (!"log2FC_TH2_vs_TH1" %in% colnames(out)) {
    out$log2FC_TH2_vs_TH1 <- log2((out$mean_expression_TH2 + 0.01) / (out$mean_expression_TH1 + 0.01))
  }
  out$max_pct_detected <- pmax(out$pct_detected_TH1, out$pct_detected_TH2, na.rm = TRUE)
  out$delta_pct_detected_TH2_minus_TH1 <- out$pct_detected_TH2 - out$pct_detected_TH1
  out$direction <- ifelse(
    out$padj < padj_cutoff & out$log2FC_TH2_vs_TH1 >= logfc_cutoff & out$max_pct_detected >= min_pct_detected,
    "Higher in TH2",
    ifelse(
      out$padj < padj_cutoff & out$log2FC_TH2_vs_TH1 <= -logfc_cutoff & out$max_pct_detected >= min_pct_detected,
      "Higher in TH1",
      "Not significant"
    )
  )
  out[order(out$padj, -abs(out$log2FC_TH2_vs_TH1)), , drop = FALSE]
}

plot_volcano <- function(dge, title, path_prefix) {
  plot_df <- dge
  plot_df$neg_log10_padj <- -log10(pmax(plot_df$padj, .Machine$double.xmin))
  plot_df$direction <- factor(plot_df$direction, levels = c("Higher in TH1", "Not significant", "Higher in TH2"))

  top_label_df <- rbind(
    head(plot_df[plot_df$direction == "Higher in TH2", ][order(plot_df[plot_df$direction == "Higher in TH2", "padj"]), ], 8),
    head(plot_df[plot_df$direction == "Higher in TH1", ][order(plot_df[plot_df$direction == "Higher in TH1", "padj"]), ], 8)
  )

  p <- ggplot(plot_df, aes(x = log2FC_TH2_vs_TH1, y = neg_log10_padj, color = direction)) +
    geom_point(size = 0.45, alpha = 0.75) +
    geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", color = "grey55", linewidth = 0.3) +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey55", linewidth = 0.3) +
    geom_text(
      data = top_label_df,
      aes(label = gene),
      size = 2.5,
      check_overlap = TRUE,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c("Higher in TH1" = "#4C78A8", "Not significant" = "#B8B8B8", "Higher in TH2" = "#E45756"),
      drop = FALSE
    ) +
    labs(x = "log2FC (TH2 / TH1)", y = "-log10 adjusted p value", color = NULL, title = title) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "top")

  ggsave(paste0(path_prefix, ".pdf"), p, width = 7, height = 5.2)
  ggsave(paste0(path_prefix, ".png"), p, width = 7, height = 5.2, dpi = 200)
  invisible(p)
}

plot_top_genes <- function(dge, title, path_prefix, n_top = 20) {
  sig_df <- dge[dge$direction != "Not significant", , drop = FALSE]
  if (nrow(sig_df) == 0) {
    return(invisible(NULL))
  }

  top_up <- head(sig_df[sig_df$direction == "Higher in TH2", ][order(sig_df[sig_df$direction == "Higher in TH2", "padj"]), ], n_top)
  top_down <- head(sig_df[sig_df$direction == "Higher in TH1", ][order(sig_df[sig_df$direction == "Higher in TH1", "padj"]), ], n_top)
  plot_df <- rbind(top_up, top_down)
  plot_df <- plot_df[order(plot_df$log2FC_TH2_vs_TH1), , drop = FALSE]
  plot_df$gene <- factor(plot_df$gene, levels = plot_df$gene)

  p <- ggplot(plot_df, aes(x = gene, y = log2FC_TH2_vs_TH1, fill = direction)) +
    geom_col(width = 0.75) +
    coord_flip() +
    scale_fill_manual(values = c("Higher in TH1" = "#4C78A8", "Higher in TH2" = "#E45756"), drop = FALSE) +
    labs(x = NULL, y = "log2FC (TH2 / TH1)", fill = NULL, title = title) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), legend.position = "top")

  height <- max(4.5, min(9, 0.2 * nrow(plot_df) + 1.5))
  ggsave(paste0(path_prefix, ".pdf"), p, width = 7, height = height)
  ggsave(paste0(path_prefix, ".png"), p, width = 7, height = height, dpi = 200)
  invisible(p)
}

proj <- loadArchRProject(path = output_dir)
cell_metadata <- as.data.frame(getCellColData(proj))
required_cols <- c(manual_col, sample_col)
missing_cols <- setdiff(required_cols, colnames(cell_metadata))
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "), ". Run script 20 first.")
}

expr_matrix <- getMatrixFromProject(ArchRProj = proj, useMatrix = "GeneExpressionMatrix")
expr <- get_first_assay(expr_matrix)
gene_names <- get_feature_names(expr_matrix)
rownames(expr) <- make.unique(gene_names)

summary_rows <- list()
for (i in seq_len(nrow(target_celltypes))) {
  output_label <- target_celltypes$output_label[[i]]
  manual_celltype <- target_celltypes$manual_celltype[[i]]

  target_cells <- rownames(cell_metadata)[cell_metadata[[manual_col]] == manual_celltype]
  target_cells <- intersect(target_cells, colnames(expr))
  samples <- as.character(cell_metadata[target_cells, sample_col])
  names(samples) <- target_cells

  if (!all(c(sample_a, sample_b) %in% samples)) {
    warning("Skipping ", output_label, ": both samples are not present.")
    next
  }

  expr_subset <- expr[, target_cells, drop = FALSE]
  keep_gene <- Matrix::rowMeans(expr_subset > 0) >= min_pct_detected
  expr_subset <- expr_subset[keep_gene, , drop = FALSE]

  message("Running DGE for ", output_label, " (", manual_celltype, "): ", ncol(expr_subset), " cells, ", nrow(expr_subset), " genes.")
  dge <- run_presto_dge(expr_subset, samples)
  dge <- add_group_summaries(dge, expr_subset, samples)

  dge$celltype_output_label <- output_label
  dge$manual_celltype <- manual_celltype
  dge$n_cells_TH1 <- sum(samples == sample_a)
  dge$n_cells_TH2 <- sum(samples == sample_b)

  all_path <- file.path(dge_dir, paste0(output_label, "_TH2_vs_TH1_rna_dge_all.csv"))
  sig_path <- file.path(dge_dir, paste0(output_label, "_TH2_vs_TH1_rna_dge_significant.csv"))
  write.csv(dge, all_path, row.names = FALSE)
  write.csv(dge[dge$direction != "Not significant", , drop = FALSE], sig_path, row.names = FALSE)

  plot_volcano(
    dge,
    title = paste0(output_label, " RNA DGE: TH2 vs TH1"),
    path_prefix = file.path(figure_dge_dir, paste0(output_label, "_TH2_vs_TH1_volcano"))
  )
  plot_top_genes(
    dge,
    title = paste0(output_label, " top RNA DGE genes"),
    path_prefix = file.path(figure_dge_dir, paste0(output_label, "_TH2_vs_TH1_top_genes"))
  )

  summary_rows[[output_label]] <- data.frame(
    celltype_output_label = output_label,
    manual_celltype = manual_celltype,
    n_cells_TH1 = sum(samples == sample_a),
    n_cells_TH2 = sum(samples == sample_b),
    n_genes_tested = nrow(dge),
    n_higher_TH2 = sum(dge$direction == "Higher in TH2"),
    n_higher_TH1 = sum(dge$direction == "Higher in TH1"),
    all_results_csv = all_path,
    significant_results_csv = sig_path,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(dge_dir, "progenitor_dge_summary.csv"), row.names = FALSE)

message("DGE tables written to: ", dge_dir)
message("DGE figures written to: ", figure_dge_dir)
