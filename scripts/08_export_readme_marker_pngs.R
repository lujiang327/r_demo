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

marker_dir <- file.path("results", "readme_figures", "markers")
ensure_dirs(marker_dir)

proj <- loadArchRProject(path = output_dir)

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

marker_genes_to_plot <- unique(c(
  "Rho", "Gnat1", "Nrl", "Arr3", "Opn1mw", "Opn1sw",
  "Vsx1", "Vsx2", "Car10", "Prkca", "Sebox", "Scgn", "Cabp5", "Grm6",
  "Tfap2a", "Gad1", "Gad2",
  "Rlbp1", "Glul", "Dkk3", "Clu", "Apoe", "Aqp4",
  "Calb1", "Onecut1", "Rbpms", "Pou4f1", "Pou4f2",
  "Ttr", "Rdh5", "Rpe65", "S100b", "Gfap", "Pax2", "Pdgfra",
  "Ctss", "C1qa", "Cldn5", "Flt1", "Kcnj8", "Pecam1", "Acta2",
  "Cdk1", "Mki67", "Top2a", "Pcna",
  "Ascl1", "Neurog2", "Insm1", "Atoh7", "Neurod1", "Otx2", "Crx", "Olig2", "Foxn4",
  "Fos", "Jun", "Stat3", "Lcn2", "Malat1"
))

expr_matrix <- getMatrixFromProject(proj, useMatrix = "GeneExpressionMatrix")
log_norm_expr <- get_first_assay(expr_matrix)
gene_names <- get_feature_names(expr_matrix)
rownames(log_norm_expr) <- make.unique(gene_names)

embedding <- getEmbedding(proj, embedding = "UMAP_Combined", returnDF = TRUE)
if (ncol(embedding) != 2) {
  stop("UMAP_Combined does not have exactly 2 columns.")
}
colnames(embedding) <- c("UMAP_1", "UMAP_2")

gene_lookup <- stats::setNames(rownames(log_norm_expr), tolower(rownames(log_norm_expr)))
summary_rows <- list()

for (gene_symbol in marker_genes_to_plot) {
  gene_key <- tolower(gene_symbol)
  if (!gene_key %in% names(gene_lookup)) {
    warning("Gene not found in expression matrix: ", gene_symbol)
    next
  }

  matrix_gene_name <- gene_lookup[[gene_key]]
  expr_values <- as.numeric(log_norm_expr[matrix_gene_name, , drop = TRUE])
  names(expr_values) <- colnames(log_norm_expr)

  common_cells <- intersect(rownames(embedding), names(expr_values))
  if (length(common_cells) == 0) {
    warning("No common cells for gene: ", gene_symbol)
    next
  }

  expr_sub <- expr_values[common_cells]
  expr_sd <- stats::sd(expr_sub)
  if (is.na(expr_sd) || expr_sd == 0) {
    expr_zscore <- rep(0, length(expr_sub))
  } else {
    expr_zscore <- (expr_sub - mean(expr_sub)) / expr_sd
  }

  color_limits <- stats::quantile(expr_zscore, probs = c(0.01, 0.99), na.rm = TRUE)
  if (!all(is.finite(color_limits)) || color_limits[1] == color_limits[2]) {
    color_limits <- range(expr_zscore, na.rm = TRUE)
  }

  df <- data.frame(
    embedding[common_cells, , drop = FALSE],
    ExpressionZ = expr_zscore
  )
  df <- df[order(df$ExpressionZ), , drop = FALSE]

  p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = ExpressionZ)) +
    geom_point(size = 0.2) +
    scale_color_gradientn(
      colors = c("#D9D9D9", "#E7B8A8", "#C15A44", "#7F0000"),
      limits = color_limits,
      oob = scales::squish,
      name = "ExpressionZ"
    ) +
    theme_classic(base_size = 11) +
    ggtitle(gene_symbol) +
    labs(x = "UMAP_1", y = "UMAP_2") +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(linewidth = 0.25, color = "black"),
      plot.title = element_text(size = 14),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )

  png_path <- file.path(marker_dir, paste0(gene_symbol, ".png"))
  ggsave(filename = png_path, plot = p, width = 3.2, height = 3.2, dpi = 150)

  summary_rows[[gene_symbol]] <- data.frame(
    gene = gene_symbol,
    png = png_path,
    min_z = min(expr_zscore),
    max_z = max(expr_zscore),
    stringsAsFactors = FALSE
  )
}

marker_png_summary <- do.call(rbind, summary_rows)
write.csv(
  marker_png_summary,
  file = file.path("results", "readme_figures", "marker_png_summary.csv"),
  row.names = FALSE
)

message("Wrote marker PNGs to: ", marker_dir)
