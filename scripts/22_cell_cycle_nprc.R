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

cell_cycle_dir <- file.path(figure_dir, "cell_cycle_nprc")
ensure_dirs(output_dir, cell_cycle_dir)

manual_col <- "ManualCelltype"
target_celltype <- "nPRC"
phase_order <- c("G1", "S", "G2M")
sample_colors <- c(TH1 = "#4C78A8", TH2 = "#E45756")
phase_colors <- c(G1 = "#B8B8B8", S = "#4C78A8", G2M = "#E45756")

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

human_to_mouse_symbol <- function(x) {
  paste0(substr(x, 1, 1), tolower(substr(x, 2, nchar(x))))
}

z_score_rows <- function(mat) {
  dense <- as.matrix(mat)
  row_means <- rowMeans(dense)
  row_sds <- apply(dense, 1, stats::sd)
  row_sds[is.na(row_sds) | row_sds == 0] <- 1
  sweep(sweep(dense, 1, row_means, "-"), 1, row_sds, "/")
}

count_proportions <- function(metadata, group_col, stage) {
  counts <- as.data.frame(table(metadata$Sample, metadata[[group_col]]), stringsAsFactors = FALSE)
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
  groups <- phase_order[phase_order %in% unique(prop_df$group)]
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

plot_stacked_phase <- function(prop_df, path_prefix) {
  plot_df <- prop_df
  plot_df$sample_id <- factor(plot_df$sample_id, levels = sample_order)
  plot_df$group <- factor(plot_df$group, levels = phase_order)

  p <- ggplot(plot_df, aes(x = sample_id, y = percent, fill = group)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.2) +
    scale_fill_manual(values = phase_colors, drop = FALSE) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(x = NULL, y = "nPRC cells (%)", fill = "Phase", title = "nPRC cell cycle phase proportions") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(paste0(path_prefix, ".pdf"), p, width = 6.5, height = 4.5)
  ggsave(paste0(path_prefix, ".png"), p, width = 6.5, height = 4.5, dpi = 180)
  invisible(p)
}

plot_grouped_phase <- function(prop_df, path_prefix) {
  plot_df <- prop_df
  plot_df$sample_id <- factor(plot_df$sample_id, levels = sample_order)
  plot_df$group <- factor(plot_df$group, levels = phase_order)

  p <- ggplot(plot_df, aes(x = group, y = percent, fill = sample_id)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.7) +
    scale_fill_manual(values = sample_colors, drop = FALSE) +
    labs(x = "Cell cycle phase", y = "nPRC cells (%)", fill = NULL, title = "nPRC cell cycle phase by sample") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "top")

  ggsave(paste0(path_prefix, ".pdf"), p, width = 6.5, height = 4.5)
  ggsave(paste0(path_prefix, ".png"), p, width = 6.5, height = 4.5, dpi = 180)
  invisible(p)
}

plot_score_scatter <- function(metadata, path_prefix) {
  plot_df <- metadata
  plot_df$CellCyclePhase <- factor(plot_df$CellCyclePhase, levels = phase_order)

  p <- ggplot(plot_df, aes(x = S.Score, y = G2M.Score, color = CellCyclePhase)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_point(size = 0.25, alpha = 0.7) +
    scale_color_manual(values = phase_colors, drop = FALSE) +
    facet_wrap(~Sample) +
    labs(x = "S score", y = "G2M score", color = "Phase", title = "nPRC cell cycle scores") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(paste0(path_prefix, ".pdf"), p, width = 7.5, height = 4.2)
  ggsave(paste0(path_prefix, ".png"), p, width = 7.5, height = 4.2, dpi = 180)
  invisible(p)
}

plot_nprc_phase_umap <- function(proj, nprc_metadata, path_prefix) {
  embedding <- getEmbedding(proj, embedding = "UMAP_Combined", returnDF = TRUE)
  if (ncol(embedding) != 2) {
    stop("UMAP_Combined does not have exactly 2 columns.")
  }
  colnames(embedding) <- c("UMAP_1", "UMAP_2")

  all_df <- data.frame(embedding, cell_id = rownames(embedding), stringsAsFactors = FALSE)
  nprc_df <- merge(
    all_df,
    data.frame(
      cell_id = rownames(nprc_metadata),
      Sample = nprc_metadata$Sample,
      CellCyclePhase = nprc_metadata$CellCyclePhase,
      stringsAsFactors = FALSE
    ),
    by = "cell_id"
  )
  nprc_df$CellCyclePhase <- factor(nprc_df$CellCyclePhase, levels = phase_order)

  p <- ggplot() +
    geom_point(data = all_df, aes(x = UMAP_1, y = UMAP_2), color = "#D9D9D9", size = 0.12, alpha = 0.45) +
    geom_point(data = nprc_df, aes(x = UMAP_1, y = UMAP_2, color = CellCyclePhase), size = 0.18, alpha = 0.9) +
    scale_color_manual(values = phase_colors, drop = FALSE) +
    labs(x = "UMAP_1", y = "UMAP_2", color = "Phase", title = "nPRC cell cycle phase on combined UMAP") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(linewidth = 0.25, color = "black")
    )

  ggsave(paste0(path_prefix, ".pdf"), p, width = 6, height = 5.5)
  ggsave(paste0(path_prefix, ".png"), p, width = 6, height = 5.5, dpi = 220)

  p_sample <- ggplot() +
    geom_point(data = all_df, aes(x = UMAP_1, y = UMAP_2), color = "#D9D9D9", size = 0.12, alpha = 0.45) +
    geom_point(data = nprc_df, aes(x = UMAP_1, y = UMAP_2, color = Sample), size = 0.18, alpha = 0.9) +
    scale_color_manual(values = sample_colors, drop = FALSE) +
    labs(x = "UMAP_1", y = "UMAP_2", color = NULL, title = "nPRC cells on combined UMAP by sample") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(linewidth = 0.25, color = "black")
    )

  sample_prefix <- sub("_phase$", "_sample", path_prefix)
  ggsave(paste0(sample_prefix, ".pdf"), p_sample, width = 6, height = 5.5)
  ggsave(paste0(sample_prefix, ".png"), p_sample, width = 6, height = 5.5, dpi = 220)

  invisible(list(phase = p, sample = p_sample))
}

proj <- loadArchRProject(path = output_dir)
cell_metadata <- as.data.frame(getCellColData(proj))

if (!manual_col %in% colnames(cell_metadata)) {
  stop("Missing `", manual_col, "`. Run script 20 before nPRC cell cycle scoring.")
}

nprc_cells <- rownames(cell_metadata)[cell_metadata[[manual_col]] == target_celltype]
if (length(nprc_cells) == 0) {
  stop("No cells annotated as ", target_celltype, ".")
}

expr_matrix <- getMatrixFromProject(ArchRProj = proj, useMatrix = "GeneExpressionMatrix")
expr <- get_first_assay(expr_matrix)
gene_names <- get_feature_names(expr_matrix)
rownames(expr) <- make.unique(gene_names)

if (requireNamespace("Seurat", quietly = TRUE)) {
  s_genes <- human_to_mouse_symbol(Seurat::cc.genes.updated.2019$s.genes)
  g2m_genes <- human_to_mouse_symbol(Seurat::cc.genes.updated.2019$g2m.genes)
} else {
  s_genes <- c("Mcm5", "Pcna", "Tyms", "Fen1", "Mcm2", "Mcm4", "Rrm1", "Ung", "Gins2", "Mcm6", "Cdca7", "Dtl", "Prim1", "Uhrf1", "Cenpu", "Hells", "Rfc2", "Rpa2", "Nasp", "Rad51ap1", "Gmnn", "Wdr76", "Slbp", "Ccne2", "Ubr7", "Pold3", "Msh2", "Atad2", "Rad51", "Rrm2", "Cdc45", "Cdc6", "Exo1", "Tipin", "Dscc1", "Blm", "Casp8ap2", "Usp1", "Clspn", "Pola1", "Chaf1b", "Brip1", "E2f8")
  g2m_genes <- c("Hmgb2", "Cdk1", "Nusap1", "Ube2c", "Birc5", "Tpx2", "Top2a", "Ndc80", "Cks2", "Nuf2", "Cks1b", "Mki67", "Tmpo", "Cenpf", "Tacc3", "Fam64a", "Smc4", "Ccnb2", "Ckap2l", "Ckap2", "Aurkb", "Bub1", "Kif11", "Anp32e", "Tubb4b", "Gtse1", "Kif20b", "Hjurp", "Cdca3", "Jpt1", "Cdc20", "Ttk", "Cdc25c", "Kif2c", "Rangap1", "Ncapd2", "Dlgap5", "Cdca2", "Cdca8", "Ect2", "Kif23", "Hmmr", "Aurka", "Psrc1", "Anln", "Lbr", "Ckap5", "Cenpe", "Ctcf", "Nek2", "G2e3", "Gas2l3", "Cbx5", "Cenpa")
}

s_genes_present <- intersect(s_genes, rownames(expr))
g2m_genes_present <- intersect(g2m_genes, rownames(expr))
score_genes <- unique(c(s_genes_present, g2m_genes_present))

if (length(s_genes_present) == 0 || length(g2m_genes_present) == 0) {
  stop("Could not find enough cell cycle genes in GeneExpressionMatrix.")
}

expr_nprc <- expr[score_genes, nprc_cells, drop = FALSE]
z_expr <- z_score_rows(expr_nprc)
s_score <- colMeans(z_expr[s_genes_present, , drop = FALSE])
g2m_score <- colMeans(z_expr[g2m_genes_present, , drop = FALSE])

phase <- ifelse(s_score <= 0 & g2m_score <= 0, "G1", ifelse(s_score >= g2m_score, "S", "G2M"))
cycling_status <- ifelse(phase == "G1", "G1_low_cycling_score", "Cycling")

nprc_metadata <- cell_metadata[nprc_cells, , drop = FALSE]
nprc_metadata$S.Score <- as.numeric(s_score[rownames(nprc_metadata)])
nprc_metadata$G2M.Score <- as.numeric(g2m_score[rownames(nprc_metadata)])
nprc_metadata$CellCyclePhase <- phase[rownames(nprc_metadata)]
nprc_metadata$CyclingStatus <- cycling_status[rownames(nprc_metadata)]

all_phase <- rep(NA_character_, nrow(cell_metadata))
names(all_phase) <- rownames(cell_metadata)
all_phase[rownames(nprc_metadata)] <- nprc_metadata$CellCyclePhase
all_s_score <- rep(NA_real_, nrow(cell_metadata))
names(all_s_score) <- rownames(cell_metadata)
all_s_score[rownames(nprc_metadata)] <- nprc_metadata$S.Score
all_g2m_score <- rep(NA_real_, nrow(cell_metadata))
names(all_g2m_score) <- rownames(cell_metadata)
all_g2m_score[rownames(nprc_metadata)] <- nprc_metadata$G2M.Score

proj <- addCellColData(
  ArchRProj = proj,
  data = all_phase,
  cells = names(all_phase),
  name = "nPRC_CellCyclePhase",
  force = TRUE
)
proj <- addCellColData(
  ArchRProj = proj,
  data = all_s_score,
  cells = names(all_s_score),
  name = "nPRC_S.Score",
  force = TRUE
)
proj <- addCellColData(
  ArchRProj = proj,
  data = all_g2m_score,
  cells = names(all_g2m_score),
  name = "nPRC_G2M.Score",
  force = TRUE
)

sample_order <- unique(nprc_metadata$Sample)
if (all(c("TH1", "TH2") %in% sample_order)) {
  sample_order <- c("TH1", "TH2")
}
nprc_metadata$Sample <- factor(nprc_metadata$Sample, levels = sample_order)
nprc_metadata$CellCyclePhase <- factor(nprc_metadata$CellCyclePhase, levels = phase_order)

phase_props <- count_proportions(nprc_metadata, "CellCyclePhase", "nPRC_cell_cycle_phase")
phase_comparison <- compare_two_samples(phase_props)
phase_counts <- as.data.frame(table(nprc_metadata$CellCyclePhase, nprc_metadata$Sample), stringsAsFactors = FALSE)
colnames(phase_counts) <- c("cell_cycle_phase", "sample_id", "n_cells")
phase_counts <- phase_counts[phase_counts$n_cells > 0, , drop = FALSE]

gene_set_summary <- data.frame(
  gene_set = c("S", "G2M"),
  n_input_genes = c(length(unique(s_genes)), length(unique(g2m_genes))),
  n_present_genes = c(length(s_genes_present), length(g2m_genes_present)),
  present_genes = c(paste(s_genes_present, collapse = ";"), paste(g2m_genes_present, collapse = ";")),
  stringsAsFactors = FALSE
)

write.csv(nprc_metadata, file.path(output_dir, "nprc_cell_cycle_metadata.csv"), row.names = TRUE)
write.csv(phase_counts, file.path(output_dir, "nprc_cell_cycle_phase_counts_by_sample.csv"), row.names = FALSE)
write.csv(phase_props, file.path(output_dir, "nprc_cell_cycle_phase_proportions_by_sample.csv"), row.names = FALSE)
write.csv(phase_comparison, file.path(output_dir, "nprc_cell_cycle_phase_comparison_TH2_vs_TH1.csv"), row.names = FALSE)
write.csv(gene_set_summary, file.path(output_dir, "nprc_cell_cycle_gene_set_summary.csv"), row.names = FALSE)

p_phase <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "nPRC_CellCyclePhase",
  embedding = "UMAP_Combined"
)
p_sample <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "UMAP_Combined",
  cells = nprc_cells
)
plotPDF(
  p_phase,
  p_sample,
  name = "nprc_cell_cycle_umap.pdf",
  ArchRProj = proj,
  addDOC = FALSE,
  width = 6,
  height = 6
)
file.copy(
  file.path(output_dir, "Plots", "nprc_cell_cycle_umap.pdf"),
  cell_cycle_dir,
  overwrite = TRUE
)

plot_stacked_phase(phase_props, file.path(cell_cycle_dir, "nprc_cell_cycle_phase_stacked"))
plot_grouped_phase(phase_props, file.path(cell_cycle_dir, "nprc_cell_cycle_phase_by_sample_grouped"))
plot_score_scatter(nprc_metadata, file.path(cell_cycle_dir, "nprc_cell_cycle_score_scatter"))
plot_nprc_phase_umap(proj, nprc_metadata, file.path(cell_cycle_dir, "nprc_cell_cycle_umap_phase"))

saveArchRProject(ArchRProj = proj, outputDirectory = output_dir, load = FALSE)

message("nPRC cell cycle metadata written to: ", file.path(output_dir, "nprc_cell_cycle_metadata.csv"))
message("nPRC cell cycle figures written to: ", cell_cycle_dir)
