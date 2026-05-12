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

annotation_dir <- file.path(figure_dir, "celltype_annotation")
ensure_dirs(annotation_dir, log_dir)

proj <- loadArchRProject(path = output_dir)
cluster_col <- "Clusters_Combined"
if (!cluster_col %in% colnames(as.data.frame(getCellColData(proj)))) {
  stop("Missing `", cluster_col, "`. Run script 04 before this annotation script.")
}

# Primary identity markers come from the prior retina scRNA-seq annotation panel.
identity_marker_sets <- list(
  Rod = c("Rho", "Gnat1", "Pde6a", "Pde6b", "Cnga1", "Nr2e3", "Nrl", "Prph2", "Rom1"),
  Cone = c("Pde6h", "Arr3", "Opn1mw", "Gnat2", "Pde6c", "Opn1sw"),
  BC = c("Vsx1", "Vsx2", "Car10", "Prkca", "Sebox", "Scgn", "Cabp5", "Grm6", "Otx2os1", "Gng13", "Nrxn3", "Gabrb3", "Trnp1", "Kcnma1", "Frmd3", "Gm4792", "Nyap2"),
  AC = c("Tfap2a", "Tfap2b", "Pax6", "Frmd5", "Nrg3", "Elavl3", "Gad1", "Gad2"),
  MG = c("Rlbp1", "Glul", "Dkk3", "Adamtsl1", "Gpr37", "Abca8a", "Rgs6", "Spc25", "Clu", "Apoe", "Aqp4"),
  HC = c("Calb1", "Onecut1", "Slc4a3", "Onecut2", "Gm45459", "C1ql1"),
  RGC = c("Nefl", "Stmn2", "Nrn1", "Pou4f1", "Pou4f2", "Sncg", "Rbpms"),
  RPE = c("Ttr", "Rdh5", "Rpe65", "Rgr", "Slc16a8"),
  Astrocyte = c("S100b", "Gfap", "Pax2", "Pdgfra", "Mlc1", "Prdx6"),
  Microglia = c("Ctss", "C1qa", "Hexb", "Trem2"),
  Endothelial = c("Cldn5", "Flt1", "Pecam1", "Ly6c1", "Ptprb"),
  Pericyte = c("Kcnj8", "Acta2", "Pdgfra")
)

state_marker_sets <- list(
  Activated_MG = c("Gfap", "Vim", "Lcn2", "Stat3", "Socs3", "Cxcl10", "Ccl2", "Serpina3n"),
  Proliferation = c("Cdk1", "Mki67", "Top2a", "Pcna", "Mcm5"),
  Neurogenic_TF = c("Ascl1", "Neurog2", "Insm1", "Atoh7", "Neurod1", "Otx2", "Crx", "Olig2", "Foxn4"),
  Neurogenic_Progenitor = c("Ascl1", "Neurog2", "Neurod1", "Atoh7", "Sox2", "Hes1", "Hes5", "Insm1", "Prdm1", "Otx2", "Crx", "Olig2", "Foxn4"),
  Stress_Response = c("Fos", "Jun", "Atf3", "Ddit3", "Hspa1a", "Hspa1b", "Gadd45b"),
  Low_Quality_Diagnostic = c("Malat1")
)
marker_sets <- c(identity_marker_sets, state_marker_sets)
identity_marker_set_names <- names(identity_marker_sets)

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

module_scores_by_cluster <- function(proj, matrix_name, marker_sets, cluster_col) {
  se <- getMatrixFromProject(ArchRProj = proj, useMatrix = matrix_name)
  mat <- get_first_assay(se)
  feature_names <- get_feature_names(se)
  rownames(mat) <- make.unique(feature_names)

  cell_metadata <- as.data.frame(getCellColData(proj))
  cell_metadata <- cell_metadata[colnames(mat), , drop = FALSE]
  clusters <- as.character(cell_metadata[[cluster_col]])

  scores <- lapply(names(marker_sets), function(module_name) {
    genes <- intersect(marker_sets[[module_name]], rownames(mat))
    if (length(genes) == 0) {
      return(data.frame())
    }

    module_score <- Matrix::colMeans(mat[genes, , drop = FALSE])
    aggregate(
      module_score,
      by = list(cluster = clusters),
      FUN = mean
    ) |>
      transform(module = module_name, n_markers_used = length(genes)) |>
      subset(select = c("cluster", "module", "x", "n_markers_used"))
  })

  out <- do.call(rbind, scores)
  colnames(out)[colnames(out) == "x"] <- paste0(matrix_name, "_mean_score")
  out
}

assign_top_modules <- function(score_df, score_col, allowed_modules) {
  score_df <- score_df[score_df$module %in% allowed_modules, , drop = FALSE]
  split_scores <- split(score_df, score_df$cluster)

  do.call(rbind, lapply(split_scores, function(df) {
    df <- df[order(df[[score_col]], decreasing = TRUE), , drop = FALSE]
    data.frame(
      cluster = df$cluster[1],
      tentative_celltype = df$module[1],
      top_score = df[[score_col]][1],
      second_module = if (nrow(df) >= 2) df$module[2] else NA_character_,
      second_score = if (nrow(df) >= 2) df[[score_col]][2] else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
}

write_marker_table <- function(marker_sets, path) {
  marker_df <- do.call(rbind, lapply(names(marker_sets), function(module_name) {
    data.frame(module = module_name, gene = marker_sets[[module_name]], stringsAsFactors = FALSE)
  }))
  write.csv(marker_df, file = path, row.names = FALSE)
}

summarize_marker_presence <- function(proj, matrix_name, marker_sets, path) {
  se <- getMatrixFromProject(ArchRProj = proj, useMatrix = matrix_name)
  mat <- get_first_assay(se)
  feature_names <- get_feature_names(se)
  rownames(mat) <- make.unique(feature_names)
  marker_df <- do.call(rbind, lapply(names(marker_sets), function(module_name) {
    genes <- marker_sets[[module_name]]
    present <- genes %in% rownames(mat)
    total_counts <- rep(NA_real_, length(genes))
    mean_detected <- rep(NA_real_, length(genes))

    if (any(present)) {
      total_counts[present] <- Matrix::rowSums(mat[genes[present], , drop = FALSE])
      mean_detected[present] <- Matrix::rowMeans(mat[genes[present], , drop = FALSE] > 0)
    }

    data.frame(
      matrix = matrix_name,
      module = module_name,
      gene = genes,
      present = present,
      total_counts_or_score = total_counts,
      fraction_cells_detected = mean_detected,
      stringsAsFactors = FALSE
    )
  }))
  write.csv(marker_df, file = path, row.names = FALSE)
  invisible(marker_df)
}

build_rna_expression_z_plots <- function(proj, genes, embedding_name = "UMAP_Combined") {
  se <- getMatrixFromProject(ArchRProj = proj, useMatrix = "GeneExpressionMatrix")
  mat <- get_first_assay(se)
  feature_names <- get_feature_names(se)
  rownames(mat) <- make.unique(feature_names)

  embedding <- getEmbedding(proj, embedding = embedding_name, returnDF = TRUE)
  if (ncol(embedding) != 2) {
    stop("Embedding `", embedding_name, "` does not have exactly 2 columns.")
  }
  colnames(embedding) <- c("UMAP_1", "UMAP_2")

  gene_lookup <- stats::setNames(rownames(mat), tolower(rownames(mat)))
  available_genes <- genes[tolower(genes) %in% names(gene_lookup)]

  z_summary <- vector("list", length(available_genes))
  plots <- vector("list", length(available_genes))
  for (i in seq_along(available_genes)) {
    gene <- available_genes[[i]]
    matrix_gene_name <- gene_lookup[[tolower(gene)]]
    expression_values <- as.numeric(mat[matrix_gene_name, , drop = TRUE])
    names(expression_values) <- colnames(mat)

    common_cells <- intersect(rownames(embedding), names(expression_values))
    if (length(common_cells) == 0) {
      warning("No common cells between embedding and expression matrix for gene: ", gene)
      next
    }

    expression_values <- expression_values[common_cells]
    expression_sd <- stats::sd(expression_values)
    expression_mean <- mean(expression_values)

    if (is.na(expression_sd) || expression_sd == 0) {
      z_values <- rep(0, length(expression_values))
    } else {
      z_values <- (expression_values - expression_mean) / expression_sd
    }

    color_limits <- stats::quantile(z_values, probs = c(0.01, 0.99), na.rm = TRUE)
    if (!all(is.finite(color_limits)) || color_limits[1] == color_limits[2]) {
      color_limits <- range(z_values, na.rm = TRUE)
    }

    plot_df <- data.frame(
      embedding[common_cells, , drop = FALSE],
      ExpressionZ = z_values
    )
    plot_df <- plot_df[order(plot_df$ExpressionZ), , drop = FALSE]

    z_summary[[i]] <- data.frame(
      gene = gene,
      matrix_gene_name = matrix_gene_name,
      expression_mean = expression_mean,
      expression_sd = expression_sd,
      min_z = min(z_values),
      max_z = max(z_values),
      stringsAsFactors = FALSE
    )

    plots[[i]] <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = ExpressionZ)) +
      geom_point(size = 0.5) +
      scale_color_gradientn(
        colors = c("#D9D9D9", "#E7B8A8", "#C15A44", "#7F0000"),
        limits = color_limits,
        oob = scales::squish,
        name = "ExpressionZ"
      ) +
      theme_classic(base_size = 12) +
      ggtitle(paste0(project_name, " - ", gene)) +
      labs(x = "UMAP_1", y = "UMAP_2") +
      theme(
        panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.3, color = "black"),
        axis.ticks = element_line(linewidth = 0.25, color = "black")
      )
  }

  keep <- !vapply(plots, is.null, logical(1))
  list(
    plots = plots[keep],
    z_summary = if (any(keep)) do.call(rbind, z_summary[keep]) else data.frame()
  )
}

write_marker_table(marker_sets, file.path(output_dir, "retinal_marker_sets.csv"))
summarize_marker_presence(
  proj = proj,
  matrix_name = "GeneExpressionMatrix",
  marker_sets = marker_sets,
  path = file.path(output_dir, "marker_presence_gene_expression_matrix.csv")
)
summarize_marker_presence(
  proj = proj,
  matrix_name = "GeneScoreMatrix",
  marker_sets = marker_sets,
  path = file.path(output_dir, "marker_presence_gene_score_matrix.csv")
)

rna_scores <- module_scores_by_cluster(
  proj = proj,
  matrix_name = "GeneExpressionMatrix",
  marker_sets = marker_sets,
  cluster_col = cluster_col
)
write.csv(rna_scores, file = file.path(output_dir, "cluster_marker_module_scores_rna.csv"), row.names = FALSE)

gene_score_scores <- module_scores_by_cluster(
  proj = proj,
  matrix_name = "GeneScoreMatrix",
  marker_sets = marker_sets,
  cluster_col = cluster_col
)
write.csv(gene_score_scores, file = file.path(output_dir, "cluster_marker_module_scores_gene_score.csv"), row.names = FALSE)

rna_assignment <- assign_top_modules(
  rna_scores,
  score_col = "GeneExpressionMatrix_mean_score",
  allowed_modules = identity_marker_set_names
)
colnames(rna_assignment)[colnames(rna_assignment) != "cluster"] <- paste0(
  "rna_",
  colnames(rna_assignment)[colnames(rna_assignment) != "cluster"]
)

gene_score_assignment <- assign_top_modules(
  gene_score_scores,
  score_col = "GeneScoreMatrix_mean_score",
  allowed_modules = identity_marker_set_names
)
colnames(gene_score_assignment)[colnames(gene_score_assignment) != "cluster"] <- paste0(
  "gene_score_",
  colnames(gene_score_assignment)[colnames(gene_score_assignment) != "cluster"]
)

cluster_annotation <- merge(rna_assignment, gene_score_assignment, by = "cluster", all = TRUE)
cluster_annotation$tentative_celltype <- ifelse(
  is.na(cluster_annotation$rna_tentative_celltype),
  cluster_annotation$gene_score_tentative_celltype,
  cluster_annotation$rna_tentative_celltype
)
cluster_annotation$mg_candidate <- cluster_annotation$rna_tentative_celltype == "MG" |
  cluster_annotation$rna_second_module == "MG" |
  cluster_annotation$gene_score_tentative_celltype == "MG" |
  cluster_annotation$gene_score_second_module == "MG"
cluster_annotation$annotation_note <- ifelse(
  cluster_annotation$tentative_celltype != "MG" & cluster_annotation$mg_candidate,
  "MG marker signal present but another identity module scored higher; inspect manually.",
  ""
)
write.csv(cluster_annotation, file = file.path(output_dir, "cluster_tentative_celltype_annotations.csv"), row.names = FALSE)

cell_metadata <- as.data.frame(getCellColData(proj))
cell_annotation <- cluster_annotation[match(cell_metadata[[cluster_col]], cluster_annotation$cluster), ]
cell_metadata$tentative_celltype <- cell_annotation$tentative_celltype
cell_metadata$mg_candidate <- cell_annotation$mg_candidate
write.csv(cell_metadata, file = file.path(output_dir, "cell_metadata_with_tentative_celltypes.csv"), row.names = TRUE)

proj <- addCellColData(
  ArchRProj = proj,
  data = cell_metadata$tentative_celltype,
  cells = rownames(cell_metadata),
  name = "TentativeCelltype",
  force = TRUE
)
proj <- addCellColData(
  ArchRProj = proj,
  data = ifelse(cell_metadata$mg_candidate, "MG_candidate", "Not_MG_candidate"),
  cells = rownames(cell_metadata),
  name = "MGCandidate",
  force = TRUE
)

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

available_rna_features <- get_feature_names(getMatrixFromProject(proj, useMatrix = "GeneExpressionMatrix"))
rna_features_to_plot <- intersect(marker_genes_to_plot, available_rna_features)
if (length(rna_features_to_plot) > 0) {
  z_score_result <- build_rna_expression_z_plots(proj, rna_features_to_plot)
  write.csv(
    z_score_result$z_summary,
    file = file.path(output_dir, "rna_marker_expression_z_summary.csv"),
    row.names = FALSE
  )

  rna_marker_pdf <- file.path(annotation_dir, "rna_marker_umaps.pdf")
  grDevices::pdf(rna_marker_pdf, width = 5, height = 5)
  for (plot_obj in z_score_result$plots) {
    print(plot_obj)
  }
  grDevices::dev.off()
}

p_celltype <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "TentativeCelltype",
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
p_mg_candidate <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "MGCandidate",
  embedding = "UMAP_Combined"
)
plotPDF(
  p_celltype,
  p_mg_candidate,
  p_cluster,
  p_sample,
  name = "tentative_celltypes_clusters_samples.pdf",
  ArchRProj = proj,
  addDOC = FALSE,
  width = 6,
  height = 6
)
file.copy(
  file.path(output_dir, "Plots", "tentative_celltypes_clusters_samples.pdf"),
  annotation_dir,
  overwrite = TRUE
)

celltype_counts <- as.data.frame(table(cell_metadata$tentative_celltype, cell_metadata$Sample))
colnames(celltype_counts) <- c("tentative_celltype", "sample_id", "n_cells")
write.csv(celltype_counts, file = file.path(output_dir, "tentative_celltype_counts_by_sample.csv"), row.names = FALSE)

saveArchRProject(ArchRProj = proj, outputDirectory = output_dir, load = FALSE)

message("Wrote cluster annotations: ", file.path(output_dir, "cluster_tentative_celltype_annotations.csv"))
message("Wrote cell metadata with tentative labels: ", file.path(output_dir, "cell_metadata_with_tentative_celltypes.csv"))
message("Annotation plots copied to: ", annotation_dir)
