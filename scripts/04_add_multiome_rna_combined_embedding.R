#!/usr/bin/env Rscript

source(file.path("config", "analysis_config.R"))
source(file.path("R", "project_helpers.R"))

suppressPackageStartupMessages({
  library(ArchR)
  library(SummarizedExperiment)
})

set.seed(random_seed)
setup_archr_session()

ensure_dirs(file.path(figure_dir, "rna_qc"), file.path(figure_dir, "combined_umap"), log_dir)

samples <- read_samples(sample_sheet)
validate_input_files(samples)

proj <- loadArchRProject(path = output_dir)
project_cells <- rownames(getCellColData(proj))
cell_metadata_before_rna <- as.data.frame(getCellColData(proj))
project_chromosomes <- as.character(GenomicRanges::seqnames(proj@genomeAnnotation$chromSizes))

prefix_rna_cell_names <- function(se_rna, sample_id) {
  cell_names <- colnames(se_rna)
  expected_prefix <- paste0(sample_id, "#")

  if (!all(startsWith(cell_names, expected_prefix))) {
    cell_names <- sub("^.*#", "", cell_names)
    colnames(se_rna) <- paste0(expected_prefix, cell_names)
  }

  se_rna
}

filter_rna_features_to_project_chromosomes <- function(se_rna, project_chromosomes) {
  rr <- rowRanges(se_rna)

  if (!inherits(rr, "GRanges")) {
    return(se_rna)
  }

  keep <- !is.na(GenomicRanges::seqnames(rr)) &
    as.character(GenomicRanges::seqnames(rr)) %in% project_chromosomes

  se_rna[keep, ]
}

import_sample_rna <- function(sample_id, h5_path, project_cells, project_chromosomes) {
  message("Importing RNA H5 for ", sample_id, ": ", h5_path)

  se_rna <- import10xFeatureMatrix(input = h5_path, names = sample_id)
  se_rna <- prefix_rna_cell_names(se_rna, sample_id)
  se_rna <- filter_rna_features_to_project_chromosomes(se_rna, project_chromosomes)

  matched_cells <- intersect(project_cells, colnames(se_rna))
  if (length(matched_cells) == 0) {
    stop("No RNA cells matched ArchR cells for sample: ", sample_id)
  }

  se_rna[, matched_cells]
}

rna_list <- lapply(seq_len(nrow(samples)), function(i) {
  import_sample_rna(
    sample_id = samples$sample_id[i],
    h5_path = samples$feature_matrix_h5[i],
    project_cells = project_cells,
    project_chromosomes = project_chromosomes
  )
})
names(rna_list) <- samples$sample_id
rna_matched_counts <- data.frame(
  stage = "matched_cells_in_rna_h5_before_project_subset",
  sample_id = names(rna_list),
  n_cells = vapply(rna_list, ncol, integer(1)),
  row.names = NULL
)

common_features <- Reduce(intersect, lapply(rna_list, rownames))
if (length(common_features) == 0) {
  stop("No shared RNA features were found across samples.")
}

rna_list <- lapply(rna_list, function(se_rna) se_rna[common_features, ])
rna_counts <- do.call(cbind, lapply(rna_list, SummarizedExperiment::assay))
se_rna_all <- SummarizedExperiment(
  assays = list(counts = rna_counts),
  rowRanges = rowRanges(rna_list[[1]])
)

matched_project_cells <- intersect(project_cells, colnames(se_rna_all))
if (length(matched_project_cells) < length(project_cells)) {
  message(
    "Subsetting ArchR project from ", length(project_cells), " to ",
    length(matched_project_cells), " cells with matched RNA."
  )
  proj <- proj[matched_project_cells]
  se_rna_all <- se_rna_all[, matched_project_cells]
}

proj <- addGeneExpressionMatrix(
  input = proj,
  seRNA = se_rna_all,
  strictMatch = TRUE,
  force = TRUE
)

p_rna_genes <- plotGroups(
  ArchRProj = proj,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "Gex_nGenes",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

p_rna_umi <- plotGroups(
  ArchRProj = proj,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "Gex_nUMI",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

plotPDF(
  p_rna_genes,
  p_rna_umi,
  name = "rna_qc_by_sample.pdf",
  ArchRProj = proj,
  addDOC = FALSE,
  width = 6,
  height = 6
)
file.copy(file.path(output_dir, "Plots", "rna_qc_by_sample.pdf"), file.path(figure_dir, "rna_qc"), overwrite = TRUE)

if (isTRUE(apply_rna_qc_filter)) {
  mito_pct <- proj$Gex_MitoRatio
  if (max(mito_pct, na.rm = TRUE) <= 1) {
    mito_pct <- mito_pct * 100
  }

  keep_cells <- proj$Gex_nGenes >= rna_min_genes &
    proj$Gex_nGenes <= rna_max_genes &
    proj$Gex_nUMI >= rna_min_umi &
    proj$Gex_nUMI <= rna_max_umi &
    mito_pct <= rna_max_mito_pct

  message("RNA QC filter keeping ", sum(keep_cells), " of ", length(keep_cells), " cells.")
  proj <- proj[keep_cells]
}

lsi_rna_args <- list(
  ArchRProj = proj,
  useMatrix = "GeneExpressionMatrix",
  depthCol = "Gex_nUMI",
  name = "LSI_RNA",
  iterations = 2,
  clusterParams = list(
    resolution = c(0.2),
    sampleCells = 10000,
    n.start = 10
  ),
  varFeatures = rna_lsi_var_features,
  firstSelection = "variable",
  filterQuantile = rna_lsi_filter_quantile,
  binarize = FALSE,
  dimsToUse = 1:30,
  force = TRUE
)

if (!is.null(rna_lsi_total_features)) {
  lsi_rna_args$totalFeatures <- as.integer(rna_lsi_total_features)
}

proj <- do.call(addIterativeLSI, lsi_rna_args)

proj <- addCombinedDims(
  ArchRProj = proj,
  reducedDims = c("IterativeLSI", "LSI_RNA"),
  name = "LSI_Combined"
)

proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "LSI_RNA",
  name = "UMAP_RNA",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine",
  force = TRUE
)

proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "LSI_Combined",
  name = "UMAP_Combined",
  nNeighbors = 30,
  minDist = 0.4,
  metric = "cosine",
  force = TRUE
)

proj <- addClusters(
  input = proj,
  reducedDims = "LSI_RNA",
  method = "Seurat",
  name = "Clusters_RNA",
  resolution = 0.8,
  force = TRUE
)

proj <- addClusters(
  input = proj,
  reducedDims = "LSI_Combined",
  method = "Seurat",
  name = "Clusters_Combined",
  resolution = 1.0,
  force = TRUE
)

p_rna_sample <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "UMAP_RNA"
)

p_combined_sample <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "UMAP_Combined"
)

p_combined_cluster <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "Clusters_Combined",
  embedding = "UMAP_Combined"
)

plotPDF(
  p_rna_sample,
  p_combined_sample,
  p_combined_cluster,
  name = "rna_combined_umap.pdf",
  ArchRProj = proj,
  addDOC = FALSE,
  width = 6,
  height = 6
)
file.copy(file.path(output_dir, "Plots", "rna_combined_umap.pdf"), file.path(figure_dir, "combined_umap"), overwrite = TRUE)

cell_col_data <- as.data.frame(getCellColData(proj))
write.csv(cell_col_data, file = file.path(output_dir, "cell_metadata_after_rna_combined.csv"), row.names = TRUE)
cell_counts_summary <- rbind(
  count_cells_by_sample(cell_metadata_before_rna, "before_rna_matching"),
  rna_matched_counts,
  count_cells_by_sample(cell_col_data, "after_rna_combined")
)
cell_counts_summary <- rbind(
  cell_counts_summary,
  aggregate(n_cells ~ stage, data = cell_counts_summary, sum) |>
    transform(sample_id = "Total") |>
    subset(select = c("stage", "sample_id", "n_cells"))
)
write.csv(
  cell_counts_summary[order(cell_counts_summary$stage, cell_counts_summary$sample_id), ],
  file = file.path(output_dir, "cell_counts_after_rna_combined.csv"),
  row.names = FALSE
)
write_filter_settings(
  path = file.path(output_dir, "filter_settings_script_04.csv"),
  settings = list(
    strict_match_rna_to_atac_cells = TRUE,
    apply_rna_qc_filter = apply_rna_qc_filter,
    rna_min_genes = rna_min_genes,
    rna_max_genes = rna_max_genes,
    rna_min_umi = rna_min_umi,
    rna_max_umi = rna_max_umi,
    rna_max_mito_pct = rna_max_mito_pct,
    rna_lsi_var_features = rna_lsi_var_features,
    rna_lsi_filter_quantile = rna_lsi_filter_quantile,
    rna_lsi_total_features = ifelse(is.null(rna_lsi_total_features), "NULL", rna_lsi_total_features),
    gene_expression_matrix_added = TRUE,
    note = "RNA QC thresholds are only applied when apply_rna_qc_filter is TRUE. Cells are always subset to those present in both ArchR and RNA H5 matrices."
  )
)

saveArchRProject(ArchRProj = proj, outputDirectory = output_dir, load = FALSE)

message("RNA integration and combined embedding saved in: ", output_dir)
message("RNA QC plots copied to: ", file.path(figure_dir, "rna_qc"))
message("Combined UMAP plots copied to: ", file.path(figure_dir, "combined_umap"))
