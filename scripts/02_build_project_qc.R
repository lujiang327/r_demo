#!/usr/bin/env Rscript

source(file.path("config", "analysis_config.R"))
source(file.path("R", "project_helpers.R"))

suppressPackageStartupMessages(library(ArchR))

set.seed(random_seed)
annotations <- setup_archr_session()

ensure_dirs(output_dir, figure_dir, file.path(figure_dir, "qc"), log_dir)

arrow_files <- readRDS(file.path(arrow_dir, "arrow_files.rds"))
samples <- read_samples(sample_sheet)

proj <- ArchRProject(
  ArrowFiles = arrow_files,
  outputDirectory = output_dir,
  copyArrows = TRUE,
  geneAnnotation = annotations$geneAnnotation,
  genomeAnnotation = annotations$genomeAnnotation
)

sample_metadata <- samples[, setdiff(colnames(samples), c("fragment_path", "feature_matrix_h5")), drop = FALSE]
rownames(sample_metadata) <- sample_metadata$sample_id
cell_names <- rownames(getCellColData(proj))
sample_ids <- proj$Sample

for (col_name in setdiff(colnames(sample_metadata), "sample_id")) {
  values <- sample_metadata[sample_ids, col_name]
  proj <- addCellColData(
    ArchRProj = proj,
    data = values,
    cells = cell_names,
    name = col_name,
    force = TRUE
  )
}

proj <- addDoubletScores(
  input = proj,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1
)

proj <- filterDoublets(
  proj,
  filterRatio = doublet_filter_ratio,
  cutEnrich = doublet_filter_cut_enrich
)

qc_dir <- file.path(figure_dir, "qc")

p_tss <- plotTSSEnrichment(proj)
plotPDF(p_tss, name = "tss_enrichment.pdf", ArchRProj = proj, addDOC = FALSE, width = 6, height = 6)
file.copy(file.path(output_dir, "Plots", "tss_enrichment.pdf"), qc_dir, overwrite = TRUE)

p_frags <- plotFragmentSizes(proj)
plotPDF(p_frags, name = "fragment_sizes.pdf", ArchRProj = proj, addDOC = FALSE, width = 6, height = 6)
file.copy(file.path(output_dir, "Plots", "fragment_sizes.pdf"), qc_dir, overwrite = TRUE)

cell_col_data <- as.data.frame(getCellColData(proj))
write.csv(cell_col_data, file = file.path(output_dir, "cell_metadata_after_qc.csv"), row.names = TRUE)

saveArchRProject(ArchRProj = proj, outputDirectory = output_dir, load = FALSE)

message("Project saved to: ", output_dir)
message("QC plots copied to: ", qc_dir)
