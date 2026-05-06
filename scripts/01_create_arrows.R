#!/usr/bin/env Rscript

source(file.path("config", "analysis_config.R"))
source(file.path("R", "project_helpers.R"))

suppressPackageStartupMessages(library(ArchR))

set.seed(random_seed)
annotations <- setup_archr_session()

ensure_dirs(arrow_dir, log_dir)

samples <- read_samples(sample_sheet)
validate_input_files(samples)

input_files <- stats::setNames(samples$fragment_path, samples$sample_id)
output_names <- file.path(arrow_dir, names(input_files))

message("Creating Arrow files for ", length(input_files), " samples...")

arrow_files <- createArrowFiles(
  inputFiles = input_files,
  sampleNames = names(input_files),
  outputNames = output_names,
  QCDir = file.path(arrow_dir, "QualityControl"),
  minTSS = min_tss,
  minFrags = min_frags,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  geneAnnotation = annotations$geneAnnotation,
  genomeAnnotation = annotations$genomeAnnotation,
  force = FALSE
)

saveRDS(arrow_files, file = file.path(arrow_dir, "arrow_files.rds"))

if (length(arrow_files) == 0) {
  stop(
    "ArchR returned zero Arrow files. Check the latest ArchRLogs/ArchR-createArrows-*.log file."
  )
}

message("Arrow files written:")
message(paste(" -", arrow_files, collapse = "\n"))
