#!/usr/bin/env Rscript

source(file.path("config", "analysis_config.R"))
source(file.path("R", "project_helpers.R"))

suppressPackageStartupMessages({
  library(ArchR)
  library(GenomicRanges)
  library(bsgenome_package, character.only = TRUE)
})

set.seed(random_seed)
setup_archr_session()

peak_dir <- file.path(figure_dir, "peaks")
ensure_dirs(output_dir, peak_dir, log_dir)

resolve_macs2_path <- function() {
  if (nzchar(peak_path_to_macs2)) {
    if (!file.exists(peak_path_to_macs2)) {
      stop("Configured MACS2_PATH does not exist: ", peak_path_to_macs2)
    }
    return(peak_path_to_macs2)
  }

  candidate_paths <- c(
    Sys.which("macs2"),
    file.path(Sys.getenv("HOME"), "miniforge3", "envs", "macs2", "bin", "macs2"),
    file.path(Sys.getenv("HOME"), "miniconda3", "envs", "macs2", "bin", "macs2"),
    file.path(Sys.getenv("HOME"), "anaconda3", "envs", "macs2", "bin", "macs2")
  )
  candidate_paths <- candidate_paths[nzchar(candidate_paths)]
  candidate_paths <- candidate_paths[file.exists(candidate_paths)]
  if (length(candidate_paths) > 0) {
    return(candidate_paths[[1]])
  }

  tryCatch(
    findMacs2(),
    error = function(e) {
      stop(
        "MACS2 was not found. Install it in the ArchR environment, for example:\n",
        "  conda install -n archr -c bioconda -c conda-forge macs2\n",
        "or set MACS2_PATH to the executable before running this script.\n",
        "Original ArchR error: ", conditionMessage(e)
      )
    }
  )
}

proj <- loadArchRProject(path = output_dir)

cell_col_data <- as.data.frame(getCellColData(proj))
if (!peak_group_by %in% colnames(cell_col_data)) {
  stop(
    "Missing peak_group_by column `", peak_group_by, "` in ArchR cell metadata. ",
    "Run scripts 03 and 04 before script 10."
  )
}

group_counts <- as.data.frame(table(cell_col_data[[peak_group_by]]), stringsAsFactors = FALSE)
colnames(group_counts) <- c("group", "n_cells")
group_counts <- group_counts[order(group_counts$group), , drop = FALSE]
write.csv(
  group_counts,
  file = file.path(output_dir, paste0("peak_calling_group_counts_", peak_group_by, ".csv")),
  row.names = FALSE
)

path_to_macs2 <- resolve_macs2_path()
message("Using MACS2 at: ", path_to_macs2)
message("Calling peaks grouped by: ", peak_group_by)

proj <- addGroupCoverages(
  ArchRProj = proj,
  groupBy = peak_group_by,
  minCells = peak_min_cells_group_coverages,
  maxCells = peak_max_cells_group_coverages,
  minReplicates = peak_min_replicates,
  maxReplicates = peak_max_replicates,
  force = peak_force_group_coverages
)

proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = peak_group_by,
  pathToMacs2 = path_to_macs2,
  genomeSize = peak_genome_size,
  reproducibility = peak_reproducibility,
  peaksPerCell = peak_peaks_per_cell,
  maxPeaks = peak_max_peaks,
  minCells = peak_min_cells,
  excludeChr = peak_exclude_chr,
  force = peak_force_peak_set
)

proj <- addPeakMatrix(
  ArchRProj = proj,
  force = TRUE
)

peak_set <- getPeakSet(proj)
peak_df <- as.data.frame(peak_set)
peak_df$peak_id <- paste0(seqnames(peak_set), ":", start(peak_set), "-", end(peak_set))
peak_df <- peak_df[, c("peak_id", setdiff(colnames(peak_df), "peak_id")), drop = FALSE]

write.csv(
  peak_df,
  file = file.path(output_dir, paste0("peak_set_", peak_group_by, ".csv")),
  row.names = FALSE
)

bed_df <- data.frame(
  chrom = as.character(seqnames(peak_set)),
  start = start(peak_set) - 1,
  end = end(peak_set),
  name = peak_df$peak_id,
  score = 0,
  strand = ".",
  stringsAsFactors = FALSE
)
write.table(
  bed_df,
  file = file.path(output_dir, paste0("peak_set_", peak_group_by, ".bed")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

peak_summary <- data.frame(
  peak_group_by = peak_group_by,
  n_peaks = length(peak_set),
  n_cells = nrow(cell_col_data),
  macs2_path = path_to_macs2,
  min_cells_group_coverages = peak_min_cells_group_coverages,
  max_cells_group_coverages = peak_max_cells_group_coverages,
  min_replicates = peak_min_replicates,
  max_replicates = peak_max_replicates,
  reproducibility = peak_reproducibility,
  genome_size = peak_genome_size,
  peaks_per_cell = peak_peaks_per_cell,
  max_peaks = peak_max_peaks,
  min_cells_peak_calling = peak_min_cells,
  excluded_chromosomes = paste(peak_exclude_chr, collapse = ";"),
  stringsAsFactors = FALSE
)
write.csv(
  peak_summary,
  file = file.path(output_dir, paste0("peak_calling_summary_", peak_group_by, ".csv")),
  row.names = FALSE
)

saveArchRProject(ArchRProj = proj, outputDirectory = output_dir, load = FALSE)

message("Peak calling complete.")
message("Peak set CSV: ", file.path(output_dir, paste0("peak_set_", peak_group_by, ".csv")))
message("Peak set BED: ", file.path(output_dir, paste0("peak_set_", peak_group_by, ".bed")))
