#!/usr/bin/env Rscript

source(file.path("config", "analysis_config.R"))
source(file.path("R", "project_helpers.R"))

samples <- read_samples(sample_sheet)

metadata_files <- c(
  after_qc_doublet_filter = file.path(output_dir, "cell_metadata_after_qc.csv"),
  after_atac_clustering = file.path(output_dir, "cell_metadata_after_clustering.csv"),
  after_rna_combined = file.path(output_dir, "cell_metadata_after_rna_combined.csv")
)
arrow_files <- c(
  list.files(arrow_dir, pattern = "[.](a|A)rrow$", full.names = TRUE),
  list.files(file.path(output_dir, "ArrowFiles"), pattern = "[.](a|A)rrow$", full.names = TRUE)
)
arrow_files <- unique(arrow_files[file.exists(arrow_files)])

existing_files <- metadata_files[file.exists(metadata_files)]
if (length(existing_files) == 0) {
  stop("No cell metadata files found in: ", output_dir)
}

metadata_counts <- do.call(
  rbind,
  lapply(names(existing_files), function(stage) {
    metadata <- read.csv(existing_files[[stage]], check.names = FALSE)
    count_cells_by_sample(metadata, stage)
  })
)

counts <- count_unique_fragment_barcodes(samples)

if (length(arrow_files) > 0) {
  counts <- rbind(
    counts,
    count_cells_in_arrow_files(arrow_files[!duplicated(basename(arrow_files))]),
    metadata_counts
  )
} else {
  counts <- rbind(counts, metadata_counts)
}

counts <- rbind(
  counts,
  aggregate(n_cells ~ stage, data = counts, sum) |>
    transform(sample_id = "Total") |>
    subset(select = c("stage", "sample_id", "n_cells"))
)

counts <- counts[order(counts$stage, counts$sample_id), ]
out_path <- file.path(output_dir, "cell_counts_summary.csv")
write.csv(counts, file = out_path, row.names = FALSE)

settings <- list(
  genome_mode = genome_mode,
  builtin_genome = builtin_genome,
  bsgenome_package = bsgenome_package,
  tenx_arc_reference_dir = tenx_arc_reference_dir,
  gene_annotation_chromosomes = paste(gene_annotation_chromosomes, collapse = ";"),
  genome_annotation_chromosomes = paste(genome_annotation_chromosomes, collapse = ";"),
  min_tss_arrow_creation = min_tss,
  min_frags_arrow_creation = min_frags,
  doublet_filter_ratio = doublet_filter_ratio,
  doublet_filter_cut_enrich = doublet_filter_cut_enrich,
  apply_rna_qc_filter = apply_rna_qc_filter,
  rna_min_genes = rna_min_genes,
  rna_max_genes = rna_max_genes,
  rna_min_umi = rna_min_umi,
  rna_max_umi = rna_max_umi,
  raw_fragment_unique_barcodes_note = "Counts unique barcodes in the raw fragment file before ArchR minTSS/minFrags filtering; this includes low-quality and non-cell barcodes."
)
settings_path <- file.path(output_dir, "analysis_filter_settings.csv")
write_filter_settings(settings_path, settings)

message("Wrote cell count summary: ", out_path)
message("Wrote analysis/filter settings: ", settings_path)
print(counts)
