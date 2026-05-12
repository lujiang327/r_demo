project_name <- "r_demo_archr"

# These samples were aligned by Cell Ranger ARC against:
# refdata-cellranger-arc-GRCm39-2024-A
#
# Use "custom_tenx_arc" once `tenx_arc_reference_dir` points to that reference
# folder. Use "builtin" with `builtin_genome <- "mm10"` only as a quick fallback.
genome_mode <- "custom_tenx_arc"
builtin_genome <- "mm10"
bsgenome_package <- "BSgenome.Mmusculus.UCSC.mm39"
tenx_arc_reference_dir <- Sys.getenv(
  "TENX_ARC_REFERENCE_DIR",
  unset = "refdata-cellranger-arc-GRCm39-2024-A"
)
blacklist_path <- NA_character_
chr_prefix <- TRUE
gene_annotation_chromosomes <- paste0("chr", c(1:19, "X", "Y"))
genome_annotation_chromosomes <- gene_annotation_chromosomes

detected_cores <- parallel::detectCores()
if (is.na(detected_cores)) {
  detected_cores <- 8
}
threads <- as.integer(max(1, detected_cores - 1))
random_seed <- 1
archr_locking <- TRUE

sample_sheet <- file.path("config", "samples.tsv")
arrow_dir <- "arrows"
output_dir <- file.path("results", project_name)
figure_dir <- file.path("results", "figures")
log_dir <- "logs"

# Loose Arrow-creation thresholds matching /Users/louis/Desktop/lab/code/atac.
# Stricter combined ATAC/RNA filtering is applied after RNA is added.
min_tss <- 1
min_frags <- 100

# Combined ATAC/RNA filtering thresholds matching /Users/louis/Desktop/lab/code/atac.
combined_filter_min_tss <- 10
combined_filter_min_frags <- 1000

# Doublet detection/removal settings. The reference /atac workflow computes its
# main QC filter without doublet removal, so removal is disabled by default here.
apply_doublet_filter <- FALSE
doublet_filter_ratio <- 1
doublet_filter_cut_enrich <- 1

# RNA QC settings for the Multiome H5 files.
apply_rna_qc_filter <- TRUE
rna_min_genes <- 1000
rna_max_genes <- 7000
rna_min_umi <- 1500
rna_max_umi <- 30000
rna_max_mito_pct <- Inf

# RNA LSI feature-selection settings. `filterQuantile` controls removal of the
# highest-count RNA features before variable-feature selection. A value of 1
# keeps highly expressed genes available for RNA LSI feature selection.
rna_lsi_var_features <- 2500
rna_lsi_filter_quantile <- 1
rna_lsi_total_features <- NULL

# Peak calling settings. MACS2 must be available on PATH, or set MACS2_PATH to
# the executable path before running script 10.
peak_group_by <- "Clusters_Combined"
peak_path_to_macs2 <- Sys.getenv("MACS2_PATH", unset = "")
peak_min_cells_group_coverages <- 40
peak_max_cells_group_coverages <- 500
peak_min_replicates <- 2
peak_max_replicates <- 5
peak_reproducibility <- "2"
peak_genome_size <- "mm"
peak_max_peaks <- 150000
peak_peaks_per_cell <- 500
peak_min_cells <- 25
peak_exclude_chr <- c("chrM", "chrY")
peak_force_group_coverages <- FALSE
peak_force_peak_set <- TRUE
