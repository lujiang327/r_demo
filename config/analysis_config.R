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

threads <- max(1, parallel::detectCores() - 1)
random_seed <- 1
archr_locking <- TRUE

sample_sheet <- file.path("config", "samples.tsv")
arrow_dir <- "arrows"
output_dir <- file.path("results", project_name)
figure_dir <- file.path("results", "figures")
log_dir <- "logs"

# Starter QC thresholds. Review QC plots before treating these as final.
min_tss <- 4
min_frags <- 1000

# Doublet detection/removal settings.
doublet_filter_ratio <- 1
doublet_filter_cut_enrich <- 1

# RNA QC settings for the Multiome H5 files. Keep filtering off until the RNA QC
# plots have been inspected.
apply_rna_qc_filter <- TRUE
rna_min_genes <- 800
rna_max_genes <- 8000
rna_min_umi <- 1200
rna_max_umi <- 30000
rna_max_mito_pct <- 10

# RNA LSI feature-selection settings. `filterQuantile` controls removal of the
# highest-count RNA features before variable-feature selection. A value of 1
# keeps highly expressed genes available for RNA LSI feature selection.
rna_lsi_var_features <- 2500
rna_lsi_filter_quantile <- 1
rna_lsi_total_features <- NULL
