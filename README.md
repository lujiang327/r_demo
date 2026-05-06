# ArchR Analysis Starter

This project is set up for an ArchR single-cell ATAC analysis using the raw
fragment files already present in this folder.

## Current Inputs

| Sample | Fragment file | RNA / feature matrix |
| --- | --- | --- |
| 15585-TH-1 | `Sample_15585-TH-1/atac_fragments.tsv.gz` | `Sample_15585-TH-1/filtered_feature_bc_matrix.h5` |
| 15585-TH-2 | `Sample_15585-TH-2/atac_fragments.tsv.gz` | `Sample_15585-TH-2/filtered_feature_bc_matrix.h5` |

## Folder Layout

- `config/`: sample manifest and shared analysis settings.
- `R/`: reusable helper functions.
- `scripts/`: numbered scripts to run in order.
- `arrows/`: generated ArchR Arrow files.
- `results/`: ArchR project, plots, and exported tables.
- `logs/`: ArchR log files.

## Run Order

Run these from the project root:

```sh
Rscript scripts/00_check_install.R
Rscript scripts/01_create_arrows.R
Rscript scripts/02_build_project_qc.R
Rscript scripts/03_reduce_cluster_umap.R
Rscript scripts/04_add_multiome_rna_combined_embedding.R
```

The first script checks/install dependencies. The later scripts create Arrow
files, build the ArchR project, generate QC plots, run dimensional reduction,
cluster cells, make ATAC/RNA/combined UMAPs, and add the 10x Multiome RNA matrix.

## Reference Setup

These samples were aligned by the core against the 10x mouse GRCm39 ARC
reference:

`refdata-cellranger-arc-GRCm39-2024-A`

Place or symlink that reference folder in the project root, or update
`TENX_ARC_REFERENCE_DIR` before running the scripts:

```sh
export TENX_ARC_REFERENCE_DIR=/path/to/refdata-cellranger-arc-GRCm39-2024-A
```

Expected files inside the reference include:

- `fasta/genome.fa`
- `fasta/genome.fa.fai`
- `genes/genes.gtf`

The scripts use Bioconductor `BSgenome.Mmusculus.UCSC.mm39` for sequence context,
the 10x FASTA index for chromosome sizes, and the 10x GTF for gene/TSS
annotation.

## Important Notes

- The configured genome is mouse GRCm39/mm39, matching the sequencing core's
  Cell Ranger ARC reference.
- The default QC cutoffs are conservative starter values. Inspect the plots in
  `results/figures/qc/` before deciding whether to tighten them.
- ArchR output can be large. Keep generated files out of cloud sync folders if
  you run into file locking or slow writes.
