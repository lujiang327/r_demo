# ArchR Walkthrough

## 1. Check Your Inputs

Each sample needs:

- `atac_fragments.tsv.gz`
- `atac_fragments.tsv.gz.tbi`
- optionally, `filtered_feature_bc_matrix.h5` for later RNA / gene expression work

The sample table lives at `config/samples.tsv`. Edit that file whenever you add
or rename samples.

## 2. Configure the Analysis

Open `config/analysis_config.R`.

Key settings:

- `genome_mode`: keep this as `custom_tenx_arc` for the Cell Ranger ARC GRCm39 reference.
- `TENX_ARC_REFERENCE_DIR`: optional environment variable pointing to `refdata-cellranger-arc-GRCm39-2024-A`.
- `bsgenome_package`: keep this as `BSgenome.Mmusculus.UCSC.mm39`.
- `min_tss`: minimum TSS enrichment during Arrow creation.
- `min_frags`: minimum unique fragments during Arrow creation.
- `threads`: number of CPU threads.

The reference folder should contain:

- `fasta/genome.fa`
- `fasta/genome.fa.fai`
- `genes/genes.gtf`

If the FASTA index is missing, create it before running ArchR:

```sh
samtools faidx refdata-cellranger-arc-GRCm39-2024-A/fasta/genome.fa
```

If the reference is not inside the project root, set:

```sh
export TENX_ARC_REFERENCE_DIR=/path/to/refdata-cellranger-arc-GRCm39-2024-A
```

## 3. Install ArchR

Run:

```sh
conda activate archr
Rscript scripts/00_check_install.R
```

This installs missing dependencies and ArchR. It can take a while on a fresh R
setup.

## 4. Create Arrow Files

Run:

```sh
Rscript scripts/01_create_arrows.R
```

Arrow files are ArchR's on-disk format. This step reads the fragment files,
computes QC metrics, builds the tile matrix, and builds gene score matrices.

## 5. Build the Project and QC

Run:

```sh
Rscript scripts/02_build_project_qc.R
```

This creates the ArchR project, adds metadata from `config/samples.tsv`, computes
doublet scores, filters doublets, and writes QC plots.

Look at:

- `results/figures/qc/tss_enrichment.pdf`
- `results/figures/qc/fragment_sizes.pdf`
- `results/r_demo_archr/cell_metadata_after_qc.csv`

## 6. LSI, Clustering, and UMAP

Run:

```sh
Rscript scripts/03_reduce_cluster_umap.R
```

This runs the standard ArchR ATAC workflow:

- iterative LSI
- Seurat clustering
- UMAP

Look at:

- `results/figures/umap/umap_sample_clusters.pdf`
- `results/r_demo_archr/cell_metadata_after_clustering.csv`

## 7. What to Do Next

## 7. Add Multiome RNA and Combined Embedding

Run:

```sh
Rscript scripts/04_add_multiome_rna_combined_embedding.R
```

This imports each `filtered_feature_bc_matrix.h5`, adds the RNA counts as
`GeneExpressionMatrix`, creates RNA QC plots, and computes:

- `LSI_RNA`
- `LSI_Combined`
- `UMAP_RNA`
- `UMAP_Combined`
- `Clusters_RNA`
- `Clusters_Combined`

Look at:

- `results/figures/rna_qc/rna_qc_by_sample.pdf`
- `results/figures/combined_umap/rna_combined_umap.pdf`
- `results/r_demo_archr/cell_metadata_after_rna_combined.csv`

## 8. What to Do Next

After the first UMAP looks reasonable, common next scripts are:

- marker genes / marker peaks by cluster
- peak calling with MACS2
- motif enrichment
- cell type annotation
- integration with the RNA counts in `filtered_feature_bc_matrix.h5`

Do not start those until QC and clustering are sensible.
