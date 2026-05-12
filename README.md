# Day 0 Retina Multiome ArchR Analysis

This repository contains the ArchR analysis workflow and current result exports for
two 10x Multiome retina samples collected at Day 0.

## Samples

| Sample | Genotype | Condition | Day | ATAC fragments | RNA feature matrix |
| --- | --- | --- | --- | --- | --- |
| TH1 | Chx10Cre | Control | Day0 | `Sample_15585-TH-1/atac_fragments.tsv.gz` | `Sample_15585-TH-1/filtered_feature_bc_matrix.h5` |
| TH2 | Chx10Cre/Insm1 KO | KO | Day0 | `Sample_15585-TH-2/atac_fragments.tsv.gz` | `Sample_15585-TH-2/filtered_feature_bc_matrix.h5` |

The raw 10x files are intentionally ignored by git. The sequencing core aligned
the samples against `refdata-cellranger-arc-GRCm39-2024-A`; the scripts use the
same GRCm39/mm39 reference for ArchR compatibility.

## Filtering

Arrow files are created with loose thresholds, then stricter combined ATAC/RNA QC
is applied after the RNA matrix is added.

| Step | Threshold |
| --- | --- |
| Arrow creation minimum TSS | `1` |
| Arrow creation minimum fragments | `100` |
| Combined QC TSS enrichment | `> 10` |
| Combined QC ATAC fragments | `> 1000` |
| RNA detected genes | `> 1000` and `< 7000` |
| RNA UMIs | `> 1500` and `< 30000` |
| RNA mitochondrial percent | no cutoff |
| Doublet removal | disabled by default; doublet scores retained |

Full settings:
[script 02](results/r_demo_archr/filter_settings_script_02.csv),
[script 04](results/r_demo_archr/filter_settings_script_04.csv).

## Cell Counts

| Stage | TH1 | TH2 | Total |
| --- | ---: | ---: | ---: |
| After loose Arrow QC, before doublet filtering | 16,233 | 16,720 | 32,953 |
| RNA H5 matched cells before project subset | 16,125 | 16,555 | 32,680 |
| After combined ATAC/RNA filtering | 15,693 | 16,054 | 31,747 |

Full tables:
[cell_counts_after_qc.csv](results/r_demo_archr/cell_counts_after_qc.csv),
[cell_counts_after_rna_combined.csv](results/r_demo_archr/cell_counts_after_rna_combined.csv),
[cell_counts_summary.csv](results/r_demo_archr/cell_counts_summary.csv).

## UMAPs

### ATAC LSI

[Full PDF](results/figures/umap/umap_sample_clusters.pdf)

![ATAC UMAP](results/readme_figures/atac_umap_sample_clusters.png)

### RNA and Combined Embeddings

[Full PDF](results/figures/combined_umap/rna_combined_umap.pdf)

![RNA and combined UMAP](results/readme_figures/rna_combined_umap.png)

## QC

### ATAC TSS Enrichment

[Full PDF](results/figures/qc/tss_enrichment.pdf)

![TSS enrichment](results/readme_figures/tss_enrichment.png)

### RNA QC

[Full PDF](results/figures/rna_qc/rna_qc_by_sample.pdf)

![RNA QC by sample](results/readme_figures/rna_qc_by_sample.png)

## Cell Type Annotation

Tentative labels are marker-score based and should be treated as first-pass
annotations for review. Current marker panels emphasize retinal identity markers
plus activation, proliferation, neurogenic transcription factor, and stress
program markers.

[Full PDF](results/figures/celltype_annotation/tentative_celltypes_clusters_samples.pdf)

![Tentative cell type labels](results/readme_figures/tentative_celltypes_clusters_samples.png)

### Tentative Cell Counts

| Tentative cell type | TH1 | TH2 |
| --- | ---: | ---: |
| AC | 9,046 | 9,074 |
| BC | 6,540 | 6,914 |
| Microglia | 15 | 18 |
| RGC | 92 | 48 |

Full outputs:
[cluster annotations](results/r_demo_archr/cluster_tentative_celltype_annotations.csv),
[cell metadata with labels](results/r_demo_archr/cell_metadata_with_tentative_celltypes.csv),
[cell type counts](results/r_demo_archr/tentative_celltype_counts_by_sample.csv).

## Marker Genes

RNA marker UMAPs are plotted directly from the `GeneExpressionMatrix` using
per-gene expression Z-scores on `UMAP_Combined`.

[Full marker PDF](results/figures/celltype_annotation/rna_marker_umaps.pdf)

### Marker Gallery

<p>
  <img src="results/readme_figures/markers/Rho.png" width="180" alt="Rho">
  <img src="results/readme_figures/markers/Gnat1.png" width="180" alt="Gnat1">
  <img src="results/readme_figures/markers/Nrl.png" width="180" alt="Nrl">
  <img src="results/readme_figures/markers/Arr3.png" width="180" alt="Arr3">
  <img src="results/readme_figures/markers/Opn1mw.png" width="180" alt="Opn1mw">
  <img src="results/readme_figures/markers/Opn1sw.png" width="180" alt="Opn1sw">
  <img src="results/readme_figures/markers/Vsx1.png" width="180" alt="Vsx1">
  <img src="results/readme_figures/markers/Vsx2.png" width="180" alt="Vsx2">
  <img src="results/readme_figures/markers/Car10.png" width="180" alt="Car10">
  <img src="results/readme_figures/markers/Tfap2a.png" width="180" alt="Tfap2a">
  <img src="results/readme_figures/markers/Gad1.png" width="180" alt="Gad1">
  <img src="results/readme_figures/markers/Gad2.png" width="180" alt="Gad2">
  <img src="results/readme_figures/markers/Rlbp1.png" width="180" alt="Rlbp1">
  <img src="results/readme_figures/markers/Glul.png" width="180" alt="Glul">
  <img src="results/readme_figures/markers/Dkk3.png" width="180" alt="Dkk3">
  <img src="results/readme_figures/markers/Clu.png" width="180" alt="Clu">
  <img src="results/readme_figures/markers/Apoe.png" width="180" alt="Apoe">
  <img src="results/readme_figures/markers/Aqp4.png" width="180" alt="Aqp4">
  <img src="results/readme_figures/markers/Calb1.png" width="180" alt="Calb1">
  <img src="results/readme_figures/markers/Onecut1.png" width="180" alt="Onecut1">
  <img src="results/readme_figures/markers/Rbpms.png" width="180" alt="Rbpms">
  <img src="results/readme_figures/markers/Pou4f1.png" width="180" alt="Pou4f1">
  <img src="results/readme_figures/markers/Pou4f2.png" width="180" alt="Pou4f2">
  <img src="results/readme_figures/markers/Ttr.png" width="180" alt="Ttr">
  <img src="results/readme_figures/markers/Rdh5.png" width="180" alt="Rdh5">
  <img src="results/readme_figures/markers/S100b.png" width="180" alt="S100b">
  <img src="results/readme_figures/markers/Gfap.png" width="180" alt="Gfap">
  <img src="results/readme_figures/markers/Ctss.png" width="180" alt="Ctss">
  <img src="results/readme_figures/markers/C1qa.png" width="180" alt="C1qa">
  <img src="results/readme_figures/markers/Cldn5.png" width="180" alt="Cldn5">
  <img src="results/readme_figures/markers/Flt1.png" width="180" alt="Flt1">
  <img src="results/readme_figures/markers/Cdk1.png" width="180" alt="Cdk1">
  <img src="results/readme_figures/markers/Mki67.png" width="180" alt="Mki67">
  <img src="results/readme_figures/markers/Top2a.png" width="180" alt="Top2a">
  <img src="results/readme_figures/markers/Pcna.png" width="180" alt="Pcna">
  <img src="results/readme_figures/markers/Ascl1.png" width="180" alt="Ascl1">
  <img src="results/readme_figures/markers/Neurog2.png" width="180" alt="Neurog2">
  <img src="results/readme_figures/markers/Insm1.png" width="180" alt="Insm1">
  <img src="results/readme_figures/markers/Atoh7.png" width="180" alt="Atoh7">
  <img src="results/readme_figures/markers/Neurod1.png" width="180" alt="Neurod1">
  <img src="results/readme_figures/markers/Otx2.png" width="180" alt="Otx2">
  <img src="results/readme_figures/markers/Crx.png" width="180" alt="Crx">
  <img src="results/readme_figures/markers/Olig2.png" width="180" alt="Olig2">
  <img src="results/readme_figures/markers/Foxn4.png" width="180" alt="Foxn4">
  <img src="results/readme_figures/markers/Fos.png" width="180" alt="Fos">
  <img src="results/readme_figures/markers/Jun.png" width="180" alt="Jun">
  <img src="results/readme_figures/markers/Stat3.png" width="180" alt="Stat3">
  <img src="results/readme_figures/markers/Lcn2.png" width="180" alt="Lcn2">
</p>

Marker tables:
[marker sets](results/r_demo_archr/retinal_marker_sets.csv),
[RNA marker presence](results/r_demo_archr/marker_presence_gene_expression_matrix.csv),
[gene-score marker presence](results/r_demo_archr/marker_presence_gene_score_matrix.csv),
[ExpressionZ summary](results/r_demo_archr/rna_marker_expression_z_summary.csv).
PNG export summary:
[marker_png_summary.csv](results/readme_figures/marker_png_summary.csv).

## Run Order

Run from the project root with the ArchR conda R:

```sh
conda activate archr
Rscript scripts/00_check_install.R
Rscript scripts/01_create_arrows.R
Rscript scripts/02_build_project_qc.R
Rscript scripts/03_reduce_cluster_umap.R
Rscript scripts/04_add_multiome_rna_combined_embedding.R
Rscript scripts/05_summarize_cell_counts.R
Rscript scripts/06_compare_clustering_parameters.R
Rscript scripts/07_annotate_retinal_celltypes.R
```

If VS Code resolves `Rscript` to system R, use:

```sh
/Users/louis/miniforge3/envs/archr/bin/Rscript scripts/07_annotate_retinal_celltypes.R
```

## Reference Setup

Set the 10x ARC reference path before running if it is not in the project root:

```sh
export TENX_ARC_REFERENCE_DIR=/path/to/refdata-cellranger-arc-GRCm39-2024-A
```

Expected files inside the reference include:

- `fasta/genome.fa`
- `fasta/genome.fa.fai`
- `genes/genes.gtf`
