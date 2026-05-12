# Day 0 Retina Multiome: ArchR Preprocessing and QC

This repository contains the current ArchR workflow and result previews for two
10x Multiome retina samples collected at Day 0.

### Samples

| Sample | Genotype | Condition | Day | ATAC fragments | RNA feature matrix |
| --- | --- | --- | --- | --- | --- |
| TH1 | Chx10Cre | Control | Day0 | `Sample_15585-TH-1/atac_fragments.tsv.gz` | `Sample_15585-TH-1/filtered_feature_bc_matrix.h5` |
| TH2 | Chx10Cre/Insm1 KO | KO | Day0 | `Sample_15585-TH-2/atac_fragments.tsv.gz` | `Sample_15585-TH-2/filtered_feature_bc_matrix.h5` |

Raw 10x output files are intentionally excluded from git. The sequencing core
aligned both samples with `refdata-cellranger-arc-GRCm39-2024-A`; the project
uses the same GRCm39/mm39 reference for ArchR compatibility.

### Filtering Thresholds

#### Cell Filtering Criteria

Cells are first imported with loose Arrow-file thresholds, then filtered after
the RNA matrix is added so ATAC and RNA QC are applied to the same cells.
After first-pass annotation, old combined clusters `C1` and `C10` were removed
and the ATAC/RNA/combined embeddings and clusters were recomputed.
ArchR renumbers clusters after reclustering, so new labels named `C1` or `C10`
can appear again; those are not the original excluded clusters.

- Arrow creation minimum TSS: `1`
- Arrow creation minimum fragments: `100`
- Combined QC TSS enrichment: `> 10`
- Combined QC ATAC fragments: `> 1000`
- RNA detected genes: `> 1000` and `< 7000`
- RNA UMIs: `> 1500` and `< 30000`
- RNA mitochondrial percent: no cutoff in the current run
- Doublet removal: disabled by default; doublet scores are retained

Full settings:
[script 02](results/r_demo_archr/filter_settings_script_02.csv) and
[script 04](results/r_demo_archr/filter_settings_script_04.csv).

### Filtering Summary

| Sample | Before combined filtering | After combined filtering | Retention |
| --- | ---: | ---: | ---: |
| TH1 | 16,233 | 15,693 | 96.7% |
| TH2 | 16,720 | 16,054 | 96.0% |
| Total | 32,953 | 31,747 | 96.3% |

### Cluster Exclusion Summary

| Sample | Before C1/C10 exclusion | Removed C1/C10 | After reclustering |
| --- | ---: | ---: | ---: |
| TH1 | 15,693 | 92 | 15,601 |
| TH2 | 16,054 | 122 | 15,932 |
| Total | 31,747 | 214 | 31,533 |

Additional count tables:
[cell_counts_after_qc.csv](results/r_demo_archr/cell_counts_after_qc.csv),
[cell_counts_after_rna_combined.csv](results/r_demo_archr/cell_counts_after_rna_combined.csv), and
[cell_counts_summary.csv](results/r_demo_archr/cell_counts_summary.csv).
Cluster exclusion tables:
[removed clusters](results/r_demo_archr/removed_clusters_before_reclustering.csv) and
[post-exclusion counts](results/r_demo_archr/cell_counts_after_cluster_exclusion_reclustered.csv).

## UMAPs

### ATAC

[Full PDF](results/figures/umap/umap_sample_clusters.pdf)

<p>
  <img src="results/readme_figures/atac_umap_sample.png" alt="ATAC UMAP by sample" width="45%">
  <img src="results/readme_figures/atac_umap_clusters.png" alt="ATAC UMAP by cluster" width="45%">
</p>

### RNA / Combined

[Full PDF](results/figures/combined_umap/rna_combined_umap.pdf)

<p>
  <img src="results/readme_figures/rna_umap_sample.png" alt="RNA UMAP by sample" width="33%">
  <img src="results/readme_figures/combined_umap_sample.png" alt="Combined UMAP by sample" width="33%">
  <img src="results/readme_figures/combined_umap_clusters.png" alt="Combined UMAP by cluster" width="33%">
</p>

## QCs

### RNA

[Full PDF](results/figures/cluster_qc/Clusters_RNA_QC.pdf)

<img src="results/figures/cluster_qc/Clusters_RNA_QC.png" alt="RNA cluster QC" width="90%">

### ATAC

[Full PDF](results/figures/cluster_qc/Clusters_ATAC_QC.pdf)

<img src="results/figures/cluster_qc/Clusters_ATAC_QC.png" alt="ATAC cluster QC" width="90%">

### Combined

[Full PDF](results/figures/cluster_qc/Clusters_Combined_QC.pdf)

<img src="results/figures/cluster_qc/Clusters_Combined_QC.png" alt="Combined cluster QC" width="90%">

### Pre-filter QC Views

<p>
  <img src="results/readme_figures/tss_enrichment.png" alt="TSS enrichment" width="45%">
  <img src="results/readme_figures/rna_qc_by_sample.png" alt="RNA QC by sample" width="45%">
</p>

Full PDFs:
[TSS enrichment](results/figures/qc/tss_enrichment.pdf) and
[RNA QC by sample](results/figures/rna_qc/rna_qc_by_sample.pdf).

## Marker Genes

RNA marker UMAPs are plotted from the `GeneExpressionMatrix` using per-gene
expression Z-scores on `UMAP_Combined`. The gallery includes retinal identity
markers, added bipolar/astrocyte/endothelial-pericyte/RPE markers, neurogenic
and stress markers, and `Malat1` as a low-quality diagnostic.

<p>
  <img src="results/readme_figures/markers/Rho.png" alt="Rho" width="33%">
  <img src="results/readme_figures/markers/Gnat1.png" alt="Gnat1" width="33%">
  <img src="results/readme_figures/markers/Nrl.png" alt="Nrl" width="33%">
  <img src="results/readme_figures/markers/Arr3.png" alt="Arr3" width="33%">
  <img src="results/readme_figures/markers/Opn1mw.png" alt="Opn1mw" width="33%">
  <img src="results/readme_figures/markers/Opn1sw.png" alt="Opn1sw" width="33%">
  <img src="results/readme_figures/markers/Vsx1.png" alt="Vsx1" width="33%">
  <img src="results/readme_figures/markers/Vsx2.png" alt="Vsx2" width="33%">
  <img src="results/readme_figures/markers/Car10.png" alt="Car10" width="33%">
  <img src="results/readme_figures/markers/Prkca.png" alt="Prkca" width="33%">
  <img src="results/readme_figures/markers/Sebox.png" alt="Sebox" width="33%">
  <img src="results/readme_figures/markers/Scgn.png" alt="Scgn" width="33%">
  <img src="results/readme_figures/markers/Cabp5.png" alt="Cabp5" width="33%">
  <img src="results/readme_figures/markers/Grm6.png" alt="Grm6" width="33%">
  <img src="results/readme_figures/markers/Tfap2a.png" alt="Tfap2a" width="33%">
  <img src="results/readme_figures/markers/Gad1.png" alt="Gad1" width="33%">
  <img src="results/readme_figures/markers/Gad2.png" alt="Gad2" width="33%">
  <img src="results/readme_figures/markers/Rlbp1.png" alt="Rlbp1" width="33%">
  <img src="results/readme_figures/markers/Glul.png" alt="Glul" width="33%">
  <img src="results/readme_figures/markers/Dkk3.png" alt="Dkk3" width="33%">
  <img src="results/readme_figures/markers/Clu.png" alt="Clu" width="33%">
  <img src="results/readme_figures/markers/Apoe.png" alt="Apoe" width="33%">
  <img src="results/readme_figures/markers/Aqp4.png" alt="Aqp4" width="33%">
  <img src="results/readme_figures/markers/Calb1.png" alt="Calb1" width="33%">
  <img src="results/readme_figures/markers/Onecut1.png" alt="Onecut1" width="33%">
  <img src="results/readme_figures/markers/Rbpms.png" alt="Rbpms" width="33%">
  <img src="results/readme_figures/markers/Pou4f1.png" alt="Pou4f1" width="33%">
  <img src="results/readme_figures/markers/Pou4f2.png" alt="Pou4f2" width="33%">
  <img src="results/readme_figures/markers/Ttr.png" alt="Ttr" width="33%">
  <img src="results/readme_figures/markers/Rdh5.png" alt="Rdh5" width="33%">
  <img src="results/readme_figures/markers/Rpe65.png" alt="Rpe65" width="33%">
  <img src="results/readme_figures/markers/S100b.png" alt="S100b" width="33%">
  <img src="results/readme_figures/markers/Gfap.png" alt="Gfap" width="33%">
  <img src="results/readme_figures/markers/Pax2.png" alt="Pax2" width="33%">
  <img src="results/readme_figures/markers/Pdgfra.png" alt="Pdgfra" width="33%">
  <img src="results/readme_figures/markers/Ctss.png" alt="Ctss" width="33%">
  <img src="results/readme_figures/markers/C1qa.png" alt="C1qa" width="33%">
  <img src="results/readme_figures/markers/Cldn5.png" alt="Cldn5" width="33%">
  <img src="results/readme_figures/markers/Flt1.png" alt="Flt1" width="33%">
  <img src="results/readme_figures/markers/Kcnj8.png" alt="Kcnj8" width="33%">
  <img src="results/readme_figures/markers/Pecam1.png" alt="Pecam1" width="33%">
  <img src="results/readme_figures/markers/Acta2.png" alt="Acta2" width="33%">
  <img src="results/readme_figures/markers/Cdk1.png" alt="Cdk1" width="33%">
  <img src="results/readme_figures/markers/Mki67.png" alt="Mki67" width="33%">
  <img src="results/readme_figures/markers/Top2a.png" alt="Top2a" width="33%">
  <img src="results/readme_figures/markers/Pcna.png" alt="Pcna" width="33%">
  <img src="results/readme_figures/markers/Ascl1.png" alt="Ascl1" width="33%">
  <img src="results/readme_figures/markers/Neurog2.png" alt="Neurog2" width="33%">
  <img src="results/readme_figures/markers/Insm1.png" alt="Insm1" width="33%">
  <img src="results/readme_figures/markers/Atoh7.png" alt="Atoh7" width="33%">
  <img src="results/readme_figures/markers/Neurod1.png" alt="Neurod1" width="33%">
  <img src="results/readme_figures/markers/Otx2.png" alt="Otx2" width="33%">
  <img src="results/readme_figures/markers/Crx.png" alt="Crx" width="33%">
  <img src="results/readme_figures/markers/Olig2.png" alt="Olig2" width="33%">
  <img src="results/readme_figures/markers/Foxn4.png" alt="Foxn4" width="33%">
  <img src="results/readme_figures/markers/Fos.png" alt="Fos" width="33%">
  <img src="results/readme_figures/markers/Jun.png" alt="Jun" width="33%">
  <img src="results/readme_figures/markers/Stat3.png" alt="Stat3" width="33%">
  <img src="results/readme_figures/markers/Lcn2.png" alt="Lcn2" width="33%">
  <img src="results/readme_figures/markers/Malat1.png" alt="Malat1" width="33%">
</p>

Marker tables:
[marker sets](results/r_demo_archr/retinal_marker_sets.csv),
[RNA marker presence](results/r_demo_archr/marker_presence_gene_expression_matrix.csv),
[gene-score marker presence](results/r_demo_archr/marker_presence_gene_score_matrix.csv), and
[ExpressionZ summary](results/r_demo_archr/rna_marker_expression_z_summary.csv).

## Marker Peaks

Peak calling is prepared in
[`scripts/10_call_peaks_add_peak_matrix.R`](scripts/10_call_peaks_add_peak_matrix.R).
The script builds group coverages by `Clusters_Combined`, calls reproducible
peaks with MACS2, adds a peak matrix, and exports the updated project. Complete
this step before interpreting differential accessibility or marker peaks.

## Run Order

Run from the project root with the ArchR conda R:

```sh
conda activate archr
Rscript scripts/00_check_install.R
Rscript scripts/01_create_arrows.R
Rscript scripts/02_build_project_qc.R
Rscript scripts/03_reduce_cluster_umap.R
Rscript scripts/04_add_multiome_rna_combined_embedding.R
Rscript scripts/12_remove_clusters_recluster.R
Rscript scripts/05_summarize_cell_counts.R
Rscript scripts/06_compare_clustering_parameters.R
Rscript scripts/07_annotate_retinal_celltypes.R
Rscript scripts/08_export_readme_marker_pngs.R
Rscript scripts/13_export_readme_umaps.R
Rscript scripts/09_celltype_proportions.R
Rscript scripts/10_call_peaks_add_peak_matrix.R
Rscript scripts/11_cluster_qc_violin_plots.R
```

If VS Code resolves `Rscript` to system R, use the conda Rscript directly:

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
