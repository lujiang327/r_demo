#!/usr/bin/env Rscript

message("Checking R packages for ArchR analysis...")

cran_packages <- c("BiocManager", "devtools")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

bioc_packages <- c(
  "BSgenome.Mmusculus.UCSC.mm39",
  "Biostrings",
  "GenomeInfoDb",
  "GenomicFeatures",
  "GenomicRanges",
  "Rsamtools",
  "rhdf5",
  "rtracklayer",
  "SummarizedExperiment"
)

missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

if (!requireNamespace("ArchR", quietly = TRUE)) {
  message("Installing ArchR from GreenleafLab/ArchR...")
  devtools::install_github(
    "GreenleafLab/ArchR",
    repos = BiocManager::repositories(),
    upgrade = "never"
  )
}

message("Done. ArchR available: ", requireNamespace("ArchR", quietly = TRUE))
