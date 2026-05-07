read_samples <- function(path = file.path("config", "samples.tsv")) {
  samples <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)

  required_cols <- c("sample_id", "fragment_path", "feature_matrix_h5")
  missing_cols <- setdiff(required_cols, colnames(samples))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in sample sheet: ", paste(missing_cols, collapse = ", "))
  }

  samples
}

validate_input_files <- function(samples) {
  input_paths <- unique(c(samples$fragment_path, paste0(samples$fragment_path, ".tbi"), samples$feature_matrix_h5))
  missing_paths <- input_paths[!file.exists(input_paths)]

  if (length(missing_paths) > 0) {
    stop("Missing input files:\n", paste(" -", missing_paths, collapse = "\n"))
  }

  invisible(TRUE)
}

ensure_dirs <- function(...) {
  dirs <- unlist(list(...), use.names = FALSE)
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
  }
  invisible(dirs)
}

find_reference_file <- function(reference_dir, candidates, label) {
  paths <- file.path(reference_dir, candidates)
  hit <- paths[file.exists(paths)][1]

  if (is.na(hit)) {
    stop(
      "Could not find ", label, " in reference directory: ", reference_dir, "\n",
      "Looked for:\n", paste(" -", paths, collapse = "\n")
    )
  }

  hit
}

first_existing_metadata_col <- function(gr, candidates) {
  cols <- colnames(S4Vectors::mcols(gr))
  hit <- candidates[candidates %in% cols][1]

  if (is.na(hit)) {
    stop(
      "None of these metadata columns were found in the GTF: ",
      paste(candidates, collapse = ", ")
    )
  }

  as.character(S4Vectors::mcols(gr)[[hit]])
}

load_bsgenome <- function(package_name) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    stop(
      "Missing BSgenome package: ", package_name, "\n",
      "Install it inside the active conda env with:\n",
      "Rscript -e 'BiocManager::install(\"", package_name, "\", ask = FALSE, update = FALSE)'"
    )
  }

  getExportedValue(package_name, package_name)
}

build_tenx_arc_genome_annotation <- function(
  reference_dir,
  bsgenome_package,
  keep_chromosomes = NULL,
  blacklist_path = NA_character_
) {
  fasta_path <- find_reference_file(
    reference_dir = reference_dir,
    candidates = c("fasta/genome.fa", "fasta/genome.fasta", "genome.fa", "genome.fasta"),
    label = "genome FASTA"
  )

  if (!file.exists(paste0(fasta_path, ".fai"))) {
    stop(
      "Missing FASTA index: ", paste0(fasta_path, ".fai"), "\n",
      "Create it with: samtools faidx ", fasta_path
    )
  }

  chrom_sizes <- Rsamtools::scanFaIndex(fasta_path)
  chrom_sizes <- filter_to_chromosomes(chrom_sizes, keep_chromosomes, "genome annotation")
  bsgenome <- load_bsgenome(bsgenome_package)

  blacklist <- GenomicRanges::GRanges()
  if (!is.na(blacklist_path) && nzchar(blacklist_path)) {
    if (!file.exists(blacklist_path)) {
      stop("Configured blacklist does not exist: ", blacklist_path)
    }
    blacklist <- rtracklayer::import(blacklist_path)
  }

  ArchR::createGenomeAnnotation(
    genome = bsgenome,
    chromSizes = chrom_sizes,
    blacklist = blacklist,
    filter = TRUE,
    filterChr = NULL
  )
}

filter_to_chromosomes <- function(gr, chromosomes, label) {
  if (is.null(chromosomes) || length(chromosomes) == 0) {
    return(gr)
  }

  chromosomes <- as.character(chromosomes)
  missing_chromosomes <- setdiff(chromosomes, seqlevels(gr))
  if (length(missing_chromosomes) > 0) {
    stop(
      "Configured ", label, " chromosomes are missing from the annotation: ",
      paste(missing_chromosomes, collapse = ", ")
    )
  }

  GenomeInfoDb::keepSeqlevels(gr, chromosomes, pruning.mode = "coarse")
}

build_tenx_arc_gene_annotation <- function(reference_dir, keep_chromosomes = NULL) {
  gtf_path <- find_reference_file(
    reference_dir = reference_dir,
    candidates = c("genes/genes.gtf", "genes.gtf"),
    label = "gene annotation GTF"
  )

  gtf <- rtracklayer::import(gtf_path)
  genes <- gtf[gtf$type == "gene"]
  exons <- gtf[gtf$type == "exon"]

  genes <- filter_to_chromosomes(genes, keep_chromosomes, "gene annotation")
  exons <- filter_to_chromosomes(exons, keep_chromosomes, "exon annotation")

  if (length(genes) == 0 || length(exons) == 0) {
    stop("The GTF must contain both gene and exon records: ", gtf_path)
  }

  gene_symbols <- first_existing_metadata_col(genes, c("gene_name", "gene", "Name", "gene_id"))
  exon_symbols <- first_existing_metadata_col(exons, c("gene_name", "gene", "Name", "gene_id"))

  S4Vectors::mcols(genes)$symbol <- gene_symbols
  S4Vectors::mcols(genes)$symbols <- gene_symbols
  S4Vectors::mcols(exons)$symbol <- exon_symbols
  S4Vectors::mcols(exons)$symbols <- exon_symbols

  tss <- GenomicRanges::promoters(genes, upstream = 0, downstream = 1)
  S4Vectors::mcols(tss)$symbol <- gene_symbols
  S4Vectors::mcols(tss)$symbols <- gene_symbols

  ArchR::createGeneAnnotation(
    genes = genes,
    exons = exons,
    TSS = tss
  )
}

setup_archr_session <- function() {
  addArchRLocking(locking = archr_locking)
  addArchRThreads(threads = threads)

  if (identical(genome_mode, "builtin")) {
    addArchRGenome(builtin_genome)
    return(list(genomeAnnotation = NULL, geneAnnotation = NULL))
  }

  if (!identical(genome_mode, "custom_tenx_arc")) {
    stop("Unsupported genome_mode: ", genome_mode)
  }

  if (!dir.exists(tenx_arc_reference_dir)) {
    stop(
      "10x ARC reference directory not found: ", tenx_arc_reference_dir, "\n",
      "Place or symlink refdata-cellranger-arc-GRCm39-2024-A here, or update ",
      "`tenx_arc_reference_dir` in config/analysis_config.R."
    )
  }

  addArchRChrPrefix(chrPrefix = chr_prefix)

  list(
    genomeAnnotation = build_tenx_arc_genome_annotation(
      reference_dir = tenx_arc_reference_dir,
      bsgenome_package = bsgenome_package,
      keep_chromosomes = genome_annotation_chromosomes,
      blacklist_path = blacklist_path
    ),
    geneAnnotation = build_tenx_arc_gene_annotation(
      reference_dir = tenx_arc_reference_dir,
      keep_chromosomes = gene_annotation_chromosomes
    )
  )
}

count_cells_by_sample <- function(cell_metadata, stage) {
  counts <- as.data.frame(table(cell_metadata$Sample), stringsAsFactors = FALSE)
  colnames(counts) <- c("sample_id", "n_cells")
  counts$stage <- stage
  counts[, c("stage", "sample_id", "n_cells")]
}

write_cell_count_summary <- function(cell_metadata, path, stage) {
  counts <- count_cells_by_sample(cell_metadata, stage)
  counts <- rbind(
    counts,
    data.frame(stage = stage, sample_id = "Total", n_cells = sum(counts$n_cells))
  )
  write.csv(counts, file = path, row.names = FALSE)
  invisible(counts)
}

write_filter_settings <- function(path, settings) {
  settings_df <- data.frame(
    setting = names(settings),
    value = unlist(settings, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  write.csv(settings_df, file = path, row.names = FALSE)
  invisible(settings_df)
}

count_cells_in_arrow_files <- function(arrow_files, stage = "after_arrow_qc_before_doublet_filter") {
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop("Package `rhdf5` is required to count cells in Arrow files.")
  }

  counts <- lapply(arrow_files, function(arrow_file) {
    sample_id <- rhdf5::h5read(arrow_file, "Metadata/Sample")
    cell_names <- rhdf5::h5read(arrow_file, "Metadata/CellNames")
    data.frame(
      stage = stage,
      sample_id = as.character(sample_id[1]),
      n_cells = length(cell_names),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, counts)
}

count_unique_fragment_barcodes <- function(samples, stage = "raw_fragment_unique_barcodes_before_archr_qc") {
  counts <- lapply(seq_len(nrow(samples)), function(i) {
    fragment_path <- samples$fragment_path[i]
    sample_id <- samples$sample_id[i]

    if (!file.exists(fragment_path)) {
      stop("Fragment file does not exist: ", fragment_path)
    }

    con <- gzfile(fragment_path, open = "rt")
    on.exit(close(con), add = TRUE)

    barcodes <- new.env(hash = TRUE, parent = emptyenv())
    repeat {
      lines <- readLines(con, n = 100000L)
      if (length(lines) == 0) {
        break
      }

      fields <- strsplit(lines, "\t", fixed = TRUE)
      sample_barcodes <- vapply(fields, `[`, character(1), 4)
      for (barcode in sample_barcodes) {
        assign(barcode, TRUE, envir = barcodes)
      }
    }

    data.frame(
      stage = stage,
      sample_id = sample_id,
      n_cells = length(ls(barcodes, all.names = TRUE)),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, counts)
}
