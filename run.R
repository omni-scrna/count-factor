#!/usr/bin/env Rscript
# Count factorization for omnibenchmark.
# Code adapted from: https://github.com/omni-scrna/scrapper/blob/main/normalize.R
#
# Output format:
#
# Implementation notes
# --------------------

suppressPackageStartupMessages({
  library(Matrix)
  library(HDF5Array)
  library(BiocSingular)
  library(SingleCellExperiment)
  library(data.table)
  library(anndataR)
  library(NewWave)
  library(glmpca)
})

# arg parsing
source("src/common/r/cli.R")
p <- arg_parser("CNTFCT module")
p <- add_base_args(p)                    # --output_dir, --name
p <- add_stage_args(p, "CNTFCT")     # the stage I/O contract
# your own method params — argparser directly (its add_argument requires `help`):
p <- add_argument(p, "--factorization_type", type = "character", help = "type of factorization")
p <- add_argument(p, "--n_dim", type = "integer", help = "number of latent dimensions")
p <- add_argument(p, "--random_seed", type = "integer", help = "seed")

args <- parse_args(p)                    # argparser's own parser

# logging
cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
cat(sprintf("LOG: command line args\n----------------------------------\n"))
for (i in 1:length(args)) {
  cat(sprintf("  %s: %s\n", names(args)[i], args[[i]]))
}
cat(sprintf("----------------------------------\n"))


run_factorization <- function(sce, args){
  # set seed
  set.seed(args$random_seed)

  cat("class:", paste(class(counts(sce)), collapse = ", "), "\n")

  if (args$factorization_type == "newwave"){
    # counts in memory
    counts(sce) <- as(counts(sce), "CsparseMatrix")
    # select latent dimensions / embeddings
    fitted <- NewWave::newFit(Y = sce, K = args$n_dim, children = 4)
    scores <- NewWave::newW(fitted)
    loadings <- t(NewWave::newAlpha(fitted))
  }

  else if (args$factorization_type == "scgbm"){
    # counts in memory
    Y <- as.matrix(counts(sce))
    out <- scGBM::gbm.sc(Y, M = args$n_dim, ncores = 4, max.iter = 1E5)
    scores <- out$scores
    loadings <- out$loadings
  }

  else if (args$factorization_type == "glmpca"){
    # counts in memory
    Y <- as.matrix(counts(sce))
    out <- glmpca::glmpca(Y, L = args$n_dim, fam = "nb")
    scores <- out$factors
    loadings <- out$loadings
    cat("glmpca components:", paste(names(out), collapse = ", "), "\n")
    cat("dim(scores):", dim(scores), "\n")
    cat("dim(loadings):", dim(loadings), "\n")
  }

  else {
    stop("Unknown factorization_type: ", args$factorization_type)
  }

  list(scores = scores,
       loadings = loadings)
}

main <- function() {
  # read cellids to subset on
  cellids <- readLines(gzfile(args$filtered_cellids))
  cat("length(cellids):", length(cellids), "\n")

  # read H5AD into SCE
  sce <- read_h5ad(args$rawdata_h5ad, as = "SingleCellExperiment")
  sce <- sce[,cellids]

  # read the filtered matrix
  sce_filt <- TENxMatrix(args$normalized_selected_h5, group = "matrix")
  rownames(sce_filt)

  # select the rownames from normalized_selected_h5
  sce <- sce[rownames(sce_filt), ]
  # filter genes with zero counts
  keep <- rowSums(counts(sce)) > 0
  sce  <- sce[keep, ]

  cat("datasets loaded: running factorization\n")
  res <- run_factorization(sce, args)

  # save embeddings: scores
  out_scores_tsv <- file.path(args$output_dir, sprintf("%s_factor_scores.tsv", args$name))
  fwrite(data.frame(cell_id = colnames(sce), res$scores), out_scores_tsv,
      sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("  wrote: %s\n", out_scores_tsv))
  
  # save embeddings: loadings with gene names
  gene_names <- rownames(sce)
  cat(sprintf("[loadings check] gene_names=%d loadings_rows=%d match=%s\n",
    length(gene_names), nrow(res$loadings),
    identical(length(gene_names), nrow(res$loadings))))
  stopifnot(length(gene_names) == nrow(res$loadings))
  
  out_loadings_tsv <- file.path(args$output_dir, sprintf("%s_factor_loadings.tsv", args$name))
  fwrite(data.frame(gene = gene_names, res$loadings), out_loadings_tsv,
      sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("  wrote: %s\n", out_loadings_tsv))
}

if (sys.nframe() == 0L) {
  main()
}