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
  library(NewWave)
  library(anndataR)
  library(SingleCellExperiment)
  library(data.table)
})

# arg parsing
source("src/common/r/cli.R")
p <- arg_parser("CNTFCT module")
p <- add_base_args(p)                    # --output_dir, --name
p <- add_stage_args(p, "CNTFCT")     # the stage I/O contract
# your own method params — argparser directly (its add_argument requires `help`):
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

# select number of latent dims
pca_dim <- args$n_dim

# read cellids to subset on
cellids <- readLines(gzfile(args$filtered_cellids))
cat("length(cellids):", length(cellids), "\n")

# read H5AD into SCE
sce <- read_h5ad(args$rawdata_h5ad, as = "SingleCellExperiment")
sce <- sce[,cellids]

# read the filtered matrix
sce_filt <- TENxMatrix(args$normalized_selected_h5, group = "matrix")
rownames(sce_filt)

# For testing: subset cells and genes
# sce <- sce[, 1:1000]
# gene_idx <- sample(nrow(sce), 200)

# select the rownames from normalized_selected_h5
sce <- sce[rownames(sce_filt), ]
# filter genes with zero counts
keep <- rowSums(counts(sce)) > 0
sce  <- sce[keep, ]

cat("datasets loaded: running factorization\n")

set.seed(args$random_seed)
sce <- NewWave::newWave(Y = sce, K = pca_dim, n_gene_disp = 100, children = 4) # children = number of cores
# select latent dimensions / embeddings
res <- reducedDim(sce, "newWave")

# save embeddings
out_embeddings_tsv <- file.path(args$output_dir, sprintf("%s_reduced_dims.tsv", args$name))
 fwrite(data.frame(cell_id = rownames(res), res), out_embeddings_tsv,
    sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  wrote: %s\n", out_embeddings_tsv))
