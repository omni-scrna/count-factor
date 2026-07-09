#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Matrix)
  library(HDF5Array)
  library(BiocSingular)
  library(NewWave)
  library(anndataR)
  library(SingleCellExperiment)
  library(data.table)
})
cat("OK\n")