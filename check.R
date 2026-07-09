#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(rhdf5)
  library(DelayedArray)
  library(Matrix)
  library(optparse)
  library(yaml)
  library(data.table)
  library(NewWave)
})
cat("OK\n")