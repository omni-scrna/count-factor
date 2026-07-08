#!/usr/bin/env Rscript

library(argparse)

# Source main functions
source("src/main.R")

# Parse command line arguments
parser <- ArgumentParser(description="OmniBenchmark module")

# Required by OmniBenchmark
parser$add_argument("--output_dir", dest="output_dir", type="character", required=TRUE,
                   help="Output directory for results")
parser$add_argument("--name", dest="name", type="character", required=TRUE,
                   help="Module name/identifier")
# Add your custom input arguments here
# Example:
# parser$add_argument("--input", dest="input", type="character", help="Input file")

args <- parser$parse_args()

cat("Output directory:", args$output_dir, "\n")
cat("Module name:", args$name, "\n")

# TODO: Implement your module logic
# Process the data using main function
process_data(args)
