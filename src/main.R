# Pharma-University Alliance Network Analysis Pipeline

cat("Pharma-University Alliance Network Analysis\n\n")

source("src/01_preprocess.R")
source("src/02_build_network.R")
source("src/03_visualize.R")
source("src/hypothesis.R")

cat("\nPipeline complete. Outputs saved to: Output/\n")