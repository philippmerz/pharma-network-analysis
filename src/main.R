# main.R - Pharma-University Alliance Network Analysis Pipeline
#
# Builds a network of pharmaceutical firms and universities based on 
# alliance data, enriched with patent counts and university rankings.

cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("  Pharma-University Alliance Network Analysis\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

source("src/01_preprocess.R")
source("src/02_build_network.R")
source("src/03_visualize.R")

cat("\n═══════════════════════════════════════════════════════════════════════════\n")
cat("  Pipeline complete. Outputs in: Output/\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")