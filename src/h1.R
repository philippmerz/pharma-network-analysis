library(tidyverse)
library(fixest)
library(modelsummary)

cat("Hypothesis Testing (H1)\n")

data <- readRDS("Output/network_data.rds")
nodes  <- data$nodes
config <- data$config

panel <- read_csv("Output/panel.csv", show_col_types = FALSE)

pharma_firms <- nodes %>%
  filter(org_type == "pharma") %>%
  select(participants, uni_alliance_count, total_alliance_count, 
         nation, degree, betweenness, eigen_centrality)

pharma_panel <- panel %>%
  filter(org_type == "pharma") %>%
  left_join(
    pharma_firms %>% select(participants, eigen_centrality),
by = "participants"
  )

cat("\n-- Descriptive Statistics --\n\n")

firm_summary <- pharma_firms %>%
  summarize(
    n_firms = n(),
    firms_with_uni_ties = sum(uni_alliance_count > 0),
    mean_uni_alliances = mean(uni_alliance_count),
    sd_uni_alliances = sd(uni_alliance_count),
    mean_total_alliances = mean(total_alliance_count),
    sd_total_alliances = sd(total_alliance_count)
  )

print(firm_summary)

patent_summary <- pharma_panel %>%
  group_by(has_uni_ties) %>%
  summarize(
    n_obs = n(),
    mean_patents = mean(patent_count),
    sd_patents = sd(patent_count),
    total_patents = sum(patent_count),
    .groups = "drop"
  )

cat("\nPatent output by university tie status:\n")
print(patent_summary)

cat("\n-- Cross-Sectional Regression (Firm Level) --\n\n")

firm_data <- pharma_panel %>%
  group_by(participants, has_uni_ties, uni_alliance_count, 
           total_alliance_count, nation, degree, betweenness) %>%
  summarize(
    total_patents = sum(patent_count),
    avg_patents_per_year = mean(patent_count),
    .groups = "drop"
  )

# Model 1: Simple OLS - uni alliance count
m1 <- lm(total_patents ~ uni_alliance_count, data = firm_data)

# Model 2: Add controls for network position
m2 <- lm(total_patents ~ uni_alliance_count + degree + betweenness, data = firm_data)

# Model 3: Add country fixed effects
m3 <- lm(total_patents ~ uni_alliance_count + degree + betweenness + nation, data = firm_data)

models_cs <- list(
  "OLS (1)" = m1,
  "OLS (2)" = m2,
  "OLS + Country FE" = m3
)

modelsummary(
  models_cs,
  stars = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
  gof_map = c("nobs", "r.squared", "adj.r.squared"),
  coef_omit = "nation",
  add_rows = tibble(
    term = "Country FE",
    `OLS (1)` = "No",
    `OLS (2)` = "No",
    `OLS + Country FE` = "Yes"
  ),
output = "markdown"
) %>% print()

cat("\n-- Panel Regression (Firm x Year) --\n\n")

# Poisson panel with fixed effects
p1 <- fepois(patent_count ~ uni_alliance_count | year, data = pharma_panel)
p2 <- fepois(patent_count ~ uni_alliance_count + degree + betweenness | year, data = pharma_panel)
p3 <- fepois(patent_count ~ uni_alliance_count + degree + betweenness | year + nation, data = pharma_panel)

models_panel <- list(
  "Poisson (1)" = p1,
  "Poisson (2)" = p2,
  "Poisson + FE" = p3
)

modelsummary(
  models_panel,
  stars = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
  gof_map = c("nobs", "r.squared.within"),
output = "markdown"
) %>% print()

cat("\n-- H1 Test Summary --\n\n")

# Extract coefficients (fixest uses coeftable(), base R uses summary()$coefficients)
coef_m2 <- summary(m2)$coefficients["uni_alliance_count", ]
coef_p3 <- coeftable(p3)["uni_alliance_count", ]

cat("H1: University alliances -> Patent output (European pharma firms)\n\n")

cat("Cross-sectional (OLS with controls):\n")
cat(sprintf("  coef = %.3f, SE = %.3f, t = %.2f, p = %.4f\n",
            coef_m2[1], coef_m2[2], coef_m2[3], coef_m2[4]))

cat("\nPanel (Poisson with year + country FE):\n")
cat(sprintf("  coef = %.3f, SE = %.3f, z = %.2f, p = %.4f\n",
            coef_p3["Estimate"], coef_p3["Std. Error"], 
            coef_p3["t value"], coef_p3["Pr(>|t|)"]))

h1_supported <- coef_m2[1] > 0 && coef_m2[4] < 0.05

cat("\n")
if (h1_supported) {
  cat("[OK] H1 SUPPORTED: Positive and significant relationship between\n")
  cat("  university alliance count and patent output.\n")
} else if (coef_m2[1] > 0) {
  cat("[~] H1 PARTIALLY SUPPORTED: Positive but not significant at p < 0.05\n")
} else {
  cat("[X] H1 NOT SUPPORTED: No positive relationship found.\n")
}

h1_results <- list(
  firm_summary = firm_summary,
  patent_summary = patent_summary,
  models_cross_sectional = models_cs,
  models_panel = models_panel,
  h1_supported = h1_supported
)

saveRDS(h1_results, "Output/h1_results.rds")

cat("\n  Saved: Output/h1_results.rds\n")
