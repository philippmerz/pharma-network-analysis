# ============================================================
# Hypothesis Testing (H2)
# H2: The QUALITY of university partners (CWUR) -> pharma patents
# Same testing approach as H1:
#   - Firm-level OLS on total patents
#   - Panel Poisson with FE (fixest::fepois)
# Only difference: regressor is tie quality (CWUR score / inv rank)
# ============================================================

library(dplyr)
library(readr)
library(igraph)
library(fixest)
library(modelsummary)
library(tibble)

cat("Hypothesis Testing (H2)\n")

# ----------------------------
# Load inputs produced by earlier steps
# ----------------------------
data  <- readRDS("Output/network_data.rds")
nodes <- data$nodes
config <- data$config

# graph might be stored as data$graph; build if needed
if (!is.null(data$graph)) {
  g <- data$graph
} else if (!is.null(data$edges)) {
  g <- graph_from_data_frame(data$edges, directed = FALSE,
                             vertices = nodes %>% transmute(name = participants))
} else {
  stop("network_data.rds must contain $graph or $edges")
}

panel <- read_csv("Output/panel.csv", show_col_types = FALSE)

# ----------------------------
# Helper: ensure we have the same columns H1 groups by
# (if panel.csv is minimal for some reason, join from nodes)
# ----------------------------
need_cols <- c("participants","has_uni_ties","uni_alliance_count","total_alliance_count",
               "nation","degree","betweenness","org_type")

missing_in_panel <- setdiff(need_cols, names(panel))

if (length(missing_in_panel) > 0) {
  cat("Note: panel.csv missing:", paste(missing_in_panel, collapse = ", "),
      "\n  -> joining these from nodes (if available)\n")
  
  # only select columns that actually exist in nodes
  avail_from_nodes <- intersect(missing_in_panel, names(nodes))
  
  if (length(avail_from_nodes) == 0) {
    stop("panel.csv is missing required columns and nodes does not contain them either.\n",
         "Check Output/panel.csv columns or rerun src/02_build_network.R.")
  }
  
  panel <- panel %>%
    left_join(nodes %>% select(participants, all_of(avail_from_nodes)),
              by = "participants")
}

# Pharma panel (same spirit as H1: start from panel + minimal join if needed)
pharma_panel <- panel %>%
  filter(org_type == "pharma")

# ----------------------------
# 1) Build tie-quality per pharma firm from university neighbors in graph
# ----------------------------
has_score <- "cwur_score" %in% names(nodes)
has_rank  <- "cwur_rank"  %in% names(nodes)
has_top   <- "is_top_university" %in% names(nodes)

if (!has_rank && !has_score) {
  stop("nodes does not contain cwur_rank or cwur_score. Recheck your preprocess CWUR merge.")
}

type_map <- setNames(nodes$org_type, nodes$participants)

get_uni_neighbors <- function(pharma_name) {
  neigh <- neighbors(g, pharma_name, mode = "all")
  neigh_names <- V(g)$name[neigh]
  neigh_names[type_map[neigh_names] == "university"]
}

pharma_names <- nodes %>%
  filter(org_type == "pharma") %>%
  pull(participants)

pharma_tie_quality <- bind_rows(lapply(pharma_names, function(p) {
  unis <- get_uni_neighbors(p)
  
  if (length(unis) == 0) {
    return(tibble(
      participants = p,
      n_uni_partners = 0L,
      tie_quality_score = NA_real_,
      tie_quality_rank_inv = NA_real_,
      share_top100 = NA_real_
    ))
  }
  
  uni_rows <- nodes %>%
    filter(participants %in% unis) %>%
    select(participants,
           any_of(c("cwur_score","cwur_rank","is_top_university")))
  
  # score: higher = better
  tq_score <- if (has_score) {
    if (all(is.na(uni_rows$cwur_score))) NA_real_ else mean(uni_rows$cwur_score, na.rm = TRUE)
  } else NA_real_
  
  # inverse rank: higher = better (since rank smaller is better)
  tq_rank_inv <- if (has_rank) {
    rinv <- ifelse(is.na(uni_rows$cwur_rank) | uni_rows$cwur_rank == 0, NA_real_, 1 / uni_rows$cwur_rank)
    if (all(is.na(rinv))) NA_real_ else mean(rinv, na.rm = TRUE)
  } else NA_real_
  
  share_top <- if (has_top) {
    if (all(is.na(uni_rows$is_top_university))) NA_real_
    else mean(as.numeric(uni_rows$is_top_university), na.rm = TRUE)
  } else NA_real_
  
  tibble(
    participants = p,
    n_uni_partners = dplyr::n_distinct(uni_rows$participants),
    tie_quality_score = tq_score,
    tie_quality_rank_inv = tq_rank_inv,
    share_top100 = share_top
  )
}))

# Choose main regressor like we discussed:
# prefer CWUR score if it exists and not all NA; else use inverse rank
use_score <- has_score && !all(is.na(pharma_tie_quality$tie_quality_score))
main_x <- if (use_score) "tie_quality_score" else "tie_quality_rank_inv"
cat("H2 main regressor:", main_x, "\n")

# Merge into pharma panel
pharma_panel_h2 <- pharma_panel %>%
  left_join(pharma_tie_quality, by = "participants")

# ----------------------------
# Descriptive stats (simple, like H1)
# ----------------------------
cat("\n-- Descriptive Statistics (H2) --\n\n")

tie_summary <- pharma_tie_quality %>%
  summarize(
    n_firms = n(),
    firms_with_uni_partners = sum(n_uni_partners > 0, na.rm = TRUE),
    mean_uni_partners = mean(n_uni_partners, na.rm = TRUE),
    mean_rank_inv = mean(tie_quality_rank_inv, na.rm = TRUE),
    mean_share_top100 = mean(share_top100, na.rm = TRUE),
    mean_score = mean(tie_quality_score, na.rm = TRUE)
  )

print(tie_summary)

# ----------------------------
# Cross-sectional regression (Firm Level) — EXACTLY like H1 structure
# ----------------------------
cat("\n-- Cross-Sectional Regression (Firm Level) --\n\n")

firm_data <- pharma_panel_h2 %>%
  group_by(participants, has_uni_ties, uni_alliance_count,
           total_alliance_count, nation, degree, betweenness) %>%
  summarize(
    total_patents = sum(patent_count),
    avg_patents_per_year = mean(patent_count),
    tie_quality_score = first(tie_quality_score),
    tie_quality_rank_inv = first(tie_quality_rank_inv),
    share_top100 = first(share_top100),
    n_uni_partners = first(n_uni_partners),
    .groups = "drop"
  )

# Model 1: Simple OLS - tie quality
m1 <- lm(as.formula(paste0("total_patents ~ ", main_x)), data = firm_data)

# Model 2: Add controls for network position
m2 <- lm(as.formula(paste0("total_patents ~ ", main_x, " + degree + betweenness")), data = firm_data)

# Model 3: Add country FE (same as H1)
m3 <- lm(as.formula(paste0("total_patents ~ ", main_x, " + degree + betweenness + nation")), data = firm_data)

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

# ----------------------------
# Panel regression (Firm x Year) — Poisson FE (same as H1)
# ----------------------------
cat("\n-- Panel Regression (Firm x Year) --\n\n")

p1 <- fepois(as.formula(paste0("patent_count ~ ", main_x, " | year")), data = pharma_panel_h2)
p2 <- fepois(as.formula(paste0("patent_count ~ ", main_x, " + degree + betweenness | year")), data = pharma_panel_h2)
p3 <- fepois(as.formula(paste0("patent_count ~ ", main_x, " + degree + betweenness | year + nation")), data = pharma_panel_h2)

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

# ----------------------------
# H2 Summary + Save (same style as H1)
# ----------------------------
cat("\n-- H2 Test Summary --\n\n")

coef_m2 <- summary(m2)$coefficients[main_x, ]
coef_p3 <- coeftable(p3)[main_x, ]

cat("H2: University tie quality -> Patent output (European pharma firms)\n\n")

cat("Cross-sectional (OLS with controls):\n")
cat(sprintf("  coef = %.3f, SE = %.3f, t = %.2f, p = %.4f\n",
            coef_m2[1], coef_m2[2], coef_m2[3], coef_m2[4]))

cat("\nPanel (Poisson with year + country FE):\n")
cat(sprintf("  coef = %.3f, SE = %.3f, z = %.2f, p = %.4f\n",
            coef_p3["Estimate"], coef_p3["Std. Error"],
            coef_p3["t value"], coef_p3["Pr(>|t|)"]))

h2_supported <- coef_m2[1] > 0 && coef_m2[4] < 0.05

cat("\n")
if (h2_supported) {
  cat("[OK] H2 SUPPORTED: Positive and significant relationship between\n")
  cat("  tie quality and patent output.\n")
} else if (coef_m2[1] > 0) {
  cat("[~] H2 PARTIALLY SUPPORTED: Positive but not significant at p < 0.05\n")
} else {
  cat("[X] H2 NOT SUPPORTED: No positive relationship found.\n")
}

h2_results <- list(
  tie_summary = tie_summary,
  pharma_tie_quality = pharma_tie_quality,
  firm_data = firm_data,
  models_cross_sectional = models_cs,
  models_panel = models_panel,
  h2_supported = h2_supported,
  main_regressor = main_x
)

saveRDS(h2_results, "Output/h2_results.rds")
cat("\n  Saved: Output/h2_results.rds\n")
