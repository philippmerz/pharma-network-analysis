
# Hypothesis Testing (H2) — QS version
# H2: QUALITY of university partners (QS) -> pharma patents
# Same testing approach as H1:
#   - Firm-level OLS on total patents
#   - Panel Poisson with FE (fixest::fepois)
#
# QS adaptation
# Includes FIXES:
#   1) Robustly joins missing controls from nodes into panel if needed
#   2) Poisson summary print uses z value + Pr(>|z|)
#   3) Robustness checks:
#        A) Linked-only sample (n_uni_partners > 0)
#        B) Rank-based regressor (tie_quality_rank_inv)

library(dplyr)
library(readr)
library(stringr)
library(stringdist)
library(fuzzyjoin)
library(igraph)
library(fixest)
library(modelsummary)
library(tibble)

cat("Hypothesis Testing (H2) — QS rankings\n")

# ----------------------------
# Helpers
# ----------------------------
clean_org_name <- function(x) {
  x <- toupper(as.character(x))
  x <- gsub("[^A-Z0-9 ]", " ", x)
  x <- gsub("\\b(BV|NV|INC|LTD|AG|SA|GMBH|SAS|SARL|SPA|PLC|LLC|CORP|CO|COMPANY)\\b", " ", x, perl = TRUE)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

# safely read a column that might have spaces/case differences
get_col <- function(df, candidates) {
  nms <- names(df)
  idx <- match(tolower(candidates), tolower(nms))
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) return(NULL)
  nms[idx[1]]
}

# ----------------------------
# 0) Load inputs produced by earlier steps
# ----------------------------
net   <- readRDS("Output/network_data.rds")
nodes <- net$nodes
config <- net$config

# graph might be stored as net$graph; build if needed
if (!is.null(net$graph)) {
  g <- net$graph
} else if (!is.null(net$edges)) {
  g <- graph_from_data_frame(
    net$edges, directed = FALSE,
    vertices = nodes %>% transmute(name = participants)
  )
} else {
  stop("Output/network_data.rds must contain $graph or $edges")
}

panel <- read_csv("Output/panel.csv", show_col_types = FALSE)

# ----------------------------
# Ensure QS columns exist in nodes; if not, build themql from data/qs_rank.csv
# ----------------------------
need_qs_cols <- c("qs_rank", "qs_score", "is_top_university")
has_qs_already <- all(need_qs_cols %in% names(nodes))

if (!has_qs_already) {
  cat("QS columns not found in nodes -> building them from data/qs_rank.csv via fuzzy match...\n")
  
  qs_path <- "data/qs_rank.csv"
  if (!file.exists(qs_path)) stop("Missing file: ", qs_path)
  
  qs_raw <- read_csv(qs_path, show_col_types = FALSE)
  
  # Identify required columns (robust to spaces/case)
  col_rank <- get_col(qs_raw, c("Rank", "rank"))
  col_inst <- get_col(qs_raw, c("institution", "Institution", "university", "University"))
  col_score <- get_col(qs_raw, c(
    "score scaled", "Score Scaled", "score_scaled",   # preferred overall score
    "cpf score", "cpf_score", "cpfscore"              # fallback if needed
  ))
  col_country <- get_col(qs_raw, c("location", "Location", "country", "Country"))
  
  if (is.null(col_rank) || is.null(col_inst)) {
    stop("qs_rank.csv must contain columns like Rank and institution. Found: ",
         paste(names(qs_raw), collapse = ", "))
  }
  if (is.null(col_score)) {
    stop("qs_rank.csv: could not find a usable score column. Expected one of: ",
         "'score scaled' or 'cpf score'. Found: ", paste(names(qs_raw), collapse = ", "))
  }
  
  qs <- qs_raw %>%
    transmute(
      qs_institution = .data[[col_inst]],
      qs_rank = suppressWarnings(as.numeric(.data[[col_rank]])),
      qs_score = suppressWarnings(as.numeric(.data[[col_score]])),
      qs_country = if (!is.null(col_country)) as.character(.data[[col_country]]) else NA_character_,
      name_clean = clean_org_name(.data[[col_inst]])
    )
  
  # ensure nodes has name_clean (needed for fuzzy match)
  if (!("name_clean" %in% names(nodes))) {
    nodes <- nodes %>% mutate(name_clean = clean_org_name(participants))
  }
  
  unis_to_match <- nodes %>%
    filter(org_type == "university") %>%
    select(participants, name_clean) %>%
    distinct()
  
  max_dist <- if (!is.null(config$fuzzy_max_dist)) config$fuzzy_max_dist else 0.10
  
  qs_matches <- stringdist_left_join(
    unis_to_match, qs,
    by = "name_clean", method = "jw", max_dist = max_dist
  ) %>%
    group_by(participants) %>%
    slice_min(qs_rank, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(participants, qs_rank, qs_score, qs_country)
  
  nodes <- nodes %>%
    left_join(qs_matches, by = "participants") %>%
    mutate(
      is_top_university = if (!is.null(config$top_uni_rank)) {
        !is.na(qs_rank) & qs_rank <= config$top_uni_rank
      } else {
        !is.na(qs_rank) & qs_rank <= 100
      }
    )
  
  if (!dir.exists("Output")) dir.create("Output")
  write_csv(qs_matches, "Output/qs_matches.csv")
  cat("Saved: Output/qs_matches.csv\n")
}

# ----------------------------
# Ensure we have the same columns H1 groups by
# (if panel.csv is minimal, join from nodes)
# ----------------------------
need_cols <- c("participants","has_uni_ties","uni_alliance_count","total_alliance_count",
               "nation","degree","betweenness","org_type","year","patent_count")

missing_in_panel <- setdiff(need_cols, names(panel))
if (length(missing_in_panel) > 0) {
  cat("Note: panel.csv missing:", paste(missing_in_panel, collapse = ", "),
      "\n  -> joining these from nodes (if available)\n")
  
  joinable <- setdiff(missing_in_panel, c("year","patent_count"))
  avail_from_nodes <- intersect(joinable, names(nodes))
  
  if (length(avail_from_nodes) > 0) {
    panel <- panel %>%
      left_join(nodes %>% select(participants, all_of(avail_from_nodes)),
                by = "participants")
  }
}

# Pharma panel
pharma_panel <- panel %>% filter(org_type == "pharma")

# ----------------------------
# 1) Build tie-quality per pharma firm from university neighbors in graph (QS)
# ----------------------------
has_score <- "qs_score" %in% names(nodes)
has_rank  <- "qs_rank"  %in% names(nodes)
has_top   <- "is_top_university" %in% names(nodes)

if (!has_rank && !has_score) {
  stop("nodes does not contain qs_rank or qs_score after QS merge. Check qs_rank.csv + matching.")
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
    select(participants, any_of(c("qs_score","qs_rank","is_top_university")))
  
  tq_score <- if (has_score) {
    if (all(is.na(uni_rows$qs_score))) NA_real_ else mean(uni_rows$qs_score, na.rm = TRUE)
  } else NA_real_
  
  tq_rank_inv <- if (has_rank) {
    rinv <- ifelse(is.na(uni_rows$qs_rank) | uni_rows$qs_rank == 0, NA_real_, 1 / uni_rows$qs_rank)
    if (all(is.na(rinv))) NA_real_ else mean(rinv, na.rm = TRUE)
  } else NA_real_
  
  share_top <- if (has_top) {
    if (all(is.na(uni_rows$is_top_university))) NA_real_ else mean(as.numeric(uni_rows$is_top_university), na.rm = TRUE)
  } else NA_real_
  
  tibble(
    participants = p,
    n_uni_partners = dplyr::n_distinct(uni_rows$participants),
    tie_quality_score = tq_score,
    tie_quality_rank_inv = tq_rank_inv,
    share_top100 = share_top
  )
}))

use_score <- has_score && !all(is.na(pharma_tie_quality$tie_quality_score))
main_x <- if (use_score) "tie_quality_score" else "tie_quality_rank_inv"
cat("H2 main regressor:", main_x, "\n")

# Merge tie quality into pharma panel
pharma_panel_h2 <- pharma_panel %>%
  left_join(pharma_tie_quality, by = "participants")

# ----------------------------
# Descriptive stats
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
# Cross-sectional regression (Firm Level) — like H1
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

m1 <- lm(as.formula(paste0("total_patents ~ ", main_x)), data = firm_data)
m2 <- lm(as.formula(paste0("total_patents ~ ", main_x, " + degree + betweenness")), data = firm_data)
m3 <- lm(as.formula(paste0("total_patents ~ ", main_x, " + degree + betweenness + nation")), data = firm_data)

models_cs <- list("OLS (1)" = m1, "OLS (2)" = m2, "OLS + Country FE" = m3)

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
# Panel regression (Firm x Year) — Poisson FE (like H1)
# ----------------------------
cat("\n-- Panel Regression (Firm x Year) --\n\n")

p1 <- fepois(as.formula(paste0("patent_count ~ ", main_x, " | year")), data = pharma_panel_h2)
p2 <- fepois(as.formula(paste0("patent_count ~ ", main_x, " + degree + betweenness | year")), data = pharma_panel_h2)
p3 <- fepois(as.formula(paste0("patent_count ~ ", main_x, " + degree + betweenness | year + nation")), data = pharma_panel_h2)

models_panel <- list("Poisson (1)" = p1, "Poisson (2)" = p2, "Poisson + FE" = p3)

modelsummary(
  models_panel,
  stars = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
  gof_map = c("nobs", "r.squared.within"),
  output = "markdown"
) %>% print()

# ----------------------------
# Robustness A: Linked-only sample (firms with ≥1 university partner)
# ----------------------------
cat("\n-- Robustness A: Linked-only sample (n_uni_partners > 0) --\n\n")

pharma_panel_h2_linked <- pharma_panel_h2 %>% filter(n_uni_partners > 0)

p1L <- fepois(as.formula(paste0("patent_count ~ ", main_x, " | year")), data = pharma_panel_h2_linked)
p2L <- fepois(as.formula(paste0("patent_count ~ ", main_x, " + degree + betweenness | year")), data = pharma_panel_h2_linked)
p3L <- fepois(as.formula(paste0("patent_count ~ ", main_x, " + degree + betweenness | year + nation")), data = pharma_panel_h2_linked)

modelsummary(
  list("Poisson linked (1)" = p1L, "Poisson linked (2)" = p2L, "Poisson linked + FE" = p3L),
  stars = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
  gof_map = c("nobs", "r.squared.within"),
  output = "markdown"
) %>% print()

# ----------------------------
# Robustness B: Rank-based regressor (inverse rank)
# ----------------------------
cat("\n-- Robustness B: Rank-based regressor (tie_quality_rank_inv) --\n\n")

p1R <- fepois(patent_count ~ tie_quality_rank_inv | year, data = pharma_panel_h2)
p2R <- fepois(patent_count ~ tie_quality_rank_inv + degree + betweenness | year, data = pharma_panel_h2)
p3R <- fepois(patent_count ~ tie_quality_rank_inv + degree + betweenness | year + nation, data = pharma_panel_h2)

modelsummary(
  list("Poisson rank (1)" = p1R, "Poisson rank (2)" = p2R, "Poisson rank + FE" = p3R),
  stars = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
  gof_map = c("nobs", "r.squared.within"),
  output = "markdown"
) %>% print()

# ----------------------------
# H2 Summary + Save
# ----------------------------
cat("\n-- H2 Test Summary --\n\n")

coef_m2 <- summary(m2)$coefficients[main_x, ]
coef_p3 <- coeftable(p3)[main_x, ]

cat("H2: University tie quality (QS) -> Patent output (European pharma firms)\n\n")

cat("Cross-sectional (OLS with controls):\n")
cat(sprintf("  coef = %.3f, SE = %.3f, t = %.2f, p = %.4f\n",
            coef_m2[1], coef_m2[2], coef_m2[3], coef_m2[4]))

cat("\nPanel (Poisson with year + country FE):\n")
cat(sprintf("  coef = %.3f, SE = %.3f, z = %.2f, p = %.4f\n",
            coef_p3["Estimate"], coef_p3["Std. Error"],
            coef_p3["z value"],  coef_p3["Pr(>|z|)"]))

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
  robustness = list(
    linked_only = list(p1L = p1L, p2L = p2L, p3L = p3L),
    rank_based  = list(p1R = p1R, p2R = p2R, p3R = p3R)
  ),
  h2_supported = h2_supported,
  main_regressor = main_x
)

saveRDS(h2_results, "Output/h2_results_qs.rds")
cat("\nSaved: Output/h2_results_qs.rds\n")
