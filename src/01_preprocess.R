# 01_preprocess.R - Load and prepare alliance, patent, and ranking data

library(tidyverse)
library(jsonlite)
library(fuzzyjoin)
library(stringdist)

# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────

config <- list(
  # Paths
  sdc_path        = "data/SDC_data_2021.rds",
  sic_path        = "data/SIC_codes_mapping.json",
  uspto_dir       = "data/USPTO",
  cwur_path       = "data/cwurData.csv",
  output_dir      = "Output",
  
  # Time periods
  alliance_years  = 2015:2021,
  patent_years    = 2015:2024,
  uspto_years     = 2000:2024,
  
  # Thresholds
  top_uni_rank    = 100,
  fuzzy_max_dist  = 0.10,
  
  # Sampling
  match_control_by_deal_count = TRUE,
  
  # Countries
  pharma_countries = c(
    "United Kingdom", "England", "Scotland", "Wales", "Northern Ireland",
    "Germany", "France"
  ),
  
  # SIC codes for pharma

  pharma_sics = c("2833", "2834", "2835", "2836")
)

if (!dir.exists(config$output_dir)) dir.create(config$output_dir)

cat("Step 1: Data Preprocessing\n")

# ─────────────────────────────────────────────────────────────────────────────
# Helper functions
# ─────────────────────────────────────────────────────────────────────────────

clean_org_name <- function(x) {
  x %>%
    str_to_upper() %>%
    str_replace_all("[^A-Z0-9 ]", " ") %>%
    str_replace_all("\\b(BV|NV|INC|LTD|AG|SA|GMBH|SAS|SARL|SPA|PLC|LLC|CORP|CO|COMPANY)\\b", " ") %>%
    str_squish()
}

# ─────────────────────────────────────────────────────────────────────────────
# Load raw data
# ─────────────────────────────────────────────────────────────────────────────

sdc_raw <- readRDS(config$sdc_path)
sic_mapping <- fromJSON(config$sic_path)
sic_lookup <- setNames(sic_mapping$Industry_Title, as.character(sic_mapping$SIC_Code))

# University detection patterns
uni_pattern <- regex(paste(c(
  "\\buniversity\\b", "\\buniv\\b", "\\buniversit", "\\bcollege\\b",
  "\\becole\\b", "\\bécole\\b", "\\bpolitecnico\\b", "\\bhochschule\\b",
  "\\buniversiteit\\b", "\\buniversidad\\b", "\\buniversità\\b",
  "\\bcnrs\\b", "\\binserm\\b", "\\bfraunhofer\\b", "\\binria\\b",
  "\\bmax.?planck\\b", "\\beth.?zurich\\b", "\\beth.?zürich\\b",
  "\\bkarolinska\\b", "\\bweizmann\\b", "\\bmit\\b"
), collapse = "|"), ignore_case = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# Classify organizations and filter alliances
# ─────────────────────────────────────────────────────────────────────────────

sdc <- sdc_raw %>%
  mutate(
    year = as.integer(substr(date_announced, 1, 4)),
    sic_label = sic_lookup[as.character(SIC_primary)],
    is_pharma = SIC_primary %in% config$pharma_sics,
    is_university = str_detect(tolower(participants), uni_pattern)
  )

sdc_filtered <- sdc %>%
  filter(
    !is.na(year),
    year %in% config$alliance_years,
    status %in% c("Completed/Signed", "Pending", "Extended", "Renegotiated")
  )

# ─────────────────────────────────────────────────────────────────────────────
# Identify sample groups
# ─────────────────────────────────────────────────────────────────────────────

pharma_target <- sdc_filtered %>%
  filter(is_pharma, participant_nation %in% config$pharma_countries)

pharma_deal_ids <- unique(pharma_target$deal_number)

uni_in_pharma_deals <- sdc_filtered %>%
  filter(deal_number %in% pharma_deal_ids, is_university)

mixed_deal_ids <- intersect(pharma_deal_ids, unique(uni_in_pharma_deals$deal_number))

# Group A: Pharma firms with university ties
pharma_with_uni <- pharma_target %>%
  filter(deal_number %in% mixed_deal_ids) %>%
  group_by(participants) %>%
  summarize(
    n_deals = n_distinct(deal_number),
    n_uni_deals = n_distinct(deal_number),
    nation = first(participant_nation),
    .groups = "drop"
  ) %>%
  mutate(has_uni_ties = TRUE)

# Group B: Universities partnered with selected pharma
selected_deal_ids <- sdc_filtered %>%
  filter(participants %in% pharma_with_uni$participants) %>%
  pull(deal_number) %>%
  unique()

universities <- sdc_filtered %>%
  filter(deal_number %in% intersect(selected_deal_ids, mixed_deal_ids), is_university) %>%
  group_by(participants) %>%
  summarize(
    n_deals = n_distinct(deal_number),
    n_uni_deals = n_distinct(deal_number),
    nation = first(participant_nation),
    .groups = "drop"
  ) %>%
  mutate(has_uni_ties = NA)

# Group C: Pharma firms without university ties (control group)
pharma_no_uni_pool <- pharma_target %>%
  filter(!participants %in% pharma_with_uni$participants) %>%
  group_by(participants) %>%
  summarize(
    n_deals = n_distinct(deal_number),
    n_uni_deals = 0L,
    nation = first(participant_nation),
    .groups = "drop"
  )

# Sample control group matched by deal count distribution
if (config$match_control_by_deal_count) {
  deal_bins <- unique(c(0, quantile(pharma_with_uni$n_deals, probs = seq(0, 1, 0.25)), Inf)) %>% sort()
  bin_labels <- paste0("Q", seq_along(deal_bins[-1]))
  
  pharma_with_uni_binned <- pharma_with_uni %>%
    mutate(deal_bin = cut(n_deals, breaks = deal_bins, labels = bin_labels, include.lowest = TRUE))
  
  bin_targets <- pharma_with_uni_binned %>%
    count(deal_bin, name = "target_n") %>%
    filter(!is.na(deal_bin))
  
  set.seed(42)
  pharma_no_uni <- pharma_no_uni_pool %>%
    mutate(deal_bin = cut(n_deals, breaks = deal_bins, labels = bin_labels, include.lowest = TRUE)) %>%
    group_by(deal_bin) %>%
    nest() %>%
    left_join(bin_targets, by = "deal_bin") %>%
    mutate(
      target_n = replace_na(target_n, 0),
      sampled = map2(data, target_n, ~ {
        if (is.null(.x) || nrow(.x) == 0 || .y == 0) tibble() 
        else slice_sample(.x, n = min(.y, nrow(.x)))
      })
    ) %>%
    select(sampled) %>%
    unnest(sampled) %>%
    ungroup() %>%
    select(-deal_bin)
  
  pharma_with_uni <- pharma_with_uni_binned %>% select(-deal_bin)
} else {
  set.seed(42)
  pharma_no_uni <- pharma_no_uni_pool %>%
    slice_sample(n = min(nrow(pharma_with_uni), nrow(pharma_no_uni_pool)))
}

pharma_no_uni <- pharma_no_uni %>% mutate(has_uni_ties = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# Combine all nodes
# ─────────────────────────────────────────────────────────────────────────────

nodes <- bind_rows(
  pharma_with_uni %>% mutate(org_type = "pharma"),
  pharma_no_uni %>% mutate(org_type = "pharma"),
  universities %>% mutate(org_type = "university")
) %>%
  distinct(participants, .keep_all = TRUE) %>%
  mutate(
    display_type = case_when(
      org_type == "university"      ~ "University/Research",
      has_uni_ties == TRUE          ~ "Pharma WITH Uni Ties",
      has_uni_ties == FALSE         ~ "Pharma WITHOUT Uni Ties",
      TRUE                          ~ "Other"
    ),
    name_clean = clean_org_name(participants)
  )

# ─────────────────────────────────────────────────────────────────────────────
# Match universities to CWUR rankings (fuzzy)
# ─────────────────────────────────────────────────────────────────────────────

cwur_raw <- read_csv(config$cwur_path, show_col_types = FALSE)

cwur <- cwur_raw %>%
  group_by(institution) %>%
  filter(year == max(year)) %>%
  ungroup() %>%
  select(institution, country, world_rank, score, year) %>%
  mutate(name_clean = clean_org_name(institution))

unis_to_match <- nodes %>%
  filter(org_type == "university") %>%
  select(participants, name_clean) %>%
  distinct()

cwur_matches <- stringdist_left_join(
  unis_to_match, cwur,
  by = "name_clean", method = "jw", max_dist = config$fuzzy_max_dist
) %>%
  group_by(participants) %>%
  slice_min(world_rank, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(
    participants, 
    cwur_institution = institution, 
    cwur_country = country,
    cwur_rank = world_rank, 
    cwur_score = score
  )

nodes <- nodes %>%
  left_join(cwur_matches, by = "participants") %>%
  mutate(is_top_university = !is.na(cwur_rank) & cwur_rank <= config$top_uni_rank)

# ─────────────────────────────────────────────────────────────────────────────
# Load and process USPTO patent data
# ─────────────────────────────────────────────────────────────────────────────

read_uspto_year <- function(year) {
  path <- file.path(config$uspto_dir, paste0(year, ".csv"))
  if (!file.exists(path)) return(tibble())
  read_csv(path, show_col_types = FALSE) %>%
    select(patent_number, grant_year, assignee) %>%
    mutate(file_year = year)
}

uspto <- map_dfr(config$uspto_years, read_uspto_year) %>%
  filter(
    grant_year %in% config$patent_years,
    !is.na(assignee), 
    assignee != ""
  ) %>%
  mutate(
    assignee = str_trim(assignee), 
    assignee_clean = clean_org_name(assignee)
  )

patents_by_assignee_year <- uspto %>%
  group_by(assignee_clean, grant_year) %>%
  summarize(patent_count = n_distinct(patent_number), .groups = "drop")

# ─────────────────────────────────────────────────────────────────────────────
# Match organizations to patents (fuzzy)
# ─────────────────────────────────────────────────────────────────────────────

orgs_to_match <- nodes %>% 
  select(participants, name_clean) %>% 
  distinct()

patent_matches <- stringdist_left_join(
  orgs_to_match, patents_by_assignee_year,
  by = c("name_clean" = "assignee_clean"), 
  method = "jw", 
  max_dist = config$fuzzy_max_dist
)

patents_by_org_year <- patent_matches %>%
  filter(!is.na(grant_year)) %>%
  group_by(participants, grant_year) %>%
  summarize(patent_count = sum(patent_count, na.rm = TRUE), .groups = "drop")

patents_by_org_total <- patents_by_org_year %>%
  group_by(participants) %>%
  summarize(
    patent_count_total = sum(patent_count), 
    patent_years = n(), 
    .groups = "drop"
  )

patents_by_org_wide <- patents_by_org_year %>%
  pivot_wider(
    id_cols = participants, 
    names_from = grant_year,
    values_from = patent_count, 
    names_prefix = "patents_", 
    values_fill = 0
  )

# ─────────────────────────────────────────────────────────────────────────────
# Save intermediate data
# ─────────────────────────────────────────────────────────────────────────────

intermediate_data <- list(
  nodes              = nodes,
  sdc_filtered       = sdc_filtered,
  patents_by_org_year   = patents_by_org_year,
  patents_by_org_total  = patents_by_org_total,
  patents_by_org_wide   = patents_by_org_wide,
  cwur_matches       = cwur_matches,
  config             = config
)

saveRDS(intermediate_data, file.path(config$output_dir, "intermediate_data.rds"))
write_csv(cwur_matches, file.path(config$output_dir, "cwur_matches.csv"))

cat("  Saved:", file.path(config$output_dir, "intermediate_data.rds"), "\n")