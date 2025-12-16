# Build alliance network and compute centrality metrics

library(tidyverse)
library(igraph)

cat("Step 2: Network Construction\n")

data <- readRDS("Output/intermediate_data.rds")

nodes             <- data$nodes
sdc_filtered      <- data$sdc_filtered
patents_by_org_year  <- data$patents_by_org_year
patents_by_org_total <- data$patents_by_org_total
patents_by_org_wide  <- data$patents_by_org_wide
config            <- data$config

relevant_deals <- sdc_filtered %>%
  filter(participants %in% nodes$participants) %>%
  pull(deal_number) %>%
  unique()

edges <- sdc_filtered %>%
  filter(participants %in% nodes$participants, deal_number %in% relevant_deals) %>%
  select(deal_number, participants) %>%
  distinct() %>%
  group_by(deal_number) %>%
  filter(n() >= 2) %>%
  summarize(orgs = list(unique(participants)), .groups = "drop") %>%
  mutate(pairs = map(orgs, ~ {
    if (length(.x) < 2) return(tibble())
    combn(.x, 2) %>% 
      t() %>% 
      as_tibble(.name_repair = ~c("from", "to"))
  })) %>%
  select(deal_number, pairs) %>%
  unnest(pairs)

edges_weighted <- edges %>% 
  count(from, to, name = "weight")

g <- graph_from_data_frame(edges_weighted, directed = FALSE, vertices = nodes$participants)
g <- simplify(g, remove.multiple = TRUE)

V(g)$degree           <- degree(g)
V(g)$betweenness      <- betweenness(g, normalized = TRUE)
V(g)$closeness        <- closeness(g, normalized = TRUE)
V(g)$eigen_centrality <- eigen_centrality(g)$vector

nodes_enriched <- nodes %>%
  left_join(patents_by_org_total, by = "participants") %>%
  left_join(patents_by_org_wide, by = "participants") %>%
  mutate(
    degree           = V(g)$degree[match(participants, V(g)$name)],
    betweenness      = V(g)$betweenness[match(participants, V(g)$name)],
    closeness        = V(g)$closeness[match(participants, V(g)$name)],
    eigen_centrality = V(g)$eigen_centrality[match(participants, V(g)$name)],
    patent_count_total = replace_na(patent_count_total, 0),
    patent_years       = replace_na(patent_years, 0)
  ) %>%
  rename(
    uni_alliance_count   = n_uni_deals, 
    total_alliance_count = n_deals
  )

panel <- expand_grid(
  participants = nodes$participants,
  year = config$patent_years
) %>%
  left_join(
    patents_by_org_year %>% rename(year = grant_year), 
    by = c("participants", "year")
  ) %>%
  left_join(
    nodes_enriched %>% 
      select(participants, org_type, display_type, has_uni_ties,
             uni_alliance_count, total_alliance_count, nation,
             degree, betweenness, is_top_university, qs_rank),
    by = "participants"
  ) %>%
  mutate(patent_count = replace_na(patent_count, 0))

summary_stats <- tibble(
  metric = c(
    "nodes", "edges", "density", "components", "avg_degree", "transitivity",
    "pharma_with_uni", "pharma_without_uni", "universities",
    "top_universities", "total_patents", "orgs_with_patents"
  ),
  value = c(
    vcount(g), 
    ecount(g), 
    round(edge_density(g), 4), 
    count_components(g),
    round(mean(degree(g)), 2), 
    round(transitivity(g, type = "global"), 4),
    sum(nodes_enriched$display_type == "Pharma WITH Uni Ties"),
    sum(nodes_enriched$display_type == "Pharma WITHOUT Uni Ties"),
    sum(nodes_enriched$display_type == "University/Research"),
    sum(nodes_enriched$is_top_university, na.rm = TRUE),
    sum(nodes_enriched$patent_count_total),
    sum(nodes_enriched$patent_count_total > 0)
  )
)

output_dir <- config$output_dir

write_csv(nodes_enriched, file.path(output_dir, "nodes.csv"))
write_csv(edges_weighted, file.path(output_dir, "edges.csv"))
write_csv(panel, file.path(output_dir, "panel.csv"))
write_csv(summary_stats, file.path(output_dir, "stats.csv"))

saveRDS(
  list(graph = g, nodes = nodes_enriched, config = config),
  file.path(output_dir, "network_data.rds")
)

cat("  Saved network data to:", output_dir, "\n")