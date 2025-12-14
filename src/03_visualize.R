# 03_visualize.R - Network visualization

library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)

# ─────────────────────────────────────────────────────────────────────────────
# Load network data
# ─────────────────────────────────────────────────────────────────────────────

cat("Step 3: Visualization\n")

data <- readRDS("Output/network_data.rds")

g      <- data$graph
nodes  <- data$nodes
config <- data$config

# ─────────────────────────────────────────────────────────────────────────────
# Configure node attributes
# ─────────────────────────────────────────────────────────────────────────────

V(g)$display_type  <- nodes$display_type[match(V(g)$name, nodes$participants)]
V(g)$patent_count  <- nodes$patent_count_total[match(V(g)$name, nodes$participants)]
V(g)$is_top_uni    <- nodes$is_top_university[match(V(g)$name, nodes$participants)]

V(g)$color_type <- case_when(
  V(g)$is_top_uni == TRUE                        ~ "Top University (CWUR)",
  V(g)$display_type == "University/Research"     ~ "University/Research",
  V(g)$display_type == "Pharma WITH Uni Ties"    ~ "Pharma WITH Uni Ties",
  V(g)$display_type == "Pharma WITHOUT Uni Ties" ~ "Pharma WITHOUT Uni Ties",
  TRUE                                           ~ "Other"
)

V(g)$node_size <- sapply(seq_along(V(g)), function(i) {
  dt <- V(g)$display_type[i]
  patents <- V(g)$patent_count[i]
  is_top <- V(g)$is_top_uni[i]
  
  if (dt %in% c("Pharma WITH Uni Ties", "Pharma WITHOUT Uni Ties")) {
    3 + 1.5 * sqrt(patents)
  } else if (isTRUE(is_top)) {
    7
  } else {
    4
  }
})

# ─────────────────────────────────────────────────────────────────────────────
# Style definitions
# ─────────────────────────────────────────────────────────────────────────────

node_colors <- c(
  "Pharma WITH Uni Ties"    = "#27AE60",
  "Pharma WITHOUT Uni Ties" = "#E74C3C",
  "University/Research"     = "#3498DB",
  "Top University (CWUR)"   = "#FFD700"
)

node_shapes <- c(
  "Pharma WITH Uni Ties"    = 16,
  "Pharma WITHOUT Uni Ties" = 16,
  "University/Research"     = 17,
  "Top University (CWUR)"   = 17
)

n_pharma_with    <- sum(nodes$display_type == "Pharma WITH Uni Ties")
n_pharma_without <- sum(nodes$display_type == "Pharma WITHOUT Uni Ties")
n_unis           <- sum(nodes$display_type == "University/Research")
n_top_unis       <- sum(nodes$is_top_university, na.rm = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# Build plot
# ─────────────────────────────────────────────────────────────────────────────

set.seed(42)
g_tidy <- as_tbl_graph(g)

alliance_range <- paste(min(config$alliance_years), max(config$alliance_years), sep = "-")
patent_range   <- paste(min(config$patent_years), max(config$patent_years), sep = "-")

p <- ggraph(g_tidy, layout = "fr") +
  geom_edge_link(alpha = 0.3, color = "gray60") +
  geom_node_point(
    aes(color = color_type, shape = color_type, size = node_size),
    alpha = 0.85
  ) +
  geom_node_text(
    aes(label = ifelse(is_top_uni == TRUE, str_trunc(name, 30), "")),
    color = "#B8860B", size = 2.5, repel = TRUE, 
    max.overlaps = 30, show.legend = FALSE
  ) +
  scale_color_manual(values = node_colors, name = "Organization Type") +
  scale_shape_manual(values = node_shapes, name = "Organization Type") +
  scale_size_continuous(
    name = "Patents",
    range = c(2, 20),
    breaks = c(3, 3 + 1.5*sqrt(10), 3 + 1.5*sqrt(100), 3 + 1.5*sqrt(500)),
    labels = c("0", "10", "100", "500")
  ) +
  labs(
    title = "Pharma-University Alliance Network",
    subtitle = glue::glue(
      "Pharma: UK + Germany + France | Unis: Global | ",
      "Alliances: {alliance_range} | Patents: {patent_range} | N={vcount(g)}\n",
      "Pharma+Uni (green \u25cf): {n_pharma_with} | ",
      "Pharma only (red \u25cf): {n_pharma_without} | ",
      "Unis (blue \u25b2): {n_unis - n_top_unis} | ",
      "Top Unis (gold \u25b2): {n_top_unis}"
    ),
    caption = glue::glue(
      "Node size = patent count (sqrt scaled). ",
      "Gold = Top {config$top_uni_rank} CWUR ranked universities."
    )
  ) +
  theme_void() +
  theme(
    legend.position  = "right",
    legend.box       = "vertical",
    plot.title       = element_text(size = 14, face = "bold"),
    plot.subtitle    = element_text(size = 9, color = "gray40"),
    plot.caption     = element_text(size = 8, color = "gray50"),
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9)
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 5)),
    shape = guide_legend(order = 1, override.aes = list(size = 5)),
    size  = guide_legend(order = 2)
  )

# ─────────────────────────────────────────────────────────────────────────────
# Save output
# ─────────────────────────────────────────────────────────────────────────────

print(p)

ggsave(
  file.path(config$output_dir, "network_plot.png"),
  plot = p, width = 16, height = 14, dpi = 300, bg = "white"
)

cat("  Saved: network_plot.png\n")