source("src/03-generalists-specialists-genome-with-gff/defaults.R")

tar_load(figure1)
tar_load(figure4)

library(tidygraph)
library(ggraph)
library(igraph)

generalists_common_coabundance_graph <- readRDS("cache/15-microbiome-rserver-hki-rdata/generalists_common_coabundance_graph.rds")

kingdoms_colors <- c(
  "Bacteria" = "#3f78a9",
  "Fungi" = "#A93F55",
  "all" = "black"
)


scale_color_kingdom <- function(drop = FALSE, na.value = "#d9d9d9", ...) {
  scale_color_manual(
    values =  purrr::simplify(kingdoms_colors),
    na.value = na.value,
    drop = drop,
    ...
  )
}

environment_group_colors <- c(
  "soil" = "#8b786d",
  "aquatic" = "#78a1bb",
  "host" = "#e0aa5a",
  "all" = "black"
)

environment_group_sets_colors <- c(
  "host-aquatic-soil" = "#000000",
  "aquatic-soil" = "#155A82",
  "host-aquatic" = "#83155F",
  "host-soil" = "#826516"
) %>%
  append(environment_group_colors)


d <- 
  generalists_common_coabundance_graph %>%
  activate(edges) %>%
  mutate(estimate_direction = estimate_direction %>% recode("positive" = "Positive", "negative" = "Negative")) %>%
  ggraph(layout = layout_in_circle(generalists_common_coabundance_graph, order = order(V(generalists_common_coabundance_graph)$taxon_group))) +
  geom_edge_link(aes(color = environment_group, alpha = involves_prevalence_group_taxon)) +
  geom_node_point(aes(color = kingdom, size = is_generalist)) +
  scale_edge_color_manual(values = environment_group_sets_colors, labels = str_to_sentence) +
  scale_color_kingdom(na.value = "black") +
  scale_size_manual(values = c(`TRUE` = 5, `FALSE` = 2), labels = c("No", "Yes")) +
  scale_edge_alpha_manual(values = c(`TRUE` = 0.5, `FALSE` = 0.05), guide = "none") +
  coord_fixed() +
  facet_edges(~estimate_direction) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_void() +
  labs(
    color = "Kingdom",
    edge_color = "Common in environments",
    size = "Generalist")
################################################################################
kingdoms_groups <- readRDS("cache/15-microbiome-rserver-hki-rdata/kingdoms_groups.rds")
sub_abundances <- readRDS("cache/15-microbiome-rserver-hki-rdata/sub_abundances.rds")
selected_generalists_specialists <- readRDS("cache/15-microbiome-rserver-hki-rdata/selected_generalists_specialists.rds")

prevalence_group_colors <- c(
  "generalist" = "#3a86ff",
  "specialist" = "#ffbe0b"
)


scale_color_prevalence_group <- function(drop = FALSE, ...) {
  scale_color_manual(
    values =  purrr::simplify(prevalence_group_colors),
    na.value = "#000000",
    drop = drop,
    ...
  )
}

symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                    symbols = c("****", "***", "**", "*", "ns"))

bc <- 
  kingdoms_groups %>%
  set_names(kingdoms_groups) %>%
  purrr::map(function(current_kingdom) {
    sub_abundances %>%
      filter(subset_value == "all" & norm_method == "tss") %>%
      pull(data) %>%
      first() %>%
      inner_join(selected_generalists_specialists) %>%
      filter(kingdom == current_kingdom & !is.na(prevalence_group)) %>%
      filter(abundance > 0.01e-2) %>%
      ggplot(aes(x = prevalence_group, y = abundance, color = prevalence_group)) +
      geom_boxplot(width = 0.5) +
      facet_grid(~kingdom) +
      scale_x_discrete(labels = str_to_sentence) +
      scale_y_log10(expand = c(0, 1)) +
      stat_compare_means(
        method = "wilcox",
        symnum.args = symnum.args,
        comparisons = list(c("generalist", "specialist"))
      ) +
      annotation_logticks(sides = "l") +
      scale_color_prevalence_group() +
      theme_bw() +
      theme(
        strip.text = element_text(size = rel(1.1)),
        axis.text = element_text(size = rel(1.1)), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        axis.ticks.x=element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
      guides(color = FALSE) +
      labs(
        x = "",
        y = "Abundance (TSS)"
      )
  })

################################################################################

layout <- "
AABBCC
DDEEEE
"

f <- wrap_plots(
  wrap_plots(
    bc[[1]] + labs(tag = "A"),
    bc[[2]]
  ),
  figure1[[3]] + labs(tag = "B"),
  figure1[[1]] + labs(tag = "C"),
  figure4 + labs(tag = "D"),
  d + labs(tag = "E"),
  ncol = 3,
  nrow = 2,
  guides = "collect") + 
  plot_layout(design = layout
  # plot_layout(widths = c(1,1,2), heights = c(0.5, 2)
  )

ggsave(filename = "plots/04-manuscript-ap/fig2-v1.3.pdf",
       plot = f,
       device = "pdf",
       width = 45,
       height = 25,
       units = "cm",
       dpi = 300)

ggsave(filename = "plots/04-manuscript-ap/fig2-v1.3.png",
       plot = f,
       device = "png",
       width = 45,
       height = 25,
       units = "cm",
       dpi = 300)

# BIG CHANGES
##############################################################################

source("src/07-amr/defaults.R")
tar_load(amr_comparisons_plot)

amr_comparisons_plot[1]

layout <- "
AABBCC
DDEEEE
"

f <- wrap_plots(
  wrap_plots(
    bc[[1]] + labs(tag = "A"),
    bc[[2]]
  ),
  figure1[[3]] + labs(tag = "B"),
  figure1[[1]] + labs(tag = "C"),
  amr_comparisons_plot[[1]] + labs(tag = "D"),
  d + labs(tag = "E"),
  ncol = 3,
  nrow = 2,
  guides = "collect") + 
  plot_layout(design = layout
              # plot_layout(widths = c(1,1,2), heights = c(0.5, 2)
  )

ggsave(filename = "plots/04-manuscript-ap/fig2-v1.4.pdf",
       plot = f,
       device = "pdf",
       width = 45,
       height = 25,
       units = "cm",
       dpi = 300)

ggsave(filename = "plots/04-manuscript-ap/fig2-v1.4.png",
       plot = f,
       device = "png",
       width = 45,
       height = 25,
       units = "cm",
       dpi = 300)
