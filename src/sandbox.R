source("src/defaults.R")

#
# dbRDA constrained by generalists and specialists --------
#

loadd(selected_generalists_specialists)
loadd(prevalences, distances, sub_abundances)

constrained_ordinations <-
  distances %>%
  inner_join(prevalences %>% rename(prevalence = data)) %>%
  inner_join(sub_abundances %>% rename(abundance = data)) %>%
  filter(norm_method == "tss" & dist_method == "bray" & dist_type == "samples") %>%
  filter(subset_value %>% str_detect("host,Fungi")) %>%
  filter(subset_name %>% str_detect("kingdom")) %>%
  separate(subset_value, sep = ",", into = c("subset", "kingdom")) %>%
  transmute(
    subset_name, subset, kingdom,
    ord = list(dist, prevalence, abundance) %>% pmap(possibly(~ {
      dist <- .x
      meta <-
        ..3 %>%
        select(sample_id, taxon, abundance) %>%
        pivot_wider(names_from = taxon, values_from = abundance, values_fill = list(abundance = 0))

      start_vars <- colnames(meta) %>% setdiff("sample_id")

      min_constrained_dbrda <- vegan::capscale(
        formula = dist ~ 1,
        data = meta %>% ungroup() %>% select(-sample_id)
      )

      max_constrained_dbrda <- vegan::capscale(
        formula = start_vars %>% paste0(collapse = "` + `") %>% paste0("dist ~ `", ., "`") %>% as.formula(),
        data = meta %>% ungroup() %>% select(-sample_id)
      )

      selected_dbrda <-
        vegan::capscale(
          formula = {
            selected_generalists_specialists %>%
              semi_join(..3) %>%
              pull(taxon) %>%
              paste0(collapse = "` + `") %>%
              paste0("dist ~ `", ., "`") %>%
              as.formula()
          },
          data = meta %>% ungroup() %>% select(-sample_id)
        )

      list(
        min_constrained_dbrda = min_constrained_dbrda,
        max_constrained_dbrda = max_constrained_dbrda,
        selected_dbrda = selected_dbrda
      )
    }, NA)),
    anova_all = ord %>% map(~ .x$selected_dbrda %>% anova(permutations = 1e2, parallel = getOption("mc.cores"))),
    anova_terms = ord %>% map(~ .x$selected_dbrda %>% anova(by = "terms", permutations = 1e2, parallel = getOption("mc.cores")))
  )


constrained_ordinations_plots <-
  constrained_ordinations %>%
  transmute(
    subset_name, subset,
    plot = list(ord, anova_all, anova_terms, kingdom) %>% pmap(~ {
      constrained_ordination <- ..1
      anova_all <- ..2
      anova_terms <- ..3
      kingdom <- ..4
      taxrank <- "genus"

      dbrda <- constrained_ordination$selected_dbrda
      terms_selected <- dbrda %>%
        as.formula() %>%
        attr("term.labels")
      terms_all <- constrained_ordination$max_constrained_dbrda %>%
        as.formula() %>%
        attr("term.labels")

      dbrda_sites_tbl <-
        dbrda %>%
        plot() %>%
        pluck("sites") %>%
        as_tibble(rownames = "sample_id") %>%
        rename(dbRDA1 = CAP1, dbRDA2 = CAP2)

      sig_terms <-
        anova_terms %>%
        tidy() %>%
        ungroup() %>%
        mutate(q.value = p.adjust(p.value)) %>%
        filter(q.value < 0.05) %>%
        pull(term)

      dbrda_terms_tbl <-
        dbrda %>%
        plot() %>%
        pluck("biplot") %>%
        as_tibble(rownames = "taxon") %>%
        rename(dbRDA1 = CAP1, dbRDA2 = CAP2) %>%
        filter(taxon %in% sig_terms) %>%
        left_join(selected_generalists_specialists)

      dbrda_annot_text <- str_glue(
        paste(
          "p = {anova_all$`Pr(>F)`[[1]] %>% sprintf(fmt = '%.2e')}",
          "R² = {dbrda %>% RsquareAdj() %>% pluck('r.squared') %*% 100 %>% sprintf(fmt = '%.2f')}%",
          # "{length(terms_selected)}/{length(terms_all)} selected",
          sep = "\n"
        )
      )

      proportion_explained <- function(dbrda, axis = "CAP1") {
        dbrda %>%
          summary() %>%
          pluck("concont", "importance") %>%
          as_tibble(rownames = "name") %>%
          filter(name == "Proportion Explained") %>%
          pluck(axis)
      }

      tibble() %>%
        ggplot(aes(dbRDA1, dbRDA2)) +
        geom_point(
          data = dbrda_sites_tbl %>% left_join(samples),
          mapping = aes(color = environment_group),
          alpha = 0.2
        ) +
        labs(color = "Environment") +
        scale_color_environment_group() +
        # Do not show alpha and small size of sample point in legend
        guides(color_new = guide_legend(override.aes = list(alpha = 1, size = 4))) +
        ggnewscale::new_scale_color() +
        geom_segment(
          # do not select kingdom here. This will overwrite facet
          data = dbrda_terms_tbl %>% select(dbRDA1, dbRDA2, prevalence_group),
          arrow = arrow(length = unit(0.1, "cm"), type = "closed"),
          mapping = aes(x = 0, y = 0, xend = dbRDA1 * 2, yend = dbRDA2 * 2, color = prevalence_group)
        ) +
        scale_color_prevalence_group() +
        annotate("text", -Inf, Inf, size = 3, label = dbrda_annot_text, hjust = "inward", vjust = "inward") +
        facet_wrap(~ str_glue("{kingdom}")) +
        labs(
          x = (dbrda %>% proportion_explained("CAP1") * 100) %>% sprintf(fmt = "dbRDA1 (%.1f%%)"),
          y = (dbrda %>% proportion_explained("CAP2") * 100) %>% sprintf(fmt = "dbRDA2 (%.1f%%)"),
          color = "Prevalence group"
        )
    })
  )

constrained_ordinations_plots %>%
  pull(plot) %>%
  wrap_plots()


#
# Shared genera -----
#
# Are there really more shared bacteria between soil and aquatic genera than
# there are genera shared between host and aquatic? Is this influenced by sample size?
#

# down sampling
n_samples <-
  sub_abundances$data[[1]] %>%
  filter(abundance > 0 & kingdom %in% kingdoms_groups & environment_group %in% environment_groups) %>%
  distinct(sample_id, environment_group) %>%
  count(environment_group) %>%
  pull(n) %>%
  min() %>%
  `*`(0.5) %>% # allow randomness even in the smallest subset
  as.integer()

tibble(trail = seq(1000)) %>%
  mutate(
    data = trail %>% map(~ {
      set.seed(.x)
      
      sub_abundances$data[[1]] %>%
        filter(
          abundance > 0 &
            kingdom %in% kingdoms_groups &
            environment_group %in% environment_groups
        ) %>%
        # sub sampling  
        nest(-c(environment_group, sample_id)) %>%
        group_by(environment_group) %>%
        sample_n(n_samples) %>%
        unnest(data) %>%
        distinct(environment_group, kingdom, taxon) %>%
        mutate(is_prevalent = TRUE) %>%
        pivot_wider(
          names_from = environment_group,
          values_from = is_prevalent,
          values_fill = list(is_prevalent = FALSE)
        ) %>%
        distinct(taxon, .keep_all = TRUE) %>%
        group_by(kingdom) %>%
        summarise(
          soil_aquatic = sum(soil & aquatic),
          soil_host = sum(soil & host),
          host_aquatic = sum(host & aquatic)
        ) %>%
        pivot_longer(c(soil_aquatic, soil_host, host_aquatic))
    })
  ) %>%
  unnest() %>%
  group_by(trail, kingdom) %>%
  arrange(-value) %>%
  summarise(sample_subset_order = paste0(name, collapse = ",")) %>%
  group_by(kingdom, sample_subset_order) %>%
  count() %>%
  group_by(kingdom) %>%
  mutate(freq = n / sum(n))

#
# Alpha diversity much lower in samples without generalits compared to samples without random taxa ? -----
#

# 6) Daniel we ned to get some additional evidence that the arbitary cut off of >40% in three biomes for a generalist in
# not critical for the main conclusions, such as the alpha diversity comparison with vs without.
# Generate the same statistic for >30% and >50% and tell us if it remains significant.

get_alpha_div <- function(permute = FALSE, seed = 1, min_prevalence_perc = 40) {
  selected_generalists <-
    common_generalists %>%
    filter(min_prevalence_perc == !!min_prevalence_perc) %>%
    pull(taxon) %>%
    unique()
  
  if (permute) {
    set.seed(seed)
    
    selected_generalists <-
      sub_abundances$data[[2]]$taxon %>%
      unique() %>%
      sample(length(selected_generalists))
  }
  
  samples_having_any_generalists <-
    sub_abundances %>%
    filter(subset_name == "kingdom" & norm_method == "tss") %>%
    transmute(kingdom = subset_value, data) %>%
    select(-kingdom) %>%
    unnest(data) %>%
    mutate(has_generalist = taxon %in% selected_generalists) %>%
    filter(abundance > 0) %>%
    group_by(sample_id, kingdom) %>%
    summarise(has_generalist = any(has_generalist))
  
  samples_having_any_generalists %>%
    inner_join(
      alphadiv %>%
        filter(subset_name == "kingdom") %>%
        select(kingdom = subset_value, data) %>%
        unnest(data) %>%
        select(-Simpson)
    ) %>%
    pivot_longer(cols = c(Chao1, Shannon), names_to = "alphadiv_metric", values_to = "alphadiv") %>%
    group_by(kingdom, has_generalist, alphadiv_metric) %>%
    summarise(median_alpha_div = sum(value))
}


expected <-
  tibble(trail = seq(100)) %>%
  expand_grid(min_prevalence_perc = c(40))

expected$data <- pbmcmapply(
  get_alpha_div,
  permute = TRUE,
  seed = expected$trail,
  min_prevalence_perc = expected$min_prevalence_perc,
  mc.cores = 10
)

expected <-
  expected %>%
  unnest() %>%
  rename(expected = alphadiv)

truth <-
  expected %>%
  distinct(min_prevalence_perc) %>%
  mutate(truth = min_prevalence_perc %>% map(~get_alpha_div(permute = FALSE, min_prevalence_perc = .))) %>%
  unnest(truth) %>%
  unnest(truth) %>%
  rename(truth = alphadiv)


expected %>%
  inner_join(truth) %>%
  group_by(min_prevalence_perc, kingdom, alphadiv_metric, has_generalist) %>%
  #filter(! has_generalist) %>%
  mutate(false_positive = expected <= truth) %>%
  summarise(p.value = sum(false_positive) / n())

expected %>%
  inner_join(truth) %>%
  group_by(min_prevalence_perc, kingdom, alphadiv_metric, has_generalist) %>%
  #filter(! has_generalist) %>%
  mutate(false_positive = expected < truth) %>%
  summarise(p.value = sum(false_positive) / n())


expected %>%
  full_join(truth, c("min_prevalence_perc", "kingdom", "has_generalist", "alphadiv_metric")) %>%
  filter(min_prevalence_perc == 40) %>%
  pivot_longer(c(expected, truth)) %>%
  ggplot(aes(name, value, color = name)) +
    geom_boxplot() +
    facet_wrap(kingdom ~ has_generalist + alphadiv_metric, scales = "free") +
    labs(x = "samples with specific taxa", y = "Median alpha diversity")


#
# network topology -----
#

min_abs_correlation <- 0.2

edges <-
  graphs %>%
  filter(subset_name %in% c("all", "environment_group") & cor_method == "sparcc") %>%
  transmute(
    subset_value,
    data = cor_graph %>% map(~ {
      .x %>%
        activate(edges) %>%
        # treat small correlation e.g. +0.001 and -0.001 the same
        # to reduce noise
        filter(q.value < 0.05 & abs(estimate) > min_abs_correlation) %>%
        as_tibble() %>%
        left_join(selected_generalists_specialists %>% rename_all(~ paste0("from_", .))) %>%
        left_join(selected_generalists_specialists %>% rename_all(~ paste0("to_", .))) %>%
        select(contains("kingdom"), contains("taxon"), estimate, contains("prevalence_group"))
    })
  ) %>%
  unnest() %>%
  mutate(
    kingdom = case_when(
      from_kingdom == "Bacteria" & to_kingdom == "Bacteria" ~ "intra bacteria",
      from_kingdom == "Fungi" & to_kingdom == "Fungi" ~ "intra fungi",
      from_kingdom != to_kingdom ~ "inter kingdom"
    ),
    prevalence_group = case_when(
      ! is.na(from_prevalence_group) ~ from_prevalence_group,
      ! is.na(to_prevalence_group) ~ to_prevalence_group,
      TRUE ~ "other",
    ) %>%
      recode(
        "generalist" = "with generalists",
        "specialist" = "with specialists"
      ),
    edge_group = paste(kingdom, prevalence_group, sep = "\n"),
    edge_sub_group = paste0(sign(estimate), prevalence_group, kingdom),
    special_taxon = case_when(
      ! is.na(from_prevalence_group) ~ from_taxon,
      ! is.na(to_prevalence_group) ~ to_taxon,
      TRUE ~ from_taxon
    ),
    other_taxon = case_when(
      to_taxon == special_taxon ~ from_taxon,
      TRUE ~ to_taxon
    )
  )

edges %>%
  ggplot(aes(edge_group, estimate)) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = - min_abs_correlation, ymax = min_abs_correlation, fill = "lightgrey") +
  geom_hline(yintercept = 0) +
  geom_tile(
    # show all combinations
    data = edges %>% distinct(edge_group, edge_sub_group, kingdom),
    mapping = aes(fill = kingdom, y = -1),
    height = 0.2
  ) +
  scale_fill_aaas() +
  scale_fill_manual(values = c(
    "intra bacteria" =  as.character(kingdoms_colors["Bacteria"]),
    "inter kingdom" = "black",
    "intra fungi" =   as.character(kingdoms_colors["Fungi"])
  )) +
  geom_boxplot(
    mapping = aes(color = prevalence_group, group = edge_sub_group),
    position = position_dodge(0),
    varwidth = TRUE
  ) +
  scale_color_manual(values = c(
    "with generalists" = as.character(prevalence_group_colors["generalist"]),
    "with specialists" = as.character(prevalence_group_colors["specialist"]),
    "other" = "black"
  )
  ) +
  stat_compare_means(
    method = "wilcox",
    comparisons = list(
      c("intra bacteria\nwith generalists", "intra bacteria\nother"),
      c("intra fungi\nwith generalists", "intra fungi\nother"),
      c("inter kingdom\nwith generalists", "inter kingdom\nother"),
      
      c("intra bacteria\nwith generalists", "intra fungi\nwith generalists"),
      c("intra bacteria\nwith generalists", "inter kingdom\nwith generalists")
    )
  ) +
  scale_y_continuous(expand = expansion(c(0, 0.3))) +
  scale_x_discrete(expand = c(0, 0)) +
  facet_wrap(~ subset_value, ncol = 1) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_blank()
  ) +
  labs(
    x = "Edge group",
    y = str_glue("SparCC correlation ($|r| > {min_abs_correlation}$)") %>% latex2exp::TeX(),
    color = "",
    fill = ""
  )

ggsave("results/correlation_pos_neg_prevalence_group.png", width = 6, height = 9)


edges %>%
  mutate(
    edge_group = edge_group %>% str_replace_all("[\n]", " "),
    edge_direction = ifelse(estimate > 0, "positive", "negative")
  ) %>%
  ggplot(aes(edge_group, fill = edge_direction)) +
    geom_bar(position = position_dodge()) +
    geom_text(stat="count", aes(label=..count.., color = edge_direction), position = position_dodge(1), hjust=-0.1) +
    scale_fill_aaas() +
    scale_color_aaas() +
    coord_flip() +
    facet_wrap(~ subset_value, ncol = 1, scales = "free_y") +
    scale_y_continuous(expand = expansion(c(0, 0.3))) +
    labs(
      x = "Edge group",
      y = "Edges"
    )

ggsave("results/correlation_pos_counts.png", width = 6, height = 9)


nodes <-
  graphs %>%
  filter(subset_name %in% c("all", "environment_group") & cor_method == "sparcc") %>%
  transmute(
    subset_value,
    data = cor_graph %>% map(~ {
      .x %>%
        activate(edges) %>%
        # treat small correlation e.g. +0.001 and -0.001 the same
        # to reduce noise
        filter(q.value < 0.05 & abs(estimate) > min_abs_correlation) %>%
        activate(nodes) %>%
        # filtering edges requires to re-calculate topology
        topologize_graph() %>%
        left_join(selected_generalists_specialists) %>%
        as_tibble()
    })
  ) %>%
  unnest() %>%
  mutate(
    prevalence_group = prevalence_group %>% replace_na("other"),
    taxon_group = paste(prevalence_group, kingdom, sep = "\n")
  ) %>%
  pivot_longer(c(degree, closeness, betweeness, hub), names_to = "node_topology_name", values_to = "node_topology_value")

nodes %>%
  nest(-c(subset_value, node_topology_name)) %>%
  mutate(
    log_transform = node_topology_name %in% c("betweeness", "degree"),
    plt = list(subset_value, node_topology_name, log_transform, data) %>% pmap(~ {
      plt <-
        ..4 %>%
        mutate(
          facet_x = ..1,
          facet_y = ..2,
          # keep all levels to combine guide
          prevalence_group = prevalence_group %>% factor(levels = prevalence_group_colors %>% names() %>% c("other"))
        ) %>%
        ggplot(aes(taxon_group, node_topology_value)) +
        geom_boxplot(aes(color = prevalence_group)) +
        scale_color_manual(values = prevalence_group_colors %>% append(c("other" = "black")), drop = FALSE) +
        stat_compare_means(
          method = "wilcox",
          comparisons = {
            # compare only those who are available
            # e.g. specialists might fail in some plots
            ..4$taxon_group %>%
              unique() %>%
              combn(2) %>%
              t() %>%
              as_tibble() %>%
              filter(V1 %>% str_detect("generalist") | V2 %>% str_detect("generalist")) %>%
              transmute(comp = list(V1, V2) %>% pmap(~ c(.x, .y))) %>%
              pull(comp)
          }
        ) +
        theme(
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()
        ) +
        facet_wrap(~ facet_x, scale = "free") +
        labs(x = "", color = "Taxon group", y = ..2)
    
      if(..3) {
        plt <-
          plt +
          scale_y_log10(expand = expansion(c(0, 0.1))) +
          annotation_logticks(sides = "l") 
      }
      
      plt
    })
  ) %>%
  pull(plt) %>%
  wrap_plots(guides = "collect", ncol = 4)


#
# Network of common correlations -----
#

node_majority_edge_group <- function(taxon, min_purity_frac = 0.5) {
  graph %>%
    activate(nodes) %>%
    filter(node_is_adjacent(taxon == {{taxon}})) %>%
    activate(edges) %>%
    as_tibble() %>%
    count(environment_group) %>%
    arrange(-n) %>%
    mutate(n = n / sum(n)) %>%
    filter(n > min_purity_frac) %>%
    pull(environment_group) %>%
    first()
}

taxa_prevalent_all_environment_groups <-
  prevalences %>%
  filter(subset_name == "environment_group") %>%
  unnest(data) %>%
  filter(prevalence_n > 0) %>%
  select(subset_value, taxon, prevalence_n) %>%
  pivot_wider(names_from = subset_value, values_from = prevalence_n) %>%
  filter(!is.na(host) & !is.na(aquatic) & !is.na(soil)) %>%
  pull(taxon)

graph <-
  graphs %>%
  filter(subset_name == "environment_group") %>%
  rename(environment_group = subset_value) %>%
  filter(cor_method == "sparcc" & norm_method == "raw") %>%
  mutate(
    cor_graph = cor_graph %>% map(~ .x %>% activate(edges) %>% as_tibble())
  ) %>%
  unnest(cor_graph) %>%
  select(from_taxon, to_taxon, environment_group, estimate, q.value) %>%
  filter(
    q.value < 0.05 &
      abs(estimate) > generalists_min_abs_correlation &
      from_taxon %in% taxa_prevalent_all_environment_groups &
      to_taxon %in% taxa_prevalent_all_environment_groups
  ) %>%
  mutate(estimate_direction = estimate %>% map_chr(~ifelse(sign(.x) == 1, "positive", "negative"))) %>%
  group_by(from_taxon, to_taxon, estimate_direction) %>%
  arrange(from_taxon, to_taxon) %>%
  filter(n() >= 2) %>% # present in all environments
  pivot_wider(names_from = environment_group, values_from = c(estimate, q.value)) %>%
  tbl_graph(edges = .) %>%
  activate(nodes) %>%
  mutate(taxon = name) %>%
  left_join(lineages) %>%
  left_join(selected_generalists) %>%
  annotate_node_attributes_in_edges() %>%
  activate(edges) %>%
  mutate(
    environment_group = case_when(
      ! is.na(estimate_host) & ! is.na(estimate_aquatic) & ! is.na(estimate_soil) ~ "host-aquatic-soil",
      ! is.na(estimate_host) & ! is.na(estimate_aquatic) ~ "host-aquatic" ,
      ! is.na(estimate_host) & ! is.na(estimate_soil) ~ "host-soil" ,
      ! is.na(estimate_aquatic) & ! is.na(estimate_soil) ~ "aquatic-soil" ,
      ! is.na(estimate_host) & ! is.na(estimate_soil) ~ "host-soil"
    ),
    kingdoms_group = case_when(
      from_kingdom == "Bacteria" & to_kingdom == "Bacteria" ~ "intra Bacteria",
      from_kingdom == "Fungi" & to_kingdom == "Fungi" ~ "intra Fungi",
      from_kingdom != to_kingdom ~ "inter kingdom"
    ),
    involves_prevalence_group_taxon = ! is.na(from_prevalence_group) | ! is.na(to_prevalence_group)
  ) %>%
  activate(nodes) %>%
  mutate(
    majority_edge_group = taxon %>% map_chr(~node_majority_edge_group(.x, 0)),
    taxon_group = paste0(majority_edge_group, kingdom)
  ) %>%
  topologize_graph() %>%
  annotate_node_attributes_in_edges()

graph %>%
  ggraph(layout = layout_in_circle(graph, order=order(V(graph)$taxon_group))) +
  geom_edge_link(aes(color = environment_group, alpha = involves_prevalence_group_taxon)) +
  geom_node_point(aes(color = kingdom, size = prevalence_group)) +
  scale_edge_color_manual(values = environment_group_sets_colors) +
  scale_color_kingdom(na.value = "black") +
  scale_size_manual(values = c("generalist" = 5), na.value = 2, guide = "none") +
  scale_edge_alpha_manual(values = c(`TRUE` = 0.5, `FALSE` = 0.05), guide = "none") +
  coord_fixed() +
  facet_edges(~estimate_direction) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  labs(
    color = "Kingdom",
    edge_color = "Common in biomes"
  )

ggsave(paste0(results_dir, "net.png"), width = 15, height = 8, dpi = 600, bg = "white")
