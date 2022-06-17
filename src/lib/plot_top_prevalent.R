plot_top_prevalent <- function(prevalence_summary, lineages) {
  tbl <-
    prevalence_summary %>%
    filter(subset_name == "kingdom") %>%
    ungroup() %>%
    select(kingdom = subset_value, data) %>%
    unnest(data) %>%
    select(taxon, group, prevalence_perc)

  # Add missing zeros
  tbl <-
    expand_grid(taxon = tbl$taxon %>% unique(), group = tbl$group %>% unique(), prevalence_perc = 0) %>%
    full_join(tbl) %>%
    group_by(taxon, group) %>%
    arrange(-prevalence_perc) %>%
    distinct() %>%
    arrange(taxon, group) %>%
    group_by(taxon, group) %>%
    arrange(-prevalence_perc) %>%
    slice(1) %>%
    ungroup()

  tbl %<>% left_join(lineages)

  top_prevalent_tbl <-
    tbl %>%
    group_by(kingdom, group) %>%
    arrange(-prevalence_perc) %>%
    slice(1:10) %>%
    select(kingdom, taxon, group) %>%
    mutate(is_top_prevalent = TRUE)

  top_prevalent_taxa <-
    top_prevalent_tbl %>%
    pull(taxon) %>%
    unique()

  span_tbl <-
    tbl %>%
    filter(taxon %in% top_prevalent_taxa) %>%
    group_by(taxon) %>%
    summarise(min = min(prevalence_perc), max = max(prevalence_perc))

  tbl <-
    tbl %>%
    left_join(span_tbl) %>%
    filter(taxon %in% top_prevalent_taxa) %>%
    mutate(phylum = phylum %>% {
      ifelse(. %in% names(phyla_colors), ., "other")
    }) %>%
    left_join(top_prevalent_tbl) %>%
    replace_na(list(is_top_prevalent = FALSE)) %>%

    # arrange taxa
    arrange(kingdom, phylum, taxon) %>%
    mutate(
      taxon = taxon %>% factor(levels = {
        lineages %>%
          filter(taxon %in% top_prevalent_taxa) %>%
          arrange(kingdom, phylum, class, order, family, taxon) %>%
          pull(taxon)
      }),
      phylum = phylum %>% factor(levels = names(phyla_colors))
    ) %>%
    
    group_by(group) %>%
    mutate(rank_top_prevalence = rank(prevalence_perc)) %>%
    ungroup() %>%

    # plot
    rename(Phylum = phylum)


  tbl %>%
    ggplot(aes(x = taxon, y = prevalence_perc)) +
    geom_tile(aes(y = -3.5, fill = Phylum), height = 5) +
    scale_fill_phyla() +
    geom_segment(
      data = tbl %>% group_by(taxon) %>% slice(1),
      mapping = aes(y = min, yend = max, xend = taxon)
    ) +
    geom_point(aes(color = group), stat = "identity", position = position_dodge(width = 0.5), size = 3) +
    scale_color_environment_group() +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    coord_flip(clip = "off") + # do not truncate edge points
    facet_wrap(~kingdom, nrow = 1, scales = "free_y") +
    theme_pub() +
    theme(
      panel.grid.major.x = element_line(),
      axis.text.y = element_text(face = "italic")
    ) +
    labs(
      y = "Prevalence (% per environment)",
      x = "Genus",
      shape = "Top prevalent",
      color = "Environment"
    )
}
