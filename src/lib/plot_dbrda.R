#' Plot multiple constrained ordinations
plot_dbrda <- function(kingdom, samples, lineages, constrained_ordination, anova_all, anova_terms, taxrank = "genus") {
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
    filter(taxon %in% sig_terms) %>%
    left_join(lineages %>% select(-taxon), by = c("taxon" = taxrank)) %>%
    rename(dbRDA1 = CAP1, dbRDA2 = CAP2) %>%
    mutate(
      length = sqrt(dbRDA1^2 + dbRDA2^2),
      phylum_color = phylum %>%
        {
          ifelse(. %in% names(phyla_colors), ., "other")
        } %>% factor(levels = names(phyla_colors))
    ) %>%
    arrange(-length)

  dbrda_annot_text <- str_glue(
    paste(
      "p = {anova_all$`Pr(>F)`[[1]] %>% sprintf(fmt = '%.2e')}",
      "R² = {dbrda %>% RsquareAdj() %>% pluck('r.squared') %*% 100 %>% sprintf(fmt = '%.2f')}%",
      #"{length(terms_selected)}/{length(terms_all)} selected",
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

  plt <-
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
      data = dbrda_terms_tbl %>% filter(!is.na(phylum_color)) %>% select(dbRDA1, dbRDA2, phylum_color),
      arrow = arrow(length = unit(0.1, "cm"), type = "closed"),
      mapping = aes(x = 0, y = 0, xend = dbRDA1 * 2, yend = dbRDA2 * 2, color = phylum_color)
    ) +
    scale_color_phyla() +
    annotate("text", -Inf, Inf, size = 3, label = dbrda_annot_text, hjust = "inward", vjust = "inward") +
    facet_wrap(~ str_glue("{kingdom}")) +
    labs(
      x = (dbrda %>% proportion_explained("CAP1") * 100) %>% sprintf(fmt = "dbRDA1 (%.1f%%)"),
      y = (dbrda %>% proportion_explained("CAP2") * 100) %>% sprintf(fmt = "dbRDA2 (%.1f%%)"),
      color = "phylum"
    )

  plt
}
