#' Plot multiple constrained ordinations
plot_dbrdas <- function(samples, lineages, stepwise_constrained_ordinations, stepwise_constrained_ordinations_anova, taxrank = "genus") {
  tbl <-
    stepwise_constrained_ordinations %>%
    inner_join(stepwise_constrained_ordinations_anova) %>%
    mutate(
      terms_selected = ord %>% map(~ .x$selected_dbrda %>%
        as.formula() %>%
        attr("term.labels")),
      terms_all = ord %>% map(~ .x$max_constrained_dbrda %>%
        as.formula() %>%
        attr("term.labels")),
      dbrda_sites_tbl = ord %>% map(~ {
        .x$selected_dbrda %>%
          plot() %>%
          pluck("sites") %>%
          as_tibble(rownames = "sample_id") %>%
          rename(dbRDA1 = CAP1, dbRDA2 = CAP2)
      }),

      sig_terms = anova_terms %>% map(~ {
        .x %>%
          tidy() %>%
          ungroup() %>%
          mutate(q.value = p.adjust(p.value)) %>%
          filter(q.value < 0.05) %>%
          pull(term)
      }),

      dbrda_terms_tbl = ord %>% map2(sig_terms, ~ {
        .x$selected_dbrda %>%
          plot() %>%
          pluck("biplot") %>%
          as_tibble(rownames = "taxon") %>%
          filter(taxon %in% .y) %>%
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
      }),

      dbrda_annot_text = list(anova_all, ord, terms_selected, terms_all) %>% pmap(function(anova_all, ord, terms_selected, terms_all) {
        str_glue(
          paste(
            "p = {anova_all$`Pr(>F)`[[1]] %>% sprintf(fmt = '%.2e')}",
            "R² = {ord$selected_dbrda %>% RsquareAdj() %>% pluck('r.squared') %*% 100 %>% sprintf(fmt = '%.2f')}%",
            "{length(terms_selected)}/{length(terms_all)} selected",
            sep = "\n"
          )
        )
      }),

      proportion_explained_axis_1 = ord %>% map_dbl(~ {
        .x$selected_dbrda %>%
          summary() %>%
          pluck("concont", "importance") %>%
          as_tibble(rownames = "name") %>%
          filter(name == "Proportion Explained") %>%
          pull(CAP1)
      }),

      proportion_explained_axis_2 = ord %>% map_dbl(~ {
        .x$selected_dbrda %>%
          summary() %>%
          pluck("concont", "importance") %>%
          as_tibble(rownames = "name") %>%
          filter(name == "Proportion Explained") %>%
          pull(CAP2)
      }),
    ) %>%
    group_by(kingdom)

  plt <-
    tibble() %>%
    ggplot(aes(dbRDA1, dbRDA2)) +

    # sample points
    geom_point(
      data = tbl %>% select(dbrda_sites_tbl) %>% unnest(dbrda_sites_tbl) %>% left_join(samples),
      mapping = aes(color = environment_group)
    ) +
    facet_wrap(~kingdom, nrow = 1) +
    labs(color = "Environment") +
    scale_color_environment_group() +

    # other kingdom arrows from dbRDA
    ggnewscale::new_scale_color() +
    geom_segment(
      data = tbl %>% select(dbrda_terms_tbl) %>% unnest(dbrda_terms_tbl, names_repair = "universal") %>% rename(kingdom = kingdom...1),
      arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
      mapping = aes(x = 0, y = 0, xend = dbRDA1 * 2, yend = dbRDA2 * 2, color = phylum_color),
      size = 0.5
    ) +
    scale_color_phyla() +

    # model info
    geom_text(
      data = tbl,
      mapping = aes(label = dbrda_annot_text),
      size = 5,
      x = -Inf,
      y = Inf,
      vjust = "inward",
      hjust = "inward"
    ) +
    labs(
      color = "Phylum",
      x = (tbl$proportion_explained_axis_1 * 100) %>% sprintf(fmt = "dbRDA1 (%.1f%%)"),
      y = (tbl$proportion_explained_axis_2 * 100) %>% sprintf(fmt = "dbRDA2 (%.1f%%)")
    )


  plt
}
