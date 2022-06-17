plot_prevalences_abundances <- function(prevalences, abundances, samples, lineages) {
  top_prevalent_taxa <-
    prevalences %>%
    filter(subset_name == "environment_group,kingdom") %>%
    unnest(data) %>%
    group_by(subset_value) %>%
    arrange(-prevalence_perc) %>%
    slice(1:10) %>%
    pull(taxon) %>%
    unique()

  bind_rows(
    prevalences %>%
      filter(subset_name == "environment_group,kingdom") %>%
      unnest(data) %>%
      filter(taxon %in% top_prevalent_taxa) %>%
      separate(col = subset_value, into = c("environment_group", "kingdom"), sep = ",") %>%
      select(environment_group, kingdom, taxon, value = prevalence_perc) %>%
      complete(taxon, environment_group, fill = list(value = 0)) %>%
      mutate(type = "prevalence\n(%)"),

    abundances %>%
      filter(norm_method == "raw") %>%
      pull(data) %>%
      first() %>%
      group_by(sample_id) %>%
      mutate(abundance = abundance %>% rank() %>% scale_min_max() * 100) %>%
      filter(taxon %in% top_prevalent_taxa) %>%
      left_join(samples %>% select(sample_id, environment_group)) %>%
      left_join(lineages %>% select(taxon, kingdom)) %>%
      filter(!is.na(environment_group)) %>%
      group_by(kingdom, taxon, environment_group) %>%
      summarise(abundance = mean(abundance)) %>%
      ungroup() %>%
      select(kingdom, taxon, environment_group, value = abundance) %>%
      complete(taxon, environment_group, fill = list(value = 0)) %>%
      mutate(type = "abundance\n(percentile)")
  ) %>%
    # Re-annotate kingdom (NA in rows added with value of 0)
    select(-kingdom) %>%
    left_join(lineages %>% select(taxon, kingdom)) %>%
    mutate(
      taxon = taxon %>% factor(levels = top_prevalent_taxa)
    ) %>%
    nest(-c(kingdom, type)) %>%
    arrange(kingdom, type) %>%
    mutate(
      type = type %>% factor(levels = c("prevalence\n(%)", "abundance\n(percentile)")),
      show_y_labels = c(FALSE, TRUE, FALSE, TRUE),
      plt = list(data, type, show_y_labels, kingdom) %>% pmap(~ {
        res <-
          ..1 %>%
          mutate(facet = ..4) %>%
          ggplot(aes(x = taxon)) +
          geom_linerange(
            data = {
              ..1 %>%
                group_by(taxon) %>%
                summarise(ymin = min(value), ymax = max(value))
            },
            mapping = aes(ymin = ymin, ymax = ymax)
          ) +
          geom_point(
            mapping = aes(y = value, color = environment_group),
            size = 2,
            position = position_dodge(width = 0.5)
          ) +
          scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, length.out = 3)) +
          scale_color_environment_group() +
          coord_flip() +
          facet_wrap(~facet) +
          theme(axis.text.y = element_text(face = "italic")) +
          guides(color = guide_legend(override.aes = list(size = 6))) +
          labs(y = ..2, x = "", color = "Environment")

        if (!..3) {
          res <- res + theme(
            axis.text.y = element_blank(),
            axis.title.y = element_blank()
          )
        }

        res
      })
    ) %>%
    arrange(kingdom, type) %>%
    pull(plt) %>%
    wrap_plots(nrow = 1, guides = "collect")
}
