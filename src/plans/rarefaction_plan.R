#' Sub plan for FEAST source tracking
get_rarefaction_plan <- function() {
  drake_plan(
    taxa_deck = target(
      # only samples passed sanity checks e.g. have both fungal and bacterial abundance
      sub_abundances$data[[1]] %>%
        left_join(lineages),
      hpc = FALSE
    ),
    rarefactions_trails = seq(10),
    rarefactions_sample_groupings = c("bioproject_id", "environment_group"),
    rarefactions = target(
      {
        tibble(sample_grouping = rarefactions_sample_groupings) %>%
          mutate(rarefaction = sample_grouping %>% map(
            possibly(~ get_rarefaction(
              sample_grouping = .x,
              rarefactions_trails = rarefactions_trails,
              samples = samples,
              alphadiv_metrics = alphadiv_metrics,
              kingdoms_groups = kingdoms_groups,
              abundances = abundances,
              lineages = lineages
            ), NA)
          ))
      },
      dynamic = map(rarefactions_sample_groupings)
    ),

    rarefactions_environment_group_plt = target(
      {
        rarefactions %>%
          filter(sample_grouping == "environment_group") %>%
          pull(rarefaction) %>%
          first() %>%
          filter(alphadiv_metric == "Observed") %>%
          # force line start at origin
          complete(environment_group, kingdom, alphadiv_metric, fill = list(alpha_diversity = 0)) %>%
          ggplot(aes(n_samples, alpha_diversity, color = environment_group, fill = environment_group)) +
          stat_smooth(method = "loess", alpha = 0.2, span = 1) +
          facet_wrap(~kingdom, scales = "free_y") +
          scale_fill_environment_group() +
          scale_color_environment_group() +
          scale_x_continuous(limits = c(1, 110), expand = c(0, 0)) +
          scale_y_continuous(expand = c(0, 0)) +
          guides(fill = FALSE, color = FALSE) +
          labs(x = "Project subset size (Samples)", y = str_glue("Observed genera"), color = "Environment")
      },
      hpc = FALSE
    ),

    rarefactions_bioproject_id_plt = target(
      {
        rarefactions %>%
          filter(sample_grouping == "bioproject_id") %>%
          pull(rarefaction) %>%
          first() %>%
          filter(alphadiv_metric == "Observed") %>%
          # annotate environment_group
          left_join(
            samples %>%
              count(bioproject_id, environment_group) %>%
              group_by(bioproject_id) %>%
              arrange(-n) %>%
              slice(1)
          ) %>%
          # force line start at origin
          complete(bioproject_id, kingdom, alphadiv_metric, fill = list(alpha_diversity = 0)) %>%
          ggplot(aes(n_samples, alpha_diversity, group = bioproject_id, color = environment_group, fill = environment_group)) +
          stat_smooth(method = "loess", alpha = 0.2, span = 1) +
          facet_wrap(~kingdom, scales = "free_y") +
          scale_fill_environment_group() +
          scale_color_environment_group() +
          scale_x_continuous(limits = c(1, 110), expand = c(0, 0)) +
          scale_y_continuous(expand = c(0, 0)) +
          guides(fill = FALSE, color = FALSE) +
          labs(x = "Project subset size (Samples)", y = str_glue("Observed genera"), color = "Environment")
      },
      hpc = FALSE
    )
  )
}
