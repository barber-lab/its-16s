plot_alphadiv <- function(alphadiv, alphadiv_metrics, samples) {
  alphadiv_tbl <-
    alphadiv %>%
    filter(subset_name == "environment_group,kingdom") %>%
    unnest(data) %>%
    separate(subset_value, into = c("environment_group", "kingdom")) %>%
    pivot_longer(alphadiv_metrics, names_to = "alphadiv_metric", values_to = "alphadiv_value")

  # alphadiv_test_tbl <-
  #   alphadiv_tbl %>%
  #   left_join(samples) %>%
  #   group_by(alphadiv_metric, kingdom) %>%
  #   nest() %>%
  #   transmute(
  #     test = data %>% map(~ lmer(alphadiv_value ~ environment_group + (1|bioproject_id), data = .x))
  #   ) %>%
  #   select(-data) %>%
  #   unnest_wider(test) %>%
  #   mutate(sig_label = significance_label(p.value))
  
  
  # lvplot does not work with kernlab
  # @see https://github.com/hadley/lvplot/issues/8
  try(detach("package:kernlab", unload=TRUE))

  alphadiv_plts <-
    alphadiv_tbl %>%
    filter(alphadiv_metric %in% c("Shannon", "Chao1")) %>%
    # plot labels
    rename(Environment = environment_group, Kingdom = kingdom, `Alpha diversity` = alphadiv_metric) %>%
    group_by(Kingdom, `Alpha diversity`) %>%
    nest() %>%
    
    # draw individual plots because we want to have a free scale for every tile and not just row/col
    transmute(
      plt = list(data, Kingdom, `Alpha diversity`) %>% pmap(~ {
        ..1 %>%
          mutate(facet_x = ..3, facet_y = str_glue("{..2}\n{..3}")) %>%
          ggplot(aes(Environment, alphadiv_value, fill = Environment)) +
          lvplot::geom_lv(color = "black") +
          scale_fill_environment_group() +
          stat_compare_means_environment_group(method = "wilcox", step.increase = 0.25) +
          scale_y_continuous(
            trans = log10_trans(),
            breaks = trans_breaks("log10", function(x) 10^x),
            labels = trans_format("log10", math_format(10^.x)),
            expand = c(0.2, 0)
          ) +
          annotation_logticks(sides = "l") +
          facet_wrap(~ facet_y) +
          guides(fill = FALSE) +
          theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank()
          ) +
          labs(y = "")
      })
    ) %>%
    arrange(`Alpha diversity`, Kingdom) %>%
    pull(plt)
  
  alphadiv_plts[[1]] <- alphadiv_plts[[1]] + labs(tag = "A")
  
  alphadiv_plts %>%
    wrap_plots(nrow = 1, byrow = FALSE, guides = "collect")
}
