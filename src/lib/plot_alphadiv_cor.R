plot_alphadiv_cor <- function(alphadiv, alphadiv_metrics, samples) {
  # visualize simpson's paradox
  alphadiv %>%
    filter(subset_name == "environment_group,kingdom") %>%
    select(-subset_name) %>%
    unnest(data) %>%
    select(c("subset_value", "sample_id"), alphadiv_metrics) %>%
    separate(subset_value, into = c("environment_group", "kingdom")) %>%
    distinct_all() %>%
    pivot_longer(alphadiv_metrics, names_to = "alphadiv_metric", values_to = "alphadiv") %>%
    pivot_wider(names_from = kingdom, values_from = alphadiv) %>%
    left_join(samples %>% select(sample_id, bioproject_id)) %>%
    filter(alphadiv_metric %in% c("Shannon", "Chao1")) %>%
    
    group_by(alphadiv_metric) %>%
    nest() %>%
    transmute(
      cor = data %>% map(~ {
        .x %>%
          select(-sample_id) %>%
          group_by(environment_group) %>%
          nest() %>%
          transmute(
            cor = data %>% map(~ correlation::correlation(.x, multilevel = TRUE))
          ) %>%
          unnest(cor)
      }),
      
      plt = list(data,  alphadiv_metric, cor) %>% pmap(~ {
        x <- ..1
        if(..2 == "Chao1") {
          x %<>%
            mutate(Bacteria = Bacteria %>% log10()) %>%
            mutate(Fungi = Fungi %>% log10())
        }
        
        multienv_bioprojects <-
          x %>%
          count(bioproject_id, environment_group) %>%
          count(bioproject_id) %>%
          filter(n > 1) %>%
          pull(bioproject_id)
        
        x %<>% filter(! bioproject_id %in% multienv_bioprojects)
        
        cors <-
          ..3 %>%
          mutate(
            r = r %>% round(2) %>% sprintf(fmt = "%.2f"),
            p = p %>% significance_label()
          ) %>%
          select(Environment = environment_group, r, p) %>%
          column_to_rownames("Environment")
        
        x %>%
          ggplot(aes(Bacteria, Fungi)) +
          stat_density2d_filled(alpha = 0.7) +
          geom_smooth(aes(group = bioproject_id, color = environment_group), method = "lm", se = FALSE, alpha = 0.2) +
          scale_color_environment_group() +
          scale_fill_brewer(palette = "Greys", direction = 1, na.value = "black") +
          scale_x_continuous(expand = c(0, 0)) +
          scale_y_continuous(expand = c(0, 0)) +
          guides(fill = FALSE, color = FALSE) +
          # annotation_custom(
          #   gridExtra::tableGrob(cors, theme = gridExtra::ttheme_minimal(base_size = 8)),
          #   xmax = max(x$Bacteria),
          #   ymin = max(x$Fungi) * 0.5
          # ) +
          labs(
            tag = "B",
            x = switch(..2, "Chao1" = TeX("Bacteria (log_{10}(Chao1))"), "Shannon" = "Bacteria (Shannon)"),
            y = switch(..2, "Chao1" = TeX("Fungi (log_{10}(Chao1))"), "Shannon" = "Fungi (Shannon)")
          )        
      })
    ) %>%
    pull(plt) %>%
    wrap_plots(nrow = 1) +
    labs(tag = "") # remove last tag
}