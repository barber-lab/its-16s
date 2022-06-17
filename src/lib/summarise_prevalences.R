summarise_prevalences <- function(prevalences, grouping = "environment_group") {
  prevalences %>%
    summarise_subsets(grouping) %>%
    mutate(taxa = data %>% map(~ {
      .x %>%
        dplyr::select(group, taxon) %>%
        group_by(group) %>%
        nest() %>%
        unnest_wider(data) %>%
        deframe() %>%
        list_modify(all = NULL)
    })) %>%
    filter(taxa %>% map_lgl(~ length(.x) == 3)) %>%
    mutate(venn_plt = taxa %>% map2(grouping, ~ {
      plt <-
        .x %>%
        ggvenn::ggvenn(
          text_size = 3,
          set_name_size = 0,
          show_percentage = FALSE
        )

      switch(grouping,
        "environment_group" = {
          plt +
            scale_fill_manual(
              values = c(
                "A" = environment_group_colors[["host"]],
                "B" = environment_group_colors[["aquatic"]],
                "C" = environment_group_colors[["soil"]]
              )
            )
        },
        "dysbalance" = {
          plt +
            scale_fill_manual(
              values = c(
                "A" = dysbalance_colors[["balanced"]],
                "B" = dysbalance_colors[["dysbalanced"]]
              )
            )
        },
        plt
      )
    })) %>%
    arrange_columns()
}
