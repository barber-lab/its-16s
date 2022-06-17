#' Sub plan for selbal supervised balance
get_selbal_plan <- function() {
  drake_plan(
    selbal_min_prevalence_perc = 10,
    
    selbal_subsets = {
      tribble(
        ~group1, ~group2,
        "aquatic", "host",
        "aquatic", "soil",
        "host", "soil"
      ) %>%
        expand_grid(
          sub_abundances %>%
            dplyr::filter(subset_name == "kingdom" & norm_method == "raw") %>%
            dplyr::select(kingdom = subset_value, abundance = data)
        )
    },

    selbal_balances = target(
      {
        abundance <- selbal_subsets$abundance[[1]]
        group1 <- selbal_subsets$group1[[1]]
        group2 <- selbal_subsets$group2[[1]]

        res <- tryCatch(
          {
            samples_counts <-
              abundance %>%
              distinct(environment_group, sample_id) %>%
              dplyr::count(environment_group, name = "n_total_samples")

            prevalent_taxa <-
              abundance %>%
              group_by(environment_group, taxon) %>%
              dplyr::count(name = "n_samples") %>%
              inner_join(samples_counts) %>%
              transmute(prevalence_perc = n_samples / n_total_samples * 100) %>%
              pivot_wider(
                names_from = environment_group,
                values_from = prevalence_perc,
                values_fill = list(prevalence_perc = 0)
              ) %>%
              # keep taxa prevalent in both groups
              filter_at(group1, ~ .x >= selbal_min_prevalence_perc) %>%
              filter_at(group2, ~ .x >= selbal_min_prevalence_perc) %>%
              pull(taxon)

            data <-
              abundance %>%
              filter(taxon %in% prevalent_taxa) %>%
              filter(environment_group %in% c(group1, group2))

            x <-
              data %>%
              select(sample_id, taxon, abundance) %>%
              pivot_wider(names_from = taxon, values_from = abundance, values_fill = list(abundance = 0)) %>%
              column_to_rownames("sample_id")

            y <-
              data %>%
              distinct(sample_id, environment_group) %>%
              pull(environment_group) %>%
              # need fixed orders for matching
              factor(levels = c(group1, group2))
            
            covar <-
              x %>%
              rownames() %>%
              tibble(sample_id = .) %>%
              left_join(samples) %>%
              transmute(
                sample_id,
                #habitat,
                lon,
                lat,
                year = year(collection_datetime),
                month = month(collection_datetime)
              ) %>%
              column_to_rownames("sample_id") %>%
              mice::mice() %>%
              mice::complete() %>%
              as.data.frame()
              # # selbal needs numeric covariates
              # dummy_cols("habitat", remove_selected_columns = TRUE)
              
            res_cv <- selbal::selbal.cv(x = x, y = y, covar = covar)
            res_global <- selbal::selbal(
              x = x,
              y = y,
              covar = covar,
              tab = TRUE,
              draw = FALSE,
              maxV = res_cv$opt.nvar
            )
            names(res_global[[1]]) <- rownames(x)

            list(cv = res_cv, global = res_global)
          },
          error = function(e) NA
        )

        selbal_subsets %>%
          select(-where(is.list)) %>%
          mutate(balance = list(res))
      },
      dynamic = map(selbal_subsets),
      hpc = TRUE
    ),
    
    selbal_plt = {
      selbal_taxa_plt <-
        selbal_balances %>%
        filter(! is.na(balance)) %>%
        mutate(
          balance_tbl = list(balance, group1, group2) %>% pmap(~ {
            ..1 %>%
              pluck("cv", "cv.tab") %>%
              as_tibble(rownames = "taxon") %>%
              select(taxon, resamples_present_perc = `%`, environment_group = Global) %>%
              mutate(environment_group = environment_group %>% recode(DEN = ..2, NUM = ..3, .default = "")) %>%
              filter(taxon != "FREQ") %>%
              type_convert()
          })
        ) %>%
        select(-balance) %>%
        unnest(balance_tbl) %>%
        left_join(lineages) %>%
        mutate(comp = paste0(group1, " vs. ", group2)) %>%
        filter(resamples_present_perc >= 90) %>%
        
        ggplot(aes(comp, taxon, fill = environment_group)) +
        geom_tile() +
        facet_wrap(~kingdom, ncol = 1, scales = "free_y", strip.position = "right") +
        scale_fill_environment_group() +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        theme(axis.text.y = element_text(face = "italic")) +
        labs(x = "Comparison", y = "Selected taxon", fill = "Environment")

      selbal_balances_unnested_tbl <-
        selbal_balances %>%
        filter(! is.na(balance)) %>%
        transmute(
          group1, group2, kingdom,
          balance = balance %>% map(~ {
            .x$cv$glm$data %>%
              as_tibble(rownames = "sample_id") %>%
              dplyr::select(sample_id, balance = V1) %>%
              inner_join(samples %>% select(sample_id, group = environment_group))
          })
        ) %>%
        unnest(balance) %>%
        mutate(comp = paste0(group1, " vs. ", group2))
      
      selbal_balances_tests_tbl <-
        selbal_balances_unnested_tbl %>%
        group_by(kingdom, comp) %>%
        do(wilcox.test(balance ~ group, data = .) %>% tidy()) %>%
        transmute(
          label = case_when(
            p.value < 1e-3 ~ "★★★",
            TRUE ~ "?"
          )
        )

      selbal_balances_plt <-
        selbal_balances_unnested_tbl %>%
        ggplot(aes(comp, balance)) +
        geom_boxplot(aes(fill = group)) +
        geom_text(
          data = selbal_balances_tests_tbl,
          mapping = aes(label = label),
          y = Inf,
          size = 5
        ) +
        facet_wrap(~kingdom, ncol = 1, strip.position = "right") +
        scale_fill_environment_group() +
        coord_cartesian(clip = "off") +
        labs(x = "Comparison", y = "Supervised balance")
      
      list(
        selbal_taxa_plt + theme(axis.title.x = element_blank()) + labs(tag = "A"),
        selbal_balances_plt + guides(fill = FALSE) + labs(tag = "B")
      ) %>%
        wrap_plots(ncol = 1) +
        plot_layout(guides = "collect")
    }
  )
}
