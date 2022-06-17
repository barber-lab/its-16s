#' Get plan to run linear mixed models
get_lmm_plan <- function() {
  drake_plan(
    lmm_abundance = {
      sub_abundances %>%
        filter(subset_name == "kingdom" & norm_method == "tss") %>%
        unnest(data) %>%
        group_by(kingdom, taxon) %>%
        nest() %>%
        transmute(
          m = data %>% map(possibly(~ {
            lme4::lmer(abundance ~ environment_group * (1|bioproject_id), data = .x) %>%
              car::Anova() %>%
              tidy()
          }, NA))
        ) %>%
        filter(! is.na(m)) %>%
        unnest(m) %>%
        mutate(q.value = p.value %>% p.adjust(method = "fdr")) %>%
        arrange(kingdom, taxon)      
    },
    
    lmm_alphadiv = {
      alphadiv %>%
        filter(subset_name == "kingdom") %>%
        rename(kingdom = subset_value) %>%
        unnest(data) %>%
        inner_join(samples) %>%
        
        # mutate collection data
        mutate(
          collection_datetime = collection_datetime %>% as.numeric()
        ) %>%
        
        pivot_longer(cols = alphadiv_metrics, names_to = "alphadiv_metric", values_to = "alphadiv") %>%
        group_by(kingdom, alphadiv_metric) %>%
        nest() %>%
        expand_grid(formula = c(
          "alphadiv ~ environment_group * (1|bioproject_id)",
          "alphadiv ~ habitat * (1|bioproject_id)"
        )) %>%
        group_by(kingdom, alphadiv_metric, formula) %>%
        transmute(
          m = data %>% map2(formula, possibly(~ {
            .y %>%
              as.formula() %>%
              lme4::lmer(data = .x) %>%
              car::Anova() %>%
              tidy()
          }, NA))
        ) %>%
        filter(! is.na(m)) %>%
        unnest(m) %>%
        mutate(q.value = p.value %>% p.adjust(method = "fdr")) %>%
        ungroup() %>%
        arrange(formula, kingdom, q.value) 
    }
  )
}