get_dgca_plan <- function() {
  drake_plan(
    dgca_norm_methods = c("raw", "tss", "clr"),

    dgca_inputMat = target(
      {
        tibble(norm_method = dgca_norm_methods) %>%
          mutate(
            inputMat = norm_method %>% map(~ {
              res <-
                sub_mats %>%
                filter(subset_value == "all") %>%
                pull(data) %>%
                first() %>%
                map(as.data.frame) %>%
                reduce(bind_cols) %>%
                as_tibble(rownames = "sample_id") %>%
                pivot_longer(-sample_id, names_to = "taxon", values_to = "abundance") %>%
                left_join(lineages)

              # normalization
              switch(.x,
                clr = {
                  res <-
                    res %>%
                    group_by(sample_id, kingdom) %>%
                    mutate(abundance = clr(abundance))
                },
                tss = {
                  res <-
                    res %>%
                    group_by(sample_id, kingdom) %>%
                    mutate(abundance = abundance / sum(abundance))
                },
                raw = {
                  # do no edit
                }
              )

              res %>%
                ungroup() %>%
                select(sample_id, taxon, abundance) %>%
                pivot_wider(names_from = taxon, values_from = "abundance") %>%

                # add meta data
                semi_join(samples) %>%
                column_to_rownames("sample_id") %>%
                as.matrix() %>%
                t()
            })
          )
      },
      dynamic = map(dgca_norm_methods)
    ),

    dgca_design = {
      dgca_inputMat$inputMat[[1]] %>%
        colnames() %>%
        tibble(sample_id = .) %>%
        left_join(samples) %>%
        model.matrix(~ environment_group - 1, .)
    },

    dgca_comparisons = {
      environment_groups %>%
        combn(2) %>%
        t() %>%
        as_tibble() %>%
        rename(from = V1, to = V2)
    },

    dgca_results = target(
      {
        dgca_comparisons %>%
          expand_grid(dgca_inputMat) %>%
          mutate(
            data = list(from, to, inputMat) %>% pmap(~ {
              DGCA::ddcorAll(
                inputMat = ..3,
                design = dgca_design,
                compare = c(..1, ..2) %>% paste("environment_group", ., sep = ""),
                corrType = "pearson",
                nPerms = 1000
              )
            })
          ) %>%
          select(-where(is.list), data)
      },
      dynamic = cross(dgca_comparisons, dgca_inputMat),
      hpc = TRUE
    ),

    dgca_differential_correlations = target(
      {
        dgca_results %>%
          mutate(
            differential_correlation = list(from, to, data) %>% pmap(function(from, to, data) {
              data %>%
                as_tibble() %>%
                mutate(from_env = from, to_env = to) %>%
                select(from_env, to_env, from_taxon = Gene1, to_taxon = Gene2, p.value = empPVals, matches("environment_group[A-z]+_cor")) %>%
                group_by(from_env, to_env) %>%
                mutate(
                  q.value = p.value %>% p.adjust(method = "fdr"),
                  is_diff = q.value < 0.05
                ) %>%
                filter(!is.na(from_taxon)) %>%
                unite(comp, c(from_env, to_env)) %>%
                # differential <=> not significant in all env comparisons
                group_by(comp, from_taxon, to_taxon) %>%
                mutate(any_differential = any(is_diff)) %>%
                pivot_longer(matches("environment_group[A-z]+_cor"), names_to = "environment_group", values_to = "estimate") %>%

                # e.g. there is no soil estimate in comparing host vs aquatic
                filter(!is.na(estimate)) %>%
                filter(any_differential) %>%
                mutate(
                  abs_estimate = abs(estimate),
                  sign_estimate = sign(estimate),
                  environment_group = environment_group %>% str_remove_all("^environment_group|_cor$")
                ) %>%
                ungroup() %>%
                select(-comp) %>%
                arrange(from_taxon, to_taxon)
            })
          ) %>%
          select(-where(is.list), differential_correlation)
      },
      dynamic = map(dgca_results)
    )
  )
}
