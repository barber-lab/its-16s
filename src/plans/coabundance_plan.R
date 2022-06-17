get_coabundance_plan <- function() {
  drake::drake_plan(
    cor_methods = c("spearman", "sparcc"),

    coabundance_min_prevalence_perc = 10,
    coabundance_n_permutations = 100,

    # filter highly correlated samples
    max_abs_samples_correlation = 1,
    samples_correlation_method = "spearman",

    spiec_easi_params = list(
      nlambda = 300,
      lambda.min.ratio = 1e-2,
      sel.criterion = "stars",
      pulsar.select = TRUE,
      pulsar.params = list(
        thresh = 0.05,
        subsample.ratio = 0.8,
        ncores = 32,
        rep.num = 20,
        seed = 1337
      )
    ),

    sparcc_params = list(
      iterations = 50,
      exclude_iterations = 20,
      bootstraps = 200
    ),
    
    jaccard_params = list(
      B = 10e3,
      min_abundance_perc = 0.1
    ),

    #
    # Analyze co abundance ----
    #

    inter_kingdom_mats = target(
      {
        sub_abundances %>%
          filter(norm_method == "raw" & !str_detect(subset_name, "kingdom")) %>%
          inner_join(prevalences %>% rename(prev_tbl = data)) %>%
          mutate(
            data = data %>% map2(prev_tbl, ~ {
              .x %<>%
                filter_samples_correlation(
                  max_abs_estimate = max_abs_samples_correlation,
                  method = samples_correlation_method
                )

              taxa <-
                .y %>%
                filter(prevalence_perc > coabundance_min_prevalence_perc) %>%
                pull(taxon)

              # only take samples with counts in both kingdoms
              inter_samples <-
                .x %>%
                filter(taxon %in% taxa) %>%
                distinct(sample_id, kingdom) %>%
                mutate(is_present = TRUE) %>%
                pivot_wider(names_from = kingdom, values_from = is_present, values_fill = list(is_present = FALSE)) %>%
                filter(Bacteria & Fungi) %>%
                pull(sample_id) %>%
                unique()

              .x %>%
                filter(
                  sample_id %in% inter_samples &
                    kingdom %in% c("Bacteria", "Fungi") &
                    taxon %in% taxa
                ) %>%
                arrange(sample_id, kingdom, taxon) %>%
                group_by(kingdom) %>%
                nest() %>%
                mutate(
                  mat = data %>% map(~ {
                    .x %>%
                      select(sample_id, taxon, abundance) %>%
                      pivot_wider(names_from = "taxon", values_from = "abundance", values_fill = list(abundance = 0)) %>%
                      as_matrix("sample_id")
                  })
                ) %>%
                pull(mat)
            })
          ) %>%
          select(-prev_tbl)
      },
      # do in parallel due to slow iterative sample correlation filter
      hpc = TRUE,
      dynamic = map(sub_abundances)
    ),

    intra_kingdom_mats = target(
      {
        sub_abundances %>%
          filter(norm_method == "raw" & str_detect(subset_name, "kingdom")) %>%
          inner_join(prevalences %>% rename(prev_tbl = data)) %>%
          mutate(
            data = data %>% map2(prev_tbl, ~ {
              .x %<>%
                filter_samples_correlation(
                  max_abs_estimate = max_abs_samples_correlation,
                  method = samples_correlation_method
                )

              taxa <-
                .y %>%
                filter(prevalence_perc > coabundance_min_prevalence_perc) %>%
                pull(taxon)

              .x %>%
                filter(taxon %in% taxa) %>%
                select(sample_id, taxon, abundance) %>%
                pivot_wider(names_from = "taxon", values_from = "abundance", values_fill = list(abundance = 0)) %>%
                as_matrix("sample_id") %>%
                list()
            })
          ) %>%
          select(-prev_tbl)
      },
      # do in parallel due to slow iterative sampel correlation filter
      hpc = TRUE,
      dynamic = map(sub_abundances)
    ),

    permuted_mats = {
      true_abundances <-
        sub_mats %>%
        filter(subset_name == "environment_group" & norm_method == "raw") %>%
        unnest(data) %>%
        mutate(data = data %>% map(~ {
          .x %>%
            as_tibble(rownames = "sample_id") %>%
            pivot_longer(-sample_id, names_to = "taxon", values_to = "abundance")
        })) %>%
        unnest(data)

      seq(coabundance_n_permutations) %>%
        map(~ {
          set.seed(.x)

          subset_value_tbl <-
            true_abundances %>%
            distinct(sample_id, subset_value) %>%
            mutate(subset_value = subset_value %>% sample())

          true_abundances %>%
            select(sample_id, taxon, abundance) %>%
            inner_join(subset_value_tbl, by = "sample_id")
        }) %>%
        enframe(name = "permutation", value = "abundance") %>%

        # re nest
        unnest(abundance) %>%
        group_by(permutation, subset_value) %>%
        nest() %>%
        rename(abundance = data) %>%
        ungroup() %>%

        # long tbl to matrix
        mutate(abundance = abundance %>% map(~ {
          .x %>%
            pivot_wider(names_from = "taxon", values_from = "abundance", values_fill = list(abundance = 0)) %>%
            as_matrix("sample_id")
          set_rownames(.$sample_id) %>%
            select(-sample_id) %>%
            as.matrix()
        }))
    },

    sub_mats = {
      bind_rows(
        inter_kingdom_mats
        #intra_kingdom_mats
      ) %>%
        filter(! subset_name %>% str_detect("habitat"))
    },

    permuted_coabundances = target(
      {
        data <- permuted_mats$abundance[[1]]

        cor_res <-
          sparcc_params %>%
          inset2("data", data) %>%
          # We don't need a p value to calculate FDR of permuted env groups
          inset2("bootstraps", 1) %>%
          inset2("threads", 1) %>%
          possibly(do.call, NA)(what = correlate_fastspar)

        if (is.na(cor_res)) {
          return(tibble())
        }

        permuted_mats %>%
          purrr::discard(is.list) %>%
          dplyr::mutate(
            cor_res = list(cor_res),
          )
      },
      dynamic = map(permuted_mats),
      hpc = TRUE
    ),

    permuted_coabundances_file = {
      permuted_coabundances %>%
        write_rds(file_out("results/permuted_coabundances.rds"), compress = "gz")
    },

    # distinct from target coabundance becacause function coabundance is under dev
    coabundances = target(
      {
        sub_mats %>%
          expand_grid(cor_method = cor_methods) %>%
          mutate(
            cor_res = data %>% map2(cor_method, possibly(~ {
              switch(.y,
                "mb" = {
                  spiec_easi_params %>%
                    inset2("data", .x) %>%
                    inset2("method", "mb") %>%
                    possibly(do.call, NA)(what = spiec.easi)
                },
                "glasso" = {
                  spiec_easi_params %>%
                    inset2("data", .x) %>%
                    inset2("method", "glasso") %>%
                    possibly(do.call, NA)(what = spiec.easi)
                },
                "sparcc" = {
                  data <- .x

                  # merge multiple matrices
                  if (length(data) > 1) {
                    data <-
                      data %>%
                      map(~ .x %>%
                        as_tibble(rownames = "sample_id") %>%
                        pivot_longer(!sample_id)) %>%
                      bind_rows() %>%
                      pivot_wider(values_fill = list(value = 0)) %>%
                      set_rownames(.$sample_id) %>%
                      select(!sample_id) %>%
                      as.matrix()
                  }

                  if (length(data) == 1 && class(data) == "list") {
                    data <- data[[1]]
                  }

                  sparcc_params %>%
                    inset2("data", data) %>%
                    possibly(do.call, NA)(what = correlate_fastspar)
                },
                "spearman" = {
                  possibly(Hmisc::rcorr, NA)(sub_mats$data[[1]] %>% reduce_inter_kingdom_mats(), type = "spearman")
                },
                "pearson" = {
                  possibly(Hmisc::rcorr, NA)(sub_mats$data[[1]] %>% reduce_inter_kingdom_mats(), type = "pearson")
                },
                error = function(e) stop("method not implemented")
              )
            }, NA))
          ) %>%
          select(-where(is.list), cor_res)
      },
      dynamic = cross(sub_mats, cor_methods),
      # use internal parallelism for memory efficiency
      hpc = FALSE
    ),

    # Interactions between taxa like in coabunadacne
    # but treating each taxa as either present or absent instead
    cooccurrences = target(
      {
        # Reduce parallelization overhead
        options(mc.cores = 10)

        sub_mats %>%
          mutate(
            data = data %>% map(reduce_inter_kingdom_mats),
            cooccurence = data %>% map(~ {
              jaccard_params %>%
                inset2("data", .x) %>%
                possibly(do.call, NA)(what = correlate_jaccard)
            })
          ) %>%
          select(-where(is.list), cooccurence)
      },
      dynamic = map(sub_mats),
      hpc = TRUE
    )
  )
}
