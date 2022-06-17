#' Sub plan for Hierarchical Modelling of Species Communities (HMSC)
get_hmsc_plan <- function() {
  drake_plan(
    #' Iterations spend on training
    #' Will be multiplied by hmsc_thin and added by hmsc_warmup_iter
    #' to get total number of iterations
    #' @see https://www.nature.com/articles/s41396-020-0732-1
    #'
    #' Thinning must be set before sampling to have small objects working good qith qs and drake
    hmsc_sampling_params = {
      tribble(
        ~hmsc_iter, ~hmsc_thin, ~hmsc_chains,
        10e3, 10, 2
      ) %>%
        mutate(
          hmsc_iter = hmsc_iter / hmsc_thin,
          hmsc_warmup_iter = as.integer(hmsc_iter / 3) * hmsc_thin
        )
    },

    hmsc_support_level = 0.05,
    hmsc_post_estimate_parameters = c("Beta", "OmegaCor"),

    # Min TSS abundance to count as occurring
    hmsc_min_prevalence_abundance = 0,
    hmsc_min_prevalence_perc = 20,
    hmsc_sample_attributes = c("bioproject_id", "habitat", "environment_group", "collection_datetime"),

    hmsc_samples_formulas = tribble(
      ~samples_formula, ~model_type,
      "~ seqdepth", "null", # null model e.g. for raw cooccurrence
      "~ seqdepth + environment_group - 1", "residual_environment_group", # full model e.g. for restricted cooccurrence
      "~ seqdepth + habitat - 1", "residual_habitat", # full model e.g. for restricted cooccurrence
      "~ seqdepth + distance_to_medoid_Bacteria + distance_to_medoid_Fungi - 1", "residual_dysbalance" # full model e.g. for restricted cooccurrence
    ),

    hmsc_traits_formulas = bind_rows(
      subsets %>%
        filter(subset_name %>% str_detect("kingdom") & subset_value %>% str_detect("Bacteria")) %>%
        select(-subset_query) %>%
        mutate(traits_formula = "~ phylum + faprotax_entry - 1"),

      subsets %>%
        filter(subset_name %>% str_detect("kingdom") & subset_value %>% str_detect("Fungi")) %>%
        select(-subset_query) %>%
        mutate(traits_formula = "~ phylum + trophicMode - 1"),

      subsets %>%
        filter(!subset_name %>% str_detect("kingdom")) %>%
        select(-subset_query) %>%
        mutate(traits_formula = "~ phylum + trophicMode + faprotax_entry - 1"),

      # Models without traits
      subsets %>%
        select(-subset_query) %>%
        mutate(traits_formula = "~1")
    ),

    hmsc_subsets = {
      subsets %>% filter(subset_name %in% c("all", "kingdom", "environment_group", "environment_group,kingdom"))
    },

    hmsc_taxa = target(
      {
        prevalent_taxa <-
          prevalences %>%
          inner_join(hmsc_subsets) %>%
          unnest(data) %>%
          filter(prevalence_perc > hmsc_min_prevalence_perc) %>%
          pull(taxon)

        hmsc_taxa <-
          augmented_abundances %>%
          filter(norm_method == "raw") %>%
          pull(data) %>%
          first() %>%
          pull(taxon) %>%
          unique() %>%
          tibble(taxon = .) %>%
          filter(taxon %in% prevalent_taxa) %>%
          mutate(clean_taxon = taxon %>% janitor::make_clean_names())

        hmsc_subsets %>%
          select(-subset_query) %>%
          mutate(hmsc_taxa = list(hmsc_taxa))
      },
      dynamic = map(hmsc_subsets)
    ),

    hmsc_all_taxa = {
      hmsc_taxa %>%
        pull(hmsc_taxa) %>%
        reduce(full_join)
    },

    hmsc_tbl = target(
      {
        filter_query <-
          hmsc_taxa %>%
          inner_join(subsets) %>%
          pull(subset_query) %>%
          first()

        res <-
          augmented_abundances %>%
          filter(norm_method == "tss") %>%
          pull(data) %>%
          first() %>%
          inner_join(samples) %>%
          inner_join(hmsc_taxa$hmsc_taxa[[1]]) %>%
          filter_at(hmsc_sample_attributes, ~ !is.na(.x)) %>%
          filter(filter_query %>% parse(text = .) %>% eval())

        hmsc_taxa %>%
          select(-where(is.list)) %>%
          mutate(hmsc_tbl = list(res))
      },
      dynamic = map(hmsc_taxa)
    ),

    hmsc_trees = target(
      {
        hmsc_tbl %>%
          mutate(
            tree = hmsc_tbl %>% map(possibly(~ {
              tree <-
                .x %>%
                distinct(clean_taxon, family, order, class, phylum, kingdom) %>%
                transmute(
                  pathString = paste(kingdom, phylum, clean_taxon, sep = "|")
                ) %>%
                data.tree::as.Node(pathDelimiter = "|") %>%
                ape::as.phylo()

              tree$edge.length <- tree$edge.length %>%
                length() %>%
                rep(1, .)

              tree
            }, NA))
          ) %>%
          select(-where(is.list), tree)
      },
      dynamic = map(hmsc_tbl)
    ),

    #' Species responses
    hmsc_Y_abundance = target(
      {
        res <-
          hmsc_tbl$hmsc_tbl[[1]] %>%
          select(sample_id, clean_taxon, abundance) %>%
          # Abundance normalization
          # pseudo counts
          mutate(abundance = round(abundance * 1e6)) %>%

          # Hurdle approach: set missing abundance to NA and log(counts) otherwise
          mutate(
            abundance = abundance %>% map_dbl(~ ifelse(.x >= hmsc_min_prevalence_abundance, log(.x), NA)),
          ) %>%

          # scale per taxon
          group_by(clean_taxon) %>%
          mutate(
            abundance = abundance %>% scale()
          ) %>%
          pivot_wider(names_from = clean_taxon, values_from = abundance, values_fill = list(abundance = NA)) %>%
          column_to_rownames("sample_id") %>%
          as.matrix()

        hmsc_tbl %>%
          select(-where(is.list)) %>%
          mutate(abundance = list(res))
      },
      dynamic = map(hmsc_tbl)
    ),

    #' Species responses
    hmsc_Y_abundance_non_hurdle = target(
      {
        res <-
          hmsc_tbl$hmsc_tbl[[1]] %>%
          select(sample_id, clean_taxon, abundance) %>%
          # Abundance normalization
          # pseudo counts
          mutate(abundance = round(abundance * 1e6)) %>%
          pivot_wider(names_from = clean_taxon, values_from = abundance, values_fill = list(abundance = 0)) %>%
          column_to_rownames("sample_id") %>%
          as.matrix()

        hmsc_tbl %>%
          select(-where(is.list)) %>%
          mutate(abundance = list(res))
      },
      dynamic = map(hmsc_tbl)
    ),

    hmsc_Y_occurrence = target(
      {
        res <-
          hmsc_tbl$hmsc_tbl[[1]] %>%
          select(sample_id, clean_taxon, abundance) %>%
          # occurrence
          mutate(abundance = abundance %>% map_dbl(~ ifelse(.x >= hmsc_min_prevalence_abundance, 1, 0))) %>%
          # Hurdle approach: set missing abundance to 0
          pivot_wider(names_from = clean_taxon, values_from = abundance, values_fill = list(abundance = 0)) %>%
          column_to_rownames("sample_id") %>%
          as.matrix()

        hmsc_tbl %>%
          select(-where(is.list)) %>%
          mutate(abundance = list(res))
      },
      dynamic = map(hmsc_tbl)
    ),

    hmsc_responses = {
      tribble(
        ~response_name, ~response_distribution, ~response_data,
        "occurrence", "probit", hmsc_Y_occurrence,
        "abundance", "normal", hmsc_Y_abundance
      ) %>%
        unnest(response_data)
    },

    #' Sample meta data
    hmsc_X = target(
      {
        res <-
          hmsc_tbl$hmsc_tbl[[1]] %>%
          group_by(sample_id) %>%
          slice(1) %>%
          ungroup() %>%
          select("sample_id" %>% c(hmsc_sample_attributes)) %>%
          mutate_if(is.character, as.factor) %>%
          mutate(
            # one hot encoding to have slopes for all factor levels
            aquatic = (environment_group == "aquatic") %>% as.integer(),
            host = (environment_group == "host") %>% as.integer(),
            soil = (environment_group == "soil") %>% as.integer(),
            collection_datetime = collection_datetime %>% as.numeric() + row_number()
          ) %>%
          # add dysbalances
          inner_join(
            dysbalances %>%
              select(sample_id, kingdom, distance_to_medoid, dysbalance) %>%
              pivot_wider(names_from = kingdom, values_from = c(distance_to_medoid, dysbalance))
          ) %>%
          # add sequencing depth
          inner_join(
            abundances %>%
              deframe() %>%
              pluck("raw") %>%
              group_by(sample_id) %>%
              summarise(seqdepth = sum(abundance))
          ) %>%
          column_to_rownames("sample_id")

        hmsc_tbl %>%
          select(-where(is.list)) %>%
          mutate(hmsc_X = list(res))
      },
      dynamic = map(hmsc_tbl)
    ),

    #' Taxa meta data
    hmsc_traits = target(
      {
        res <-
          hmsc_tbl$hmsc_tbl[[1]] %>%
          group_by(clean_taxon) %>%
          slice(1) %>%
          ungroup() %>%
          select(clean_taxon, taxon) %>%
          left_join(lineages) %>%
          left_join(traits) %>%
          group_by(taxon) %>%
          slice(1) %>%
          ungroup() %>%
          mutate_all(~ ifelse(is.na(.x), "NA", .x)) %>%
          select(-taxon) %>%
          column_to_rownames("clean_taxon")

        hmsc_tbl %>%
          select(-where(is.list)) %>%
          mutate(hmsc_traits = list(res))
      },
      dynamic = map(hmsc_tbl)
    ),

    hmsc_random_levels = target(
      {
        res <- list(
          bioproject_id = HmscRandomLevel(units = hmsc_X$hmsc_X[[1]]$bioproject_id)
        )

        hmsc_X %>%
          select(-where(is.list)) %>%
          mutate(random_levels = list(res))
      },
      dynamic = map(hmsc_X)
    ),

    hmsc_models_data = {
      hmsc_responses %>%
        left_join(hmsc_traits_formulas) %>%
        full_join(hmsc_samples_formulas, by = character()) %>%
        full_join(hmsc_sampling_params, by = character()) %>%
        inner_join(hmsc_X) %>%
        inner_join(hmsc_random_levels) %>%
        inner_join(hmsc_traits) %>%
        inner_join(hmsc_trees) %>%
        inner_join(hmsc_subsets) %>%

        # Controlling for env and habitat in habitat subsets does not make sense
        filter(!(subset_name %>% str_detect("habitat") & model_type != "null")) %>%
        filter(!(subset_name %>% str_detect("environment_group") & model_type == "residual_environment_group"))
    },

    hmsc_models = target(
      {
        # Resitrict Linear Algebra solver threads to minimize overhead
        #' @see https://github.com/hmsc-r/HMSC/issues/23
        RhpcBLASctl::blas_set_num_threads(4)

        paste0("log/", id_chr(), ".log") %>% sink()

        res <- tryCatch(
          {
            Hmsc::Hmsc(
              Y = hmsc_models_data$abundance[[1]],
              XData = hmsc_models_data$hmsc_X[[1]],
              XFormula = hmsc_models_data$samples_formula[[1]] %>% as.formula(),
              studyDesign = hmsc_models_data$hmsc_X[[1]],
              ranLevels = hmsc_models_data$random_levels[[1]],
              TrData = hmsc_models_data$hmsc_traits[[1]],
              TrFormula = hmsc_models_data$traits_formula[[1]] %>% as.formula(),
              phyloTree = hmsc_models_data$tree[[1]],
              distr = hmsc_models_data$response_distribution[[1]]
            ) %>%
              sampleMcmc(
                samples = hmsc_models_data$hmsc_iter[[1]], nChains = hmsc_models_data$hmsc_chains[[1]],
                thin = hmsc_models_data$hmsc_thin[[1]], transient = hmsc_models_data$hmsc_warmup_iter[[1]]
              )
          },
          error = function(e) NA
        )

        hmsc_models_data %>%
          select(-where(is.list)) %>%
          mutate(model = list(res))
      },
      dynamic = map(hmsc_models_data),
      hpc = TRUE
    ),

    hmsc_post_estimates = target(
      {
        hmsc_models %>%
          expand_grid(parameter = hmsc_post_estimate_parameters) %>%
          mutate(
            post_estimate = model %>% map2(parameter, possibly(function(model, parameter) {
              model %>%
                # this can take many minutes
                getPostEstimate(parameter) %>%
                # for mean, support, and supportNeg
                map(function(m) {
                  switch(parameter,
                    "OmegaCor" = {
                      rownames(m) <- model$spNames
                      colnames(m) <- model$spNames
                      m %>%
                        as_tibble(rownames = "from_clean_taxon") %>%
                        pivot_longer(-from_clean_taxon, names_to = "to_clean_taxon")
                    },
                    "Beta" = {
                      rownames(m) <- model$covNames

                      m %>%
                        as_tibble(rownames = "covariate") %>%
                        pivot_longer(-covariate, names_to = "clean_taxon")
                    },
                    # generic result
                    {
                      m %>%
                        as_tibble(rownames = "from") %>%
                        pivot_longer(-from, names_to = "to")
                    }
                  )
                }) %>%
                enframe() %>%
                unnest(value) %>%
                pivot_wider(names_from = name, values_from = value)
            }, NA))
          ) %>%
          select(-where(is.list), post_estimate)
      },
      dynamic = cross(hmsc_models, hmsc_post_estimate_parameters),
      hpc = TRUE
    ),

    hmsc_evaluations = target(
      {
        # Resitrict Linear Algebra solver threads to minimize overhead
        #' @see https://github.com/hmsc-r/HMSC/issues/23
        RhpcBLASctl::blas_set_num_threads(1)

        hmsc_models %>%
          mutate(
            model = model %>% map(possibly(~ {
              # factors are required to calculate predictions
              # TODO: Refactor this to target hmsc_X
              res <- .x
              res$studyDesign <- res$studyDesign %>% mutate_if(is.character, as.factor)
              res
            }, NA)),
            gelman_diag = model %>% map(possibly(~ {
              .x %>%
                convertToCodaObject() %>%
                pluck("Beta") %>%
                gelman.diag(multivariate = FALSE)
            }, NA)),
            median_gelman_diag = gelman_diag %>% map_dbl(possibly(~ {
              .x %>%
                pluck("psrf") %>%
                as_tibble() %>%
                pull("Point est.") %>%
                median()
            }, NA)),
            pred = model %>% map(possibly(computePredictedValues, NA)),
            eval = pred %>% map2(model, possibly(~ {
              .x %>%
                evaluateModelFit(hM = .y) %>%
                map_dbl(~ .x %>% mean(na.rm = TRUE)) %>%
                enframe() %>%
                pivot_wider()
            }, NA))
          ) %>%
          unnest(eval) %>%
          select(-where(is.list))
      },
      dynamic = map(hmsc_models),
      hpc = TRUE
    ),

    hmsc_correlations = target(
      {
        hmsc_post_estimates %>%
          filter(parameter == "OmegaCor") %>%
          mutate(
            correlation = post_estimate %>% map(possibly(~ {
              res <-
                .x %>%
                filter(abs(mean) > 0.02) %>%
                mutate(
                  is_significant = support < hmsc_support_level | supportNeg < hmsc_support_level,
                  result = case_when(
                    mean > 0 & is_significant ~ "positive",
                    mean < 0 & is_significant ~ "negative"
                  )
                ) %>%
                # discard uncredible correlations
                filter(!is.na(result))

              taxa_pairs <-
                res$from_clean_taxon %>%
                union(res$to_clean_taxon) %>%
                combn(2) %>%
                t() %>%
                as_tibble() %>%
                rename(from_clean_taxon = V1, to_clean_taxon = V2)

              res <-
                res %>%
                # Only keep significant correlations between two different taxa
                inner_join(taxa_pairs) %>%
                filter(!is.na(result) & from_clean_taxon != to_clean_taxon) %>%
                left_join(hmsc_all_taxa %>% rename(from_clean_taxon = clean_taxon, from_taxon = taxon)) %>%
                left_join(hmsc_all_taxa %>% rename(to_clean_taxon = clean_taxon, to_taxon = taxon)) %>%
                left_join(lineages %>% select(from_taxon = taxon, from_kingdom = kingdom)) %>%
                left_join(lineages %>% select(to_taxon = taxon, to_kingdom = kingdom)) %>%
                mutate(
                  kingdom_group = case_when(
                    from_kingdom == "Bacteria" & to_kingdom == "Bacteria" ~ "within bacteria",
                    from_kingdom == "Fungi" & to_kingdom == "Fungi" ~ "within fungi",
                    from_kingdom != to_kingdom ~ "inter kingdom"
                  )
                )

              res
            }, NA))
          ) %>%
          select(-where(is.list), correlation)
      },
      hpc = FALSE,
    ),

    #' Phylogenetic signal in species niches
    #' wehere 0 is none and 1 is full influence
    hmsc_eval_phylogenetic_signals = target(
      {
        hmsc_models %>%
          mutate(
            phylogenetic_signal = model %>% map(possibly(~ {
              .x %>%
                convertToCodaObject() %>%
                pluck("Rho") %>%
                simplify() %>%
                as.double() %>%
                mean(na.rm = TRUE)
            }, NA))
          ) %>%
          select(-where(is.list), phylogenetic_signal)
      },
      dynamic = map(hmsc_models)
    ),

    hmsc_trace_plt = {
      param_set <- "Gamma"

      least_credible_params <-
        hmsc_model %>%
        convertToCodaObject() %>%
        pluck(param_set) %>%
        bayesplot::mcmc_intervals_data() %>%
        mutate(span = abs(hh - ll)) %>%
        arrange(-span) %>%
        pull(parameter) %>%
        head(12) %>%
        as.character()

      hmsc_model %>%
        convertToCodaObject() %>%
        pluck(param_set) %>%
        bayesplot::mcmc_trace(pars = least_credible_params) +
        labs(title = str_glue("MCMC Trace of least credible {param_set} parameters"))
    },

    hmsc_trait_samples_plt = {
      # extract longest parameter
      str_extract_samples_param <- ~ .x %>%
        str_extract(hmsc_model$X %>% colnames()) %>%
        discard(~ is.na(.x)) %>%
        last()

      str_extract_traits_param <- function(.x) {
        trait_groups_regex <-
          hmsc_model %>%
          pluck("TrFormula") %>%
          as.character() %>%
          last() %>%
          str_replace_all("\\+", "|") %>%
          str_remove_all("[ ]") %>%
          paste0("(", ., ")")

        .x %>%
          str_extract(hmsc_model$Tr %>% colnames()) %>%
          discard(~ is.na(.x)) %>%
          last() %>%
          str_remove(trait_groups_regex)
      }

      hmsc_model <- hmsc_models$model[[4]]

      # Do 90% of the posterior distribution have the same direction  of the effect?
      credibilty <-
        hmsc_model %>%
        convertToCodaObject() %>%
        pluck("Gamma") %>%
        bayesplot::mcmc_intervals_data(point_est = "mean", prob = 0.90, prob_outer = 0.95) %>%
        transmute(
          # extract longest parameter
          trait_param = parameter %>% map_chr(str_extract_traits_param),
          samples_param = parameter %>% map_chr(str_extract_samples_param),

          credibility_label = case_when(
            sign(ll) == sign(hh) ~ "**",
            sign(l) == sign(h) ~ "*",
            TRUE ~ ""
          )
        )

      hmsc_model %>%
        getPostEstimate(parName = "Gamma") %>%
        pluck("mean") %>%
        set_colnames(hmsc_model$Tr %>% colnames()) %>%
        set_rownames(hmsc_model$X %>% colnames()) %>%
        as_tibble(rownames = "samples_param") %>%
        pivot_longer(-samples_param, names_to = "trait_param", values_to = "estimate") %>%
        mutate(
          samples_param = samples_param %>% map_chr(str_extract_samples_param),
          trait_group = trait_param %>% str_extract("^(trophicMode|phylum|faprotax_entry)"),
          trait_param = trait_param %>% map_chr(str_extract_traits_param)
        ) %>%
        left_join(credibilty) %>%
        filter_at(c("samples_param", "trait_param", "trait_group"), ~ !is.na(.x)) %>%
        # filter(samples_param != "(Intercept)") %>%
        ggplot(aes(samples_param, trait_param)) +
        geom_tile(aes(fill = estimate)) +
        geom_text(aes(label = credibility_label), size = 8) +
        facet_grid(trait_group ~ "", scales = "free_y", space = "free") +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        scale_fill_gradient2(limits = c(-0.6, 0.6), low = muted("blue"), high = muted("red")) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(
          x = "Sample group",
          y = "Taxon trait",
          fill = TeX("$\\hat{\\gamma}$")
        )
    },

    hmsc_envirornment_group_kingdom_correlation_plot = {
      hmsc_correlations <-
        hmsc_models %>%
        filter(subset_name == "environment_group,kingdom" & samples_formula %>% str_detect("habitat")) %>%
        separate(subset_value, into = c("environment_group", "kingdom"), sep = ",") %>%
        distinct(response_name, environment_group, kingdom, model_type, model) %>%
        mutate(
          correlations = model %>% map(~ {
            hmsc_model <- .x
            post_estimate <- hmsc_model %>% Hmsc::getPostEstimate("OmegaCor", thin = 1)

            support <-
              post_estimate %>%
              pluck("support") %>%
              {
                .x <- .
                rownames(.x) <- hmsc_model$spNames
                colnames(.x) <- hmsc_model$spNames
                .x
              } %>%
              as_tibble(rownames = "from_clean") %>%
              pivot_longer(-from_clean, names_to = "to_clean", values_to = "support")

            post_estimate %>%
              pluck("mean") %>%
              {
                .x <- .
                rownames(.x) <- hmsc_model$spNames
                colnames(.x) <- hmsc_model$spNames
                .x
              } %>%
              as_tibble(rownames = "from_clean") %>%
              pivot_longer(-from_clean, names_to = "to_clean", values_to = "estimate") %>%
              left_join(support)
          })
        ) %>%
        select(-model) %>%
        unnest(correlations) %>%
        pivot_wider(
          names_from = model_type,
          values_from = c(estimate, support),
          values_fil = list(estimate = 0)
        )

      #' Keep onyl one record per pair and arrange taxa by similarity
      triangulize_taxa <- function(data, cluster_by = "abundance", model_type = "full") {
        estimate_name <- paste0("estimate_", model_type)

        ordered_taxa <-
          data %>%
          filter(response_name == cluster_by) %>%
          select_at(c("from_clean", "to_clean", estimate_name)) %>%
          pivot_wider(names_from = to_clean, values_from = estimate_name) %>%
          column_to_rownames("from_clean") %>%
          as.matrix() %>%
          dist() %>%
          hclust() %>%
          {
            .x <- .
            .x$labels[.x$order]
          }

        taxa_pairs <-
          ordered_taxa %>%
          combn(2) %>%
          t() %>%
          as_tibble() %>%
          rename(from_clean = V1, to_clean = V2) %>%
          mutate_at(c("from_clean", "to_clean"), ~ factor(.x, levels = ordered_taxa))

        data %>%
          mutate_at(c("from_clean", "to_clean"), ~ factor(.x, levels = ordered_taxa)) %>%
          arrange(from_clean, to_clean) %>%
          inner_join(taxa_pairs)
      }

      hmsc_correlations %>%
        # Significance filter
        # mutate(
        #   estimate_full = estimate_full * (support_full < hmsc_support_level | support_full > 1 - hmsc_support_level),
        #   estimate_null = estimate_null * (support_null < hmsc_support_level | support_null > 1 - hmsc_support_level),
        # ) %>%
        nest(-c("environment_group", "kingdom")) %>%
        mutate(
          plt = data %>% map(~ {
            .x %>%
              triangulize_taxa(model_type = "full") %>%
              pivot_wider(names_from = "response_name", values_from = c("estimate_null", "estimate_full", "support_null", "support_full")) %>%
              ggasym::asymmetrise(from_clean, to_clean) %>%
              ggplot(aes(from_clean, to_clean)) +
              ggasym::geom_asymmat(aes(
                fill_tl = estimate_full_abundance,
                fill_br = estimate_full_occurrence,
              ),
              fill_diag = "black"
              ) +
              ggasym::scale_fill_tl_gradientn(colors = vik_colors, limits = c(-1, 1)) +
              ggasym::scale_fill_br_gradientn(colors = cork_colors, limits = c(-1, 1)) +
              coord_fixed() +
              theme_void() +
              theme(
                axis.text = element_blank(),
                axis.ticks = element_blank(),

                strip.text = element_text(size = 20),
                axis.title = element_text(size = 20),

                plot.margin = unit(c(0, 5, 0, 0), "mm")
              ) +
              labs(x = "", y = "", fill_tl = "Coabundance (r)", fill_br = "Cooccurrence (r)")
          })
        ) %>%
        rename(Environment = environment_group, Kingdom = kingdom) %>%
        wrap_plots_grid(Kingdom ~ Environment, plot_column = "plt")
    },

    hmsc_euler_alpha = 0.7,
    hmsc_euler_legend_plt = {
      (
        habitats %>%
          mutate(habitat = habitat %>% factor(names(habitat_colors))) %>%
          filter(!is.na(habitat)) %>%
          ggplot(aes(habitat, fill = habitat)) +
          geom_bar(alpha = hmsc_euler_alpha) +
          scale_fill_habitat()
      ) %>%
        cowplot::get_legend()
    },

    hmsc_habitat_separated_kingdom_euler_plt = {
      correlations <-
        hmsc_post_estimates %>%
        filter(parameter == "OmegaCor" & subset_name == "habitat,kingdom") %>%
        separate(subset_value, into = c("habitat", "kingdom"), sep = ",") %>%
        rename(correlation = post_estimate) %>%
        select(-where(is.list), correlation) %>%
        unnest(correlation) %>%
        mutate(
          is_significant = support < hmsc_support_level | supportNeg < hmsc_support_level,
          result = case_when(
            mean > 0 & is_significant ~ "positive",
            mean < 0 & is_significant ~ "negative"
          )
        )

      taxa_pairs <-
        correlations$from_clean_taxon %>%
        union(correlations$to_clean_taxon) %>%
        combn(2) %>%
        t() %>%
        as_tibble() %>%
        rename(from_clean_taxon = V1, to_clean_taxon = V2)

      correlations <-
        correlations %>%
        # Only keep significant correlations between two different taxa
        inner_join(taxa_pairs) %>%
        filter(!is.na(result) & from_clean_taxon != to_clean_taxon) %>%
        select(response_name, habitat, from_clean_taxon, to_clean_taxon, result, kingdom) %>%
        left_join(habitats)

      correlations %>%
        nest(-c(response_name, environment_group, kingdom)) %>%
        mutate(
          plt = data %>% map(possibly(~ {
            res <-
              .x %>%
              unite("correlation", c(from_clean_taxon, to_clean_taxon, result)) %>%
              distinct(habitat, correlation) %>%
              nest(-habitat) %>%
              mutate(habitat = habitat %>% factor(levels = habitats$habitat) %>% map_chr(abbreviate)) %>%
              deframe() %>%
              map(~ .x %>% pluck("correlation"))

            # Note: Ensure distinct abbreviations
            colors <- habitat_colors
            names(colors) <- names(colors) %>% map_chr(abbreviate)
            colors <- colors[names(res)]

            res %>%
              eulerr::euler(shape = "circle") %>%
              plot(quantities = TRUE, labels = TRUE, fills = colors, alpha = hmsc_euler_alpha)
          }, NA))
        ) %>%
        filter(!is.na(kingdom)) %>%
        nest(-response_name) %>%
        mutate(
          plt = data %>% map2(response_name, ~ {
            plt <-
              .x %>%
              filter(!is.na(plt)) %>%
              rename(Environment = environment_group, `Kingdom group` = kingdom) %>%
              wrap_plots_grid(formula = formula(Environment ~ `Kingdom group`), plot_column = "plt")


            (plt | hmsc_euler_legend_plt) +
              plot_layout(widths = c(5, 1)) +
              plot_annotation(title = str_glue("Hmsc co-{.y}"))
          })
        ) %>%
        select(response_name, plt) %>%
        deframe()
    },

    hmsc_habitat_inter_intra_kingdom_euler_plts = {
      correlations <-
        hmsc_post_estimates %>%
        filter(parameter == "OmegaCor" & subset_name == "habitat") %>%
        rename(correlation = post_estimate, habitat = subset_value) %>%
        select(-where(is.list), correlation) %>%
        unnest(correlation) %>%
        mutate(
          is_significant = support < hmsc_support_level | supportNeg < hmsc_support_level,
          result = case_when(
            mean > 0 & is_significant ~ "positive",
            mean < 0 & is_significant ~ "negative"
          )
        )

      taxa_pairs <-
        correlations$from_clean_taxon %>%
        union(correlations$to_clean_taxon) %>%
        combn(2) %>%
        t() %>%
        as_tibble() %>%
        rename(from_clean_taxon = V1, to_clean_taxon = V2)

      correlations <-
        correlations %>%
        # Only keep significant correlations between two different taxa
        inner_join(taxa_pairs) %>%
        filter(!is.na(result) & from_clean_taxon != to_clean_taxon) %>%
        select(response_name, habitat, from_clean_taxon, to_clean_taxon, result) %>%

        # annotate
        left_join(habitats) %>%
        left_join(hmsc_all_taxa %>% rename(from_clean_taxon = clean_taxon, from_taxon = taxon)) %>%
        left_join(hmsc_all_taxa %>% rename(to_clean_taxon = clean_taxon, to_taxon = taxon)) %>%
        left_join(lineages %>% select(from_taxon = taxon, from_kingdom = kingdom)) %>%
        left_join(lineages %>% select(to_taxon = taxon, to_kingdom = kingdom)) %>%
        mutate(
          kingdom_group = case_when(
            from_kingdom == "Bacteria" & to_kingdom == "Bacteria" ~ "within bacteria",
            from_kingdom == "Fungi" & to_kingdom == "Fungi" ~ "within fungi",
            from_kingdom != to_kingdom ~ "inter kingdom"
          )
        )

      correlations %>%
        nest(-c(response_name, environment_group, kingdom_group)) %>%
        mutate(
          plt = data %>% map(possibly(~ {
            res <-
              .x %>%
              unite("correlation", c(from_clean_taxon, to_clean_taxon, result)) %>%
              distinct(habitat, correlation) %>%
              nest(-habitat) %>%
              mutate(habitat = habitat %>% factor(levels = habitats$habitat) %>% map_chr(abbreviate)) %>%
              deframe() %>%
              map(~ .x %>% pluck("correlation"))

            # Note: Ensure distinct abbreviations
            colors <- habitat_colors
            names(colors) <- names(colors) %>% map_chr(abbreviate)
            colors <- colors[names(res)]

            res %>%
              eulerr::euler(shape = "circle") %>%
              plot(quantities = TRUE, labels = TRUE, fills = colors, alpha = hmsc_euler_alpha)
          }, NA))
        ) %>%
        filter(!is.na(kingdom_group)) %>%
        nest(-response_name) %>%
        mutate(
          plt = data %>% map2(response_name, ~ {
            plt <-
              .x %>%
              filter(!is.na(plt)) %>%
              rename(Environment = environment_group, `Kingdom group` = kingdom_group) %>%
              wrap_plots_grid(formula = formula(Environment ~ `Kingdom group`), plot_column = "plt")


            (plt | hmsc_euler_legend_plt) +
              plot_layout(widths = c(5, 1)) +
              plot_annotation(title = str_glue("Hmsc co-{.y}"))
          })
        ) %>%
        select(response_name, plt) %>%
        deframe()
    },

    hmsc_environment_group_inter_intra_kingdom_euler_plts = {
      correlations <-
        hmsc_post_estimates %>%
        filter(parameter == "OmegaCor" & subset_name == "environment_group") %>%
        rename(correlation = post_estimate, environment_group = subset_value) %>%
        distinct(response_name, correlation, environment_group, .keep_all = TRUE) %>%
        select(-where(is.list), correlation) %>%
        unnest(correlation) %>%
        mutate(
          is_significant = support < hmsc_support_level | supportNeg < hmsc_support_level,
          result = case_when(
            mean > 0 & is_significant ~ "positive",
            mean < 0 & is_significant ~ "negative"
          )
        )

      taxa_pairs <-
        correlations$from_clean_taxon %>%
        union(correlations$to_clean_taxon) %>%
        combn(2) %>%
        t() %>%
        as_tibble() %>%
        rename(from_clean_taxon = V1, to_clean_taxon = V2)

      correlations <-
        correlations %>%
        # Only keep significant correlations between two different taxa
        inner_join(taxa_pairs) %>%
        filter(!is.na(result) & from_clean_taxon != to_clean_taxon) %>%
        select(response_name, environment_group, from_clean_taxon, to_clean_taxon, result) %>%

        # annotate
        left_join(hmsc_all_taxa %>% rename(from_clean_taxon = clean_taxon, from_taxon = taxon)) %>%
        left_join(hmsc_all_taxa %>% rename(to_clean_taxon = clean_taxon, to_taxon = taxon)) %>%
        left_join(lineages %>% select(from_taxon = taxon, from_kingdom = kingdom)) %>%
        left_join(lineages %>% select(to_taxon = taxon, to_kingdom = kingdom)) %>%
        mutate(
          kingdom_group = case_when(
            from_kingdom == "Bacteria" & to_kingdom == "Bacteria" ~ "within bacteria",
            from_kingdom == "Fungi" & to_kingdom == "Fungi" ~ "within fungi",
            from_kingdom != to_kingdom ~ "inter kingdom"
          )
        )

      plt <-
        correlations %>%
        nest(-c(response_name, kingdom_group)) %>%
        mutate(
          plt = data %>% map(possibly(~ {
            res <-
              .x %>%
              unite("correlation", c(from_clean_taxon, to_clean_taxon, result)) %>%
              distinct(environment_group, correlation) %>%
              group_by(environment_group) %>%
              nest() %>%
              deframe() %>%
              map(~ .x %>% pluck("correlation"))

            res %>%
              eulerr::euler(shape = "circle") %>%
              # ensure order of colors
              plot(
                quantities = TRUE, labels = FALSE,
                fills = environment_group_colors[names(res)], alpha = hmsc_euler_alpha
              )
          }, NA))
        ) %>%
        filter(!is.na(kingdom_group)) %>%
        rename(`Response` = response_name, `Kingdom group` = kingdom_group) %>%
        wrap_plots_grid(Response ~ `Kingdom group`, plot_column = "plt")

      (plt | hmsc_euler_legend_plt) +
        plot_layout(widths = c(5, 1))
    },

    hmsc_var_part_plt = {
      hmsc_varpart <-
        hmsc_models %>%
        filter(
          samples_formula %>% str_detect("habitat") &
            model_type %in% c("full", "null") &
            traits_formula == "~ kingdom + trophicMode + phylum + faprotax_entry - 1"
        ) %>%
        transmute(
          response_name,
          varpart = model %>% map(function(hmsc_model) {
            hmsc_model_habitats <-
              hmsc_model %>%
              pluck("X") %>%
              colnames() %>%
              map_chr(~ .x %>% str_remove("^habitat")) %>%
              discard(~ .x == "seqdepth") %>%
              tibble(habitat = .) %>%
              left_join(habitats) %>%
              nest(-environment_group) %>%
              mutate(environment_group_id = row_number() + 1) %>%
              unnest(data)

            # Taxa vs samples
            hmsc_vp <-
              hmsc_model %>%
              computeVariancePartitioning(
                group = {
                  hmsc_model_habitats %>%
                    pull(environment_group_id) %>%
                    c(1, .)
                },
                groupnames = c("seqdepth", hmsc_model_habitats$environment_group %>% unique())
              ) %>%
              pluck("vals") %>%
              as_tibble(rownames = "variable") %>%
              pivot_longer(-variable, names_to = "clean_taxon", values_to = "var_proportion") %>%
              mutate(
                variable = variable %>% str_remove("^environment_group") %>% factor(levels = environment_groups)
              ) %>%
              left_join(hmsc_taxa) %>%
              left_join(lineages) %>%

              # pool subgroups
              group_by(kingdom, taxon, variable) %>%
              summarise(var_proportion = sum(var_proportion)) %>%

              # re scale, only look at fixed effects
              filter(!is.na(variable)) %>%
              group_by(taxon) %>%
              mutate(var_proportion = var_proportion / sum(var_proportion))
          })
        )

      hmsc_varpart %>%
        unnest(varpart) %>%
        select(response_name, variable, taxon, var_proportion) %>%
        pivot_wider(names_from = variable, values_from = var_proportion) %>%
        left_join(lineages) %>%
        # annotate outliers
        dplyr::mutate(
          # label = taxon %>% map_chr(~ ifelse(.x %in% extreme_taxa, .x, NA)),
          phylum = phylum %>% map_chr(~ ifelse(.x %in% names(phyla_colors), .x, "other"))
        ) %>%
        ggplot(aes(x = soil, y = aquatic, z = host)) +
        geom_point(aes(color = phylum)) +
        # geom_point(data = tibble(soil = 1/3, aquatic = 1/3, host = 1/3), color = "lightgrey") +
        # geom_text(aes(label = label), size = 3) +
        scale_color_phyla() +
        ggtern::coord_tern() +
        facet_grid(response_name ~ kingdom) +
        ggtern::theme_custom(
          col.T = environment_group_colors[["aquatic"]],
          col.R = environment_group_colors[["host"]],
          col.L = environment_group_colors[["soil"]],
          col.grid.minor = "grey"
        ) +
        ggplot2::theme(
          panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 18),
          plot.title = element_text(size = 20)
        ) +
        labs(
          # title = "Partitioning variance",
          # subtitle = str_glue(
          #   hmsc_model %>% capture.output() %>% paste0(collapse = "\n"),
          #   hmsc_model$XFormula %>% format(),
          #   hmsc_model$TrFormula %>% format(),
          #   .sep = "\n"
          # ),
          color = "Phylum"
        )

      # majority_variables <-
      #   hmsc_vp %>%
      #   group_by(taxon) %>%
      #   filter(! is.na(variable)) %>%
      #   arrange(-var_proportion) %>%
      #   slice(1) %>%
      #   distinct(taxon, majority_variable = variable)
      #
      # hmsc_vp %>%
      #   ggplot(aes(variable, var_proportion, fill = variable)) +
      #   geom_boxplot() +
      #   scale_fill_environment_group() +
      #   scale_y_continuous(expand = c(0, 0))
      #
      # majority_variables %>%
      #   left_join(lineages) %>%
      #   ungroup() %>%
      #   count(kingdom, majority_variable) %>%
      #   ggplot(aes(kingdom, n)) +
      #     geom_bar(aes(fill = majority_variable), stat = "identity") +
      #     scale_fill_environment_group()
    },

    hmsc_pca_taxa_habitats_plt = {
      hmsc_models %>%
        filter(subset_name == "environment_group,kingdom" & model_type == "residual_habitat" & response_name == "abundance") %>%
        expand_grid(parameter = "Beta") %>% # only this needed for pca plot
        mutate(
          post_estimate = model %>% map2(parameter, possibly(function(model, parameter) {
            model %>%
              # this can take many minutes
              getPostEstimate(parameter) %>%
              # for mean, support, and supportNeg
              map(function(m) {
                switch(parameter,
                  "OmegaCor" = {
                    rownames(m) <- model$spNames
                    colnames(m) <- model$spNames
                    m %>%
                      as_tibble(rownames = "from_clean_taxon") %>%
                      pivot_longer(-from_clean_taxon, names_to = "to_clean_taxon")
                  },
                  "Beta" = {
                    rownames(m) <- model$covNames

                    m %>%
                      as_tibble(rownames = "covariate") %>%
                      pivot_longer(-covariate, names_to = "clean_taxon")
                  },
                  # generic result
                  {
                    m %>%
                      as_tibble(rownames = "from") %>%
                      pivot_longer(-from, names_to = "to")
                  }
                )
              }) %>%
              enframe() %>%
              unnest(value) %>%
              pivot_wider(names_from = name, values_from = value)
          }, NA))
        ) %>%
        select(-where(is.list), post_estimate) %>%
        separate(subset_value, into = c("environment_group", "kingdom"), sep = ",") %>%
        # plot
        mutate(
          plt = list(post_estimate, environment_group, kingdom, subset_value) %>% pmap(function(post_estimate, environment_group, kingdom, subset_value) {
            pca <-
              post_estimate %>%
              select(covariate, clean_taxon, mean) %>%
              pivot_wider(names_from = clean_taxon, values_from = mean) %>%
              column_to_rownames("covariate") %>%
              prcomp()

            auto_plt <- pca %>% autoplot(loadings = TRUE)


            tibble(
              facet_x = environment_group,
              facet_y = kingdom
            ) %>%
              ggplot(aes(PC1, PC2)) +
              geom_point(
                mapping = aes(color = phylum),
                data = {
                  pca$rotation %>%
                    as_tibble(rownames = "clean_taxon") %>%
                    left_join(hmsc_all_taxa) %>%
                    left_join(lineages) %>%
                    mutate(phylum = phylum %>% map_chr(~ ifelse(.x %in% names(phyla_colors), .x, NA)))
                }
              ) +
              scale_color_phyla() +
              new_scale_color() +
              geom_point(
                mapping = aes(color = habitat),
                data = {
                  auto_plt %>%
                    pluck("data") %>%
                    as_tibble(rownames = "covariate") %>%
                    mutate(
                      habitat = covariate %>%
                        str_extract(habitats$habitat %>% paste0(collapse = "|")) %>%
                        factor(levels = names(habitat_colors))
                    )
                },
                shape = 24,
                size = 2
              ) +
              facet_grid(facet_x ~ facet_y) +
              scale_color_habitat() +
              labs(
                x = auto_plt$labels$x,
                y = auto_plt$labels$y
              )
          })
        ) %>%
        arrange(subset_value, response_name) %>%
        pull(plt) %>%
        wrap_plots(ncol = 2, guides = "collect")
    },

    #
    # Do the bacterial traits depend on the fungal traits in inter kingdom correlations in soil?
    #
    hmsc_soil_correlation_inter_kingdom_traits_dependencies = {
      hmsc_correlations %>%
        filter(subset_name == "habitat") %>%
        rename(habitat = subset_value) %>%
        left_join(habitats) %>%
        filter(environment_group == "soil") %>%
        unnest(correlation) %>%
        filter(kingdom_group == "inter kingdom") %>%

        # Unify, because the fungus can be either from_taxon or to_taxon
        transmute(
          habitat,
          result,
          fungus_taxon = case_when(
            from_kingdom == "Fungi" ~ from_taxon,
            to_kingdom == "Fungi" ~ to_taxon
          ),
          bacteria_taxon = case_when(
            from_kingdom == "Bacteria" ~ from_taxon,
            to_kingdom == "Bacteria" ~ to_taxon
          )
        ) %>%
        left_join(
          traits %>% select(fungus_taxon = taxon, fungus_guild = guild, fungus_trophic_mode = trophicMode)
        ) %>%
        left_join(
          traits %>% select(bacteria_taxon = taxon, bacteria_trait = faprotax_entry)
        ) %>%
        distinct(result, fungus_taxon, fungus_guild, fungus_trophic_mode, bacteria_taxon, bacteria_trait) %>%
        filter(
          !is.na(fungus_guild) &
            !is.na(fungus_trophic_mode) &
            !is.na(bacteria_trait)
        ) %>%
        nest(-c(result, fungus_trophic_mode)) %>%
        mutate(
          test = data %>% map(~ {
            .x %>%
              count(fungus_guild, bacteria_trait) %>%
              pivot_wider(names_from = bacteria_trait, values_from = n, values_fill = list(n = 0)) %>%
              column_to_rownames("fungus_guild") %>%
              as.matrix() %>%
              chisq.test() %>%
              tidy()
          })
        ) %>%
        unnest(test) %>%
        ungroup() %>%
        mutate(q.value = p.value %>% p.adjust()) %>%
        select(-where(is.list))
    },

    hmsc_soil_correlations_plt = {
      correlations <-
        hmsc_post_estimates %>%
        filter(parameter == "OmegaCor" & subset_name == "habitat" & model_type == "null") %>%
        rename(correlation = post_estimate, habitat = subset_value) %>%
        select(-where(is.list), correlation) %>%
        unnest(correlation) %>%
        mutate(
          is_significant = support < hmsc_support_level | supportNeg < hmsc_support_level,
          result = case_when(
            mean > 0 & is_significant ~ "positive",
            mean < 0 & is_significant ~ "negative"
          )
        )

      taxa_pairs <-
        correlations$from_clean_taxon %>%
        union(correlations$to_clean_taxon) %>%
        combn(2) %>%
        t() %>%
        as_tibble() %>%
        rename(from_clean_taxon = V1, to_clean_taxon = V2)

      correlations <-
        correlations %>%
        # Only keep significant correlations between two different taxa
        inner_join(taxa_pairs) %>%
        filter(!is.na(result) & from_clean_taxon != to_clean_taxon) %>%
        select(response_name, habitat, from_clean_taxon, to_clean_taxon, result, mean, support) %>%

        # annotate
        left_join(habitats) %>%
        left_join(hmsc_all_taxa %>% rename(from_clean_taxon = clean_taxon, from_taxon = taxon)) %>%
        left_join(hmsc_all_taxa %>% rename(to_clean_taxon = clean_taxon, to_taxon = taxon)) %>%
        left_join(lineages %>% select(from_taxon = taxon, from_kingdom = kingdom)) %>%
        left_join(lineages %>% select(to_taxon = taxon, to_kingdom = kingdom)) %>%
        mutate(
          kingdom_group = case_when(
            from_kingdom == "Bacteria" & to_kingdom == "Bacteria" ~ "within bacteria",
            from_kingdom == "Fungi" & to_kingdom == "Fungi" ~ "within fungi",
            from_kingdom != to_kingdom ~ "inter kingdom"
          )
        ) %>%
        # restrict to soil habitats
        left_join(habitats) %>%
        filter(environment_group == "soil")

      topology_plts <-
        correlations %>%
        nest(-c(response_name, habitat)) %>%
        mutate(
          graph = data %>% map(~ {
            .x %>%
              rename(from = from_taxon, to = to_taxon, estimate = mean) %>%
              as_tbl_graph()
          }),

          topology = graph %>% map(~ {
            .x %>%
              activate(nodes) %>%
              as_tibble() %>%
              pivot_longer(
                cols = c(degree, closeness, betweeness),
                names_to = "topology_name",
                values_to = "topology_value"
              )
          })
        ) %>%
        unnest(topology) %>%
        nest(-topology_name) %>%
        mutate(
          plt = data %>% map2(topology_name, ~ {
            .x %>%
              ggplot(aes(habitat, topology_value, color = habitat)) +
              geom_boxplot() +
              stat_compare_means(comparisons = list(
                c("conifer forests", "temperate forests"),
                c("tropical forests", "temperate forests"),
                c("conifer forests", "tropical forests")
              )) +
              scale_color_habitat() +
              theme(
                panel.grid.major.x = element_blank(),
                axis.text.x = element_blank()
              ) +
              labs(y = .y)
          })
        ) %>%
        pull(plt) %>%
        wrap_plots(nrow = 1, guides = "collect")

      euler_plts <-
        correlations %>%
        filter(response_name == "abundance" & habitat %in% {
          habitats %>%
            filter(environment_group == "soil") %>%
            pull(habitat)
        }) %>%
        nest(-c(kingdom_group)) %>%
        arrange(kingdom_group) %>%
        mutate(
          plt = data %>% map2(kingdom_group, possibly(~ {
            res <-
              .x %>%
              unite("correlation", c(from_clean_taxon, to_clean_taxon, result)) %>%
              distinct(habitat, correlation) %>%
              group_by(habitat) %>%
              nest() %>%
              deframe() %>%
              map(~ .x %>% pluck("correlation"))

            res %>%
              eulerr::euler(shape = "circle") %>%
              # ensure order of colors
              plot(labels = FALSE, fills = habitat_colors[names(res)], main = list(label = .y, cex = 0.7), quantities = list(cex = 0.7), alpha = hmsc_euler_alpha)
          }, NA))
        ) %>%
        pull(plt) %>%
        wrap_plots(ncol = 1)

      (euler_plts | topology_plts) +
        plot_layout(widths = c(1, 4)) +
        plot_annotation(tag_levels = "A")
    },

    hmsc_environment_group_correlations_plt = target(
      {
        correlations <-
          hmsc_post_estimates %>%
          filter(parameter == "OmegaCor" & subset_name == "environment_group,kingdom" & model_type == "null" & response_name == "abundance") %>%
          separate(subset_value, into = c("environment_group", "kingdom"), sep = ",") %>%
          rename(correlation = post_estimate) %>%
          select(-where(is.list), correlation) %>%
          unnest(correlation) %>%
          mutate(
            is_significant = support < hmsc_support_level | supportNeg < hmsc_support_level,
            result = case_when(
              mean > 0 & is_significant ~ "positive",
              mean < 0 & is_significant ~ "negative"
            )
          )

        taxa_pairs <-
          correlations$from_clean_taxon %>%
          union(correlations$to_clean_taxon) %>%
          combn(2) %>%
          t() %>%
          as_tibble() %>%
          rename(from_clean_taxon = V1, to_clean_taxon = V2)

        correlations <-
          correlations %>%
          # Only keep significant correlations between two different taxa
          inner_join(taxa_pairs) %>%
          filter(!is.na(result) & from_clean_taxon != to_clean_taxon) %>%
          select(response_name, environment_group, kingdom, from_clean_taxon, to_clean_taxon, result, mean, support) %>%

          # annotate
          left_join(hmsc_all_taxa %>% rename(from_clean_taxon = clean_taxon, from_taxon = taxon)) %>%
          left_join(hmsc_all_taxa %>% rename(to_clean_taxon = clean_taxon, to_taxon = taxon)) %>%
          left_join(lineages %>% select(from_taxon = taxon, from_kingdom = kingdom)) %>%
          left_join(lineages %>% select(to_taxon = taxon, to_kingdom = kingdom)) %>%
          mutate(
            kingdom_group = case_when(
              from_kingdom == "Bacteria" & to_kingdom == "Bacteria" ~ "within bacteria",
              from_kingdom == "Fungi" & to_kingdom == "Fungi" ~ "within fungi",
              from_kingdom != to_kingdom ~ "inter kingdom"
            )
          )

        topology_plts <-
          correlations %>%
          nest(-c(response_name, environment_group, kingdom)) %>%
          mutate(
            graph = data %>% map(~ {
              .x %>%
                rename(from = from_taxon, to = to_taxon, estimate = mean) %>%
                as_tbl_graph()
            }),

            topology = graph %>% map(~ {
              .x %>%
                activate(nodes) %>%
                as_tibble() %>%
                pivot_longer(
                  cols = c(degree, closeness, betweeness),
                  names_to = "topology_name",
                  values_to = "topology_value"
                )
            })
          ) %>%
          unnest(topology) %>%
          nest(-c(topology_name, kingdom)) %>%
          mutate(
            plt = list(data, topology_name, kingdom) %>% pmap(~ {
              ..1 %>%
                mutate(environment_group = environment_group %>% factor(levels = c("aquatic", "host", "soil"))) %>%
                ggplot(aes(environment_group, topology_value, color = environment_group)) +
                geom_boxplot() +
                stat_compare_means_environment_group() +
                scale_color_environment_group() +
                theme(
                  panel.grid.major.x = element_blank(),
                  axis.text.x = element_blank()
                ) +
                facet_wrap(~ str_glue("{..3}")) +
                labs(y = ..2, x = "Environment", color = "Environment")
            })
          ) %>%
          arrange(topology_name, kingdom) %>%
          pull(plt) %>%
          wrap_plots(ncol = 2, guides = "collect")

        euler_plts <-
          correlations %>%
          nest(-c(kingdom)) %>%
          arrange(kingdom) %>%
          mutate(
            plt = data %>% map2(kingdom, possibly(~ {
              res <-
                .x %>%
                unite("correlation", c(from_clean_taxon, to_clean_taxon, result)) %>%
                distinct(environment_group, correlation) %>%
                group_by(environment_group) %>%
                nest() %>%
                deframe() %>%
                map(~ .x %>% pluck("correlation"))

              res %>%
                eulerr::euler(shape = "circle") %>%
                # ensure order of colors
                plot(
                  labels = FALSE, fills = environment_group_colors[names(res)],
                  main = list(label = str_glue("{.y}"), cex = 0.7),
                  quantities = list(cex = 0.7), alpha = hmsc_euler_alpha
                )
            }, NA))
          ) %>%
          pull(plt) %>%
          wrap_plots(ncol = 1)

        (euler_plts | topology_plts) +
          plot_layout(widths = c(1, 4)) +
          plot_annotation(tag_levels = "A")
      },
      hpc = FALSE
    ),

    hmsc_residual_correlations_plt = target(
      {
        hmsc_correlations <-
          hmsc_models %>%
          filter(subset_name == "environment_group,kingdom") %>%
          mutate(
            correlation = model %>% map(get_hmsc_correlations)
          ) %>%
          select(-where(is.list), correlation)

        influence_plt <-
          hmsc_correlations %>%
          unnest(correlation) %>%
          transmute(
            response_name, subset_value, model_type, from_clean_taxon, to_clean_taxon,
            result = case_when(
              mean > 0 & (support < hmsc_support_level | supportNeg > hmsc_support_level) ~ "positive",
              mean < 0 & (support < hmsc_support_level | supportNeg > hmsc_support_level) ~ "negative",
              TRUE ~ "n.s."
            )
          ) %>%
          pivot_wider(
            names_from = model_type,
            values_from = result
          ) %>%
          mutate(habitat_influence = residual_habitat != null) %>%
          separate(subset_value, into = c("environment_group", "kingdom")) %>%
          group_by(response_name, environment_group, kingdom, habitat_influence) %>%
          count() %>%
          pivot_wider(names_from = habitat_influence, values_from = n, values_fill = list(n = 0)) %>%
          mutate(perc_habitat_influenced = `TRUE` / (`TRUE` + `FALSE`) * 100) %>%
          ggplot(aes(environment_group, perc_habitat_influenced, fill = environment_group)) +
          geom_bar(stat = "identity") +
          scale_fill_environment_group() +
          scale_y_continuous(expand = c(0, 0)) +
          theme(
            axis.text.x = element_blank(),
            panel.grid.major.x = element_blank()
          ) +
          facet_grid("" ~ kingdom + response_name) +
          labs(
            x = "Environment",
            fill = "Environment",
            y = "Interactions influenced by habitat (%)"
          )

        influence_plt
      },
      hpc = FALSE
    )
  )
}
