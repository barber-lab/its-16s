get_ml_plan <- function() {
  drake::drake_plan(
    ml_min_prevalence_perc = 20,
    ml_n_tuning_folds = 3,

    tuning_metrics = tribble(
      ~mode, ~tuning_metric,
      "classification", "roc_auc",
      "regression", "rmse"
    ),

    model_specs = {
      tribble(
        ~spec_name, ~spec,

        "null", null_model() %>%
          set_engine("parsnip"),

        "linear_reg", linear_reg(
          penalty = tune(),
          mixture = tune()
        ) %>%
          set_engine("glmnet"),

        "rand_forest", rand_forest(
          mtry = tune(),
          trees = 1000, # not worth to tune. Just needs to be high enough
          min_n = tune()
        ) %>%
          set_engine("ranger", importance = "impurity"),

        "nearest_neighbor", nearest_neighbor(
          neighbors = tune(),
          weight_func = tune(),
          dist_power = tune()
        ) %>%
          set_engine("kknn"),

        "svm rbf", svm_rbf(
          cost = tune(),
          rbf_sigma = tune()
        ) %>%
          set_engine("kernlab"),

        "mlp", mlp(
          hidden_units = 100 # must be low enough
        ) %>%
          set_engine("nnet")
      ) %>%
        expand_grid(outcomes %>% select(-where(is.list))) %>%
        mutate(spec = spec %>% map2(mode, ~ set_mode(.x, .y))) %>%
        relocate(spec, .after = last_col()) %>%

        # remove unsupported combinations
        anti_join(tribble(
          ~spec_name, ~mode,
          "linear_reg", "classification"
        ))
    },

    outcomes = {
      expand_grid(
        outcome_kingdom = kingdoms_groups,
        predictor_kingdom = kingdoms_groups,
        mode = c("classification", "regression")
      ) %>%
        mutate(
          outcome = mode %>% map2_chr(outcome_kingdom, ~ {
            switch(.x,
              "classification" = .y %>% tolower() %>% paste0("dysbalance_", .),
              "regression" =  .y %>% tolower() %>% paste0("distance_to_medoid_", .)
            )
          }),
          predictors = predictor_kingdom %>% map(~ switch(.x,
            "Fungi" = all_predictors$Fungi,
            "Bacteria" = all_predictors$Bacteria
          ))
        )
    },

    features_meta_tbl = {
      abundances %>%
        filter(norm_method == "raw") %>%
        pull(data) %>%
        first() %>%
        pull(taxon) %>%
        unique() %>%
        tibble(taxon = .) %>%
        mutate(feature = taxon %>% janitor::make_clean_names())
    },

    all_predictors = {
      features_meta_tbl %>%
        left_join(lineages) %>%
        group_by(kingdom) %>%
        nest() %>%
        mutate(data = data %>% map(~ .x %>%
          pull(feature) %>%
          intersect(features_all %>% colnames()))) %>%
        deframe()
    },

    features_all = {
      prevalent_taxa <-
        prevalences %>%
        filter(subset_name == "kingdom") %>%
        unnest(data) %>%
        filter(prevalence_perc >= ml_min_prevalence_perc) %>%
        pull(taxon) %>%
        unique()

      abundances %>%
        filter(norm_method == "raw") %>%
        pull(data) %>%
        first() %>%
        inner_join(
          dysbalances %>%
            select(sample_id, kingdom, distance_to_medoid, dysbalance) %>%
            pivot_wider(
              names_from = kingdom,
              values_from = c(distance_to_medoid, dysbalance),
            ) %>%
            filter(!is.na(dysbalance_Bacteria) & !is.na(dysbalance_Fungi)) %>%
            ungroup()
        ) %>%
        filter(taxon %in% prevalent_taxa) %>%
        ungroup() %>%
        pivot_wider(
          names_from = taxon,
          values_from = abundance,
          values_fill = list(abundance = 0)
        ) %>%
        mutate_if(is.character, as.factor) %>%
        # some methods e.g. ranger needs names w/o spaces
        janitor::clean_names()
    },

    # test hold out to report performance
    training_splits = target(
      {
        res <-
          features_all %>%
          rsample::initial_split(prop = 0.8, strata = outcomes$outcome[[1]])

        outcomes %>%
          select(-where(is.list)) %>%
          mutate(split = list(res))
      },
      dynamic = map(outcomes)
    ),

    recipes = target(
      {
        outcomes %>%
          left_join(training_splits) %>%
          mutate(
            recipe = list(outcome, split, predictors, mode) %>% pmap(function(outcome, split, predictors, mode) {
              split %>%
                training() %>%
                recipe(formula = paste0(outcome, "~.") %>% as.formula(), data = .) %>%
                update_role(all_predictors(), new_role = "etc") %>%
                update_role(predictors, new_role = "predictor") %>%
                update_role(sample_id, new_role = "id variable") %>%
                step_zv(all_predictors()) %>%
                {
                  .x <- .

                  if (mode == "classification") {
                    .x %>% themis::step_smote(outcome)
                  } else {
                    .x
                  }
                } %>%
                step_center(all_predictors(), -all_outcomes()) %>%
                step_scale(all_predictors(), -all_outcomes()) %>%
                step_corr(all_predictors(), method = "spearman", threshold = 0.95) %>%
                recipeselectors::step_select_boruta(all_predictors(), outcome = outcome)
              # recipeselectors::step_select_roc(all_predictors(), outcome = outcome, top_p = tune())
            })
          ) %>%
          select(!where(is.list), recipe)
      },
      trigger = trigger(condition = TRUE)
    ),

    workflows = target({
      model_specs %>%
        inner_join(recipes) %>%
        mutate(
          workflow = recipe %>% map2(spec, ~ workflows::workflow() %>%
            add_recipe(.x) %>%
            add_model(.y))
        ) %>%
        select(-where(is.list), workflow)
    }),

    tuneable_workflows = {
      workflows %>%
        filter(
          workflow %>% map_lgl(~ {
            .x %>%
              parameters() %>%
              as.data.frame() %>%
              nrow() > 0
          })
        )
    },

    untuneable_workflows = {
      # Do not consider list columns because filtering join gets confused
      workflows %>% anti_join(tuneable_workflows %>% select(-where(is.list)))
    },

    tuneable_workflows_by = {
      tuneable_workflows %>%
        nrow() %>%
        seq()
    },

    untuneable_workflows_by = {
      untuneable_workflows %>%
        nrow() %>%
        seq()
    },

    untuneable_models = target(
      {
        untuneable_workflows %>%
          inner_join(training_splits) %>%
          inner_join(tuning_metrics) %>%
          mutate(
            model = workflow %>% map2(
              split,
              possibly(~ .x %>% last_fit(split = .y), NA)
            )
          ) %>%
          select(-where(is.list), model)
      },
      dynamic = group(untuneable_workflows, .by = untuneable_workflows_by),
      hpc = TRUE
    ),

    tunings = target(
      {
        tuneable_workflows %>%
          inner_join(training_splits) %>%
          inner_join(tuning_metrics) %>%
          mutate(
            resample = split %>% map2(outcome, possibly(~ {
              .x %>%
                training() %>%
                vfold_cv(v = ml_n_tuning_folds, strata = .y)
            }, NA)),
            model = list(workflow, resample, tuning_metric) %>%
              pmap(possibly(function(workflow, resample, tuning_metric) {
                workflow %>%
                  tune_bayes(
                    resamples = resample,
                    param_info = workflow %>% parameters() %>% dials::finalize(x = resample),
                    metrics = tuning_metric %>% switch("roc_auc" = metric_set(roc_auc), "rmse" = metric_set(rmse)),
                    control = control_bayes(verbose = TRUE, time_limit = 4 * 60),
                    iter = 50
                  )
              }, NA))
          ) %>%
          select(-where(is.list), model)
      },
      dynamic = group(tuneable_workflows, .by = tuneable_workflows_by),
      hpc = TRUE
    ),

    tuned_models = target(
      {
        tunings %>%
          inner_join(workflows) %>%
          inner_join(training_splits) %>%
          inner_join(tuning_metrics) %>%
          mutate(
            model = list(model, workflow, split, tuning_metric) %>% pmap(possibly(function(tuning, workflow, split, tuning_metric) {
              best_tuning <-
                tuning %>%
                select_best(tuning_metric)

              workflow %>%
                finalize_workflow(best_tuning) %>%
                last_fit(split)
            }, NA))
          ) %>%
          select(-where(is.list), model)
      },
      dynamic = map(tunings),
      hpc = TRUE
    ),

    models = {
      tuned_models %>% bind_rows(untuneable_models)
    },

    best_models = target(
      {
        models %>%
          unnest(model) %>%
          unnest(.metrics) %>%
          select(-where(is.list)) %>%
          filter(tuning_metric == .metric) %>%
          group_by(mode, outcome, outcome_kingdom, predictor_kingdom) %>%
          # select lowest value but highest auc is the best one
          mutate(performance = .estimate %>% map2_dbl(tuning_metric, ~ {
            switch(.y, "rmse" = .x, "roc_auc" = .x * -1)
          })) %>%
          arrange(performance) %>%
          slice(1) %>%
          select(spec_name) %>% # select only meta data columns specifiying model
          inner_join(models) %>%
          ungroup()
      },
      trigger = trigger(change = models)
    ),

    # feature_importance_batches_data = target({
    #   n_batches <- (parallel::detectCores() / nrow(best_models)) %>% ceiling()
    #
    #   # paralllelize: Split model to different feature sets
    #   best_models %>%
    #     inner_join(outcomes) %>%
    #     unnest(predictors) %>%
    #     select(-where(is.list)) %>%
    #     mutate(batch = row_number() %% n_batches) %>%
    #     group_by(across(-predictors)) %>%
    #     summarise(predictors = predictors %>% list()) %>%
    #     ungroup() %>%
    #     inner_join(best_models) %>%
    #     inner_join(training_splits)
    # },
    # # explicitly prevent loading big target to all workers
    # hpc = FALSE
    # ),

    feature_importance_batches = target(
      {
        best_models %>%
          inner_join(outcomes) %>%
          inner_join(training_splits) %>%
          mutate(
            feature_importance = list(model, split, predictors) %>% pmap(possibly(~ {
              browser()

              vip::vi(
                object = ..1$.workflow[[1]] %>% pull_workflow_fit(),
                method = "firm", feature_names = ..3, train = ..2 %>% training(),
                parallel = TRUE # ... forwarded to pdp::partial if method == "firm"
              )
            }, NA))
          ) %>%
          select(-where(is.list), feature_importance)
      },
      dynamic = map(best_models),
      hpc = TRUE
    ),

    fig_dysbalance_feature_importance = {
      best_models %>%
        filter(mode == "regression" & spec_name == "rand_forest") %>%
        unnest(model) %>%
        transmute(
          predictor_kingdom,
          outcome_kingdom,
          feature_importance = .workflow %>% map(~ .x$fit$fit$fit$variable.importance %>%
            enframe() %>%
            arrange(-value) %>%
            head(5)
            )
        ) %>%
        unnest(feature_importance) %>%
        rename(feature = name, importance = value) %>%
        left_join(features_meta_tbl) %>%
        left_join(lineages) %>%
        mutate(
          phylum = phylum %>% map_chr(~ ifelse(.x %in% names(phyla_colors), .x, NA)),
          facet = outcome_kingdom %>% recode("Bacteria" = "Bacterial dysbalance", "Fungi" = "Fungal dysbalance")
        ) %>%
        ggplot(aes(predictor_kingdom, importance, group = taxon, fill = phylum)) +
        geom_col(color = "black", alpha = 0.6) +
        geom_text(aes(label = taxon), position = position_stack(vjust = 0.5), size = 2.5) +
        facet_wrap(~facet, scales = "free") +
        scale_fill_phyla() +
        scale_y_continuous(expand = c(0, 0)) +
        theme(
          legend.position = "none",
          panel.grid.major.x = element_blank()
        ) +
        labs(
          x = "Predictor kingdom",
          y = "Gini feature importance"
        )
    }
  )
}
