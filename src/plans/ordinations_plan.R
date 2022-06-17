get_ordinations_plan <- function() {
  drake_plan(
    dbrda_max_n = 50,
    max_clusters = 10,
    ord_methods = "PCoA",
    ord_distances = "bray",

    univariate_constrained_ordinations = target(
      {
        tbl <-
          list(
            sub_abundances %>% rename(sub_abundance = data),
            distances,
            prevalences %>% rename(prevalences = data)
          ) %>%
          reduce(inner_join) %>%
          filter(
            dist_type == "samples" &
              norm_method == "tss" &
              subset_name == "kingdom" &
              subset_value == kingdoms_groups
          )

        if (nrow(tbl) != 1) stop("univariate_constrained_ordinations unanbigious")

        taxa <-
          tbl$prevalences[[1]] %>%
          arrange(-prevalence_perc) %>%
          pull(taxon) %>%
          head(dbrda_max_n)

        other_kingdom <- kingdoms_groups %>% {
          ifelse(. == "Fungi", "Bacteria", "Fungi")
        }

        meta <- get_samples_abundance_meta(sub_abundance = tbl$sub_abundance[[1]], other_kingdom = other_kingdom)
        dist <- tbl$dist[[1]]

        selected_samples <-
          meta %>%
          pull(sample_id) %>%
          unique() %>%
          intersect(
            dist %>% as.matrix() %>% rownames()
          )

        dist %<>% usedist::dist_subset(selected_samples)
        meta %<>% filter(sample_id %in% selected_samples)


        tibble(kingdom = kingdoms_groups) %>%
          expand_grid(taxon = taxa) %>%
          mutate(
            res = taxon %>% map(~ vegan::capscale(
              formula = str_glue("dist ~ {.x}") %>% as.formula(),
              data = meta,
              distance = dist
            ) %>% RsquareAdj())
          ) %>%
          unnest_wider(res)
      },
      dynamic = map(kingdoms_groups)
    ),


    preselected_constrained_ordination_subsets = target(
      {
        subsets <-
          constrained_ordinations_subsets %>%
          filter(subset_name %>% str_detect("kingdom")) %>%
          mutate(other_kingdom = subset_value %>% str_extract("Fungi|Bacteria") %>% recode("Fungi" = "Bacteria", "Bacteria" = "Fungi"))

        generalists_preselections <-
          selected_generalists %>%
          left_join(lineages %>% select(taxon, other_kingdom = kingdom)) %>%
          rename(taxon_group = prevalence_group) %>%
          nest(-c(other_kingdom, taxon_group)) %>%
          rename(constrain_taxa = data) %>%
          mutate(constrain_taxa = constrain_taxa %>% map(~ .x[[1]])) %>%
          inner_join(subsets) %>%
          select(!where(is.list), where(is.list)) %>%
          unite(taxon_group, taxon_group, other_kingdom)

        # use specialists of the corresponding habitat
        specialists_preselections <-
          selected_specialists %>%
          left_join(bioproject_habitats) %>%
          left_join(lineages) %>%
          mutate(
            other_kingdom = kingdom %>% recode("Bacteria" = "Fungi", "Fungi" = "Bacteria"),
            subset_name = "habitat,kingdom",
            subset_value = paste0(habitat, ",", other_kingdom),
            taxon_group = paste0("specialists_", kingdom),
            constrain_taxa = taxon
          ) %>%
          select(taxon_group, subset_name, subset_value, constrain_taxa) %>%
          nest(constrain_taxa) %>%
          rename(constrain_taxa = data) %>%
          mutate(constrain_taxa = constrain_taxa %>% map(~ .x[[1]]))

        # generalists and specialists together
        generalists_specialists_preselections <-
          generalists_preselections %>%
          bind_rows(specialists_preselections) %>%
          unnest(cols = c(constrain_taxa)) %>%
          filter(subset_name == "habitat,kingdom") %>%
          mutate(taxon_group = taxon_group %>% str_replace("generalist|specialists", "generalists_specialists")) %>%
          nest(constrain_taxa) %>%
          rename(constrain_taxa = data) %>%
          mutate(constrain_taxa = constrain_taxa %>% map(~ .x[[1]]))

        list(
          generalists_preselections,
          specialists_preselections,
          generalists_specialists_preselections
        ) %>% bind_rows()
      },
      hpc = FALSE
    ),

    # Unlike stepwise constrained ordinations
    # constraining variables were preselected insted of using ordistep
    preselected_constrained_ordinations = target(
      {
        preselected_constrained_ordination_subsets %>%
          inner_join(
            distances %>%
              filter(norm_method == "tss" & dist_method == "bray" & dist_type == "samples"),
            by = c("subset_name", "subset_value")
          ) %>%
          inner_join(
            sub_abundances %>%
              filter(norm_method == "tss") %>%
              rename(abundance = data),
            by = c("subset_name", "subset_value", "norm_method")
          ) %>%
          group_by(taxon_group, subset_name, subset_value) %>%
          transmute(
            ord = list(constrain_taxa, phy, dist, abundance, subset_value) %>% pmap(possibly(
              function(constrain_taxa, phy, dist, abundance, subset_value) {
                meta <-
                  sub_abundances$data[[2]] %>%
                  filter(sample_id %in% rownames(phy@sam_data) & taxon %in% constrain_taxa) %>%
                  # pool counts
                  group_by(sample_id, genus) %>%
                  summarise(abundance = sum(abundance)) %>%
                  select(sample_id, genus, abundance) %>%
                  pivot_wider(names_from = genus, values_from = abundance, values_fill = list(abundance = 0))

                # taxa must be abundant
                constrain_taxa <- constrain_taxa %>% intersect(colnames(meta))

                selected_samples <-
                  meta %>%
                  pull(sample_id) %>%
                  unique() %>%
                  intersect(dist %>% as.matrix() %>% rownames())

                dist %<>% usedist::dist_subset(selected_samples)
                meta %<>% filter(sample_id %in% selected_samples)

                min_constrained_dbrda <- vegan::capscale(
                  formula = dist ~ 1,
                  data = meta %>% ungroup() %>% select(-sample_id)
                )

                preselected_formula <-
                  constrain_taxa %>%
                  paste0(collapse = "` + `") %>%
                  paste0("dist ~ `", ., "`") %>%
                  as.formula()

                preselected_constrained_dbrda <- vegan::capscale(
                  formula = preselected_formula,
                  data = meta %>% ungroup() %>% select(-sample_id)
                )

                # Get contribution e.g. R2

                # Do adonis because capscale does not provide R2
                preselected_adonis <- vegan::adonis(
                  formula = preselected_formula,
                  data = meta,
                  parallel = 1 # avoid overhead
                )

                univariate_adonis <-
                  constrain_taxa %>%
                  map(~ vegan::adonis(
                    formula = str_glue("dist ~ `{.x}`") %>% as.formula(),
                    data = meta,
                    parallel = 1 # avoid overhead
                  ))

                contribution <-
                  bind_rows(
                    univariate_adonis %>%
                      map(~ .x$aov.tab %>%
                        tidy() %>%
                        head(1)) %>%
                      bind_rows() %>%
                      mutate(type = "univariate"),

                    preselected_adonis$aov.tab %>% tidy() %>%
                      mutate(type = "multivariate") %>%
                      filter(!term %in% c("Residuals", "Total"))
                  )

                list(
                  min_constrained_dbrda = min_constrained_dbrda,
                  preselected_constrained_dbrda = preselected_constrained_dbrda,
                  contribution = contribution
                )
              }, NA
            ))
          )
      },
      # avoid overhead
      dynamic = map(preselected_constrained_ordination_subsets)
    ),

    preselected_constrained_ordinations_plt = target(
      {
        preselected_constrained_ordinations %>%
          filter(taxon_group %>% str_detect("generalists_specialists")) %>%
          filter(subset_name == "habitat,kingdom") %>%
          separate(subset_value, into = c("habitat", "kingdom"), sep = ",") %>%
          unnest_wider(ord) %>%
          unnest(contribution) %>%
          mutate(term = term %>% str_remove_all("[`]")) %>%
          left_join(selected_generalists_specialists, by = c("term" = "taxon")) %>%
          left_join(habitats) %>%
          mutate(
            term_color = prevalence_group %>% map_chr(~ prevalence_group_colors[.x] %>% replace_na("black")),
            taxon = str_glue("<i style='color:{term_color}'>{term}</i>"),
            R2 = R2 * 100 # percentages
          ) %>%
          group_by(term) %>%
          # both models must be present
          filter("univariate" %in% type & "multivariate" %in% type) %>%
          group_by(kingdom) %>%
          arrange(prevalence_group) %>%
          mutate(taxon = taxon %>% fct_inorder()) %>%
          ggplot(aes(taxon, R2)) +
          geom_linerange(aes(ymin = 0, ymax = R2), color = "lightgrey") +
          geom_point(aes(color = environment_group), stat = "identity", position = "dodge", size = 2) +
          scale_color_environment_group() +
          facet_grid(kingdom ~ type, scales = "free_y", space = "free") +
          scale_y_continuous(expand = c(0, 0)) +
          coord_flip() +
          theme(
            axis.text.y = element_markdown(),
            panel.grid.major.y = element_blank()
          ) +
          labs(
            title = "Preselected constrained ordinations",
            color = "Biome of the habitat",
            x = "Explanatory genus",
            y = "partial R² (%)"
          )
      },
      hpc = FALSE
    ),

    # Ordination of these kingdom constrained by the other kingdom
    # values must be available subset_values and must contain a kingdom
    constrained_ordinations_subsets = target({
      subsets %>%
        filter(subset_value %>% str_detect(kingdoms_groups %>% paste0(collapse = "|"))) %>%
        select(subset_name, subset_value) %>%
        filter(subset_name %in% c("kingdom", "environment_group,kingdom", "habitat,kingdom"))
    }),

    stepwise_constrained_ordinations = target(
      {
        kingdom <- constrained_ordinations_subsets$subset_value[[1]] %>% str_extract(kingdoms_groups %>% paste0(collapse = "|"))
        other_kingdom <- switch(kingdom, Bacteria = "Fungi", Fungi = "Bacteria")

        subset_values <-
          constrained_ordinations_subsets$subset_value[[1]] %>%
          c(
            constrained_ordinations_subsets$subset_value[[1]] %>% str_replace(kingdom, other_kingdom)
          )

        data <-
          tibble(
            subset_value = subset_values
          ) %>%
          inner_join(
            distances %>%
              filter(norm_method == "tss" & dist_method == "bray" & dist_type == "samples"),
            by = c("subset_value")
          ) %>%
          inner_join(
            prevalences %>%
              rename(prevalence = data),
            by = c("subset_name", "subset_value")
          ) %>%
          inner_join(
            sub_abundances %>%
              filter(norm_method == "tss") %>%
              rename(abundance = data),
            by = c("subset_name", "subset_value", "norm_method")
          )

        if (nrow(data) != 2) {
          warning("Data is incomplete")
          return(tibble())
        }

        start_vars <-
          data %>%
          filter(subset_value %>% str_detect(other_kingdom)) %>%
          pull(prevalence) %>%
          first() %>%
          arrange(-prevalence_perc) %>%
          pull(taxon) %>%
          unique() %>%
          head(dbrda_max_n)

        dist <-
          data %>%
          filter(
            norm_method == "tss" & dist_method == "bray" & dist_type == "samples" &
              subset_value %>% str_detect(kingdom)
          ) %>%
          pull(dist) %>%
          first()

        meta <-
          data %>%
          filter(subset_value %>% str_detect(other_kingdom) & norm_method == "tss") %>%
          pull(abundance) %>%
          first() %>%
          # pool counts
          group_by(sample_id, genus) %>%
          summarise(abundance = sum(abundance)) %>%
          select(sample_id, genus, abundance) %>%
          pivot_wider(names_from = genus, values_from = abundance, values_fill = list(abundance = 0))

        selected_samples <-
          meta %>%
          pull(sample_id) %>%
          unique() %>%
          intersect(
            dist %>% as.matrix() %>% rownames()
          )

        dist %<>% usedist::dist_subset(selected_samples)
        meta %<>% filter(sample_id %in% selected_samples)

        min_constrained_dbrda <- vegan::capscale(
          formula = dist ~ 1,
          data = meta %>% ungroup() %>% select(-sample_id)
        )

        max_constrained_dbrda <- vegan::capscale(
          formula = start_vars %>% paste0(collapse = "` + `") %>% paste0("dist ~ `", ., "`") %>% as.formula(),
          data = meta %>% ungroup() %>% select(-sample_id)
        )

        RhpcBLASctl::blas_set_num_threads(1)
        selected_dbrda <- vegan::ordistep(
          object = min_constrained_dbrda,
          scope = max_constrained_dbrda %>% as.formula(),
          direction = "forward",
          permutations = 100,
          steps = dbrda_max_n
        )

        res <- list(
          min_constrained_dbrda = min_constrained_dbrda,
          max_constrained_dbrda = max_constrained_dbrda,
          selected_dbrda = selected_dbrda
        )

        constrained_ordinations_subsets %>%
          mutate(
            kingdom = kingdom,
            ord = list(res)
          )
      },
      dynamic = map(constrained_ordinations_subsets),
      hpc = TRUE
    ),

    stepwise_constrained_ordinations_contribution_table = {
      stepwise_constrained_ordinations_contribution %>%
        select(-kingdom) %>%
        unnest(contribution) %>%
        select(-df, -kingdom) %>%
        arrange(subset_name, subset_value, type, -R2) %>%
        select(subset_name, subset_value, type, everything()) %>%
        left_join(selected_generalists_specialists)
    },

    stepwise_constrained_ordinations_contribution = target(
      {
        stepwise_constrained_ordinations %>%
          mutate(
            contribution = list(ord, kingdom, subset_value) %>% pmap(possibly(function(ord, kingdom, this_subset_value) {
              other_kingdom <- kingdom %>% recode("Bacteria" = "Fungi", "Fungi" = "Bacteria")
              other_subset_value <- this_subset_value %>% str_replace(kingdom, other_kingdom)

              dist <-
                distances %>%
                filter(
                  norm_method == "tss" & dist_method == "bray" &
                    dist_type == "samples" & subset_value == this_subset_value
                ) %>%
                pull(dist) %>%
                first()

              meta <-
                sub_abundances %>%
                filter(subset_value == other_subset_value & norm_method == "tss") %>%
                pull(data) %>%
                first() %>%
                # pool counts
                group_by(sample_id, genus) %>%
                summarise(abundance = sum(abundance)) %>%
                select(sample_id, genus, abundance) %>%
                pivot_wider(names_from = genus, values_from = abundance, values_fill = list(abundance = 0))

              selected_samples <-
                meta %>%
                pull(sample_id) %>%
                unique() %>%
                intersect(dist %>% as.matrix() %>% rownames())

              dist %<>% usedist::dist_subset(selected_samples)
              meta %<>% filter(sample_id %in% selected_samples)

              step_selected_formula <- ord$selected_dbrda %>% as.formula()

              # Do adonis because capscale does not provide R2
              step_selected_adonis <- vegan::adonis(
                formula = step_selected_formula,
                data = meta,
                parallel = 1 # avoid overhead
              )

              univariate_adonis <-
                step_selected_formula %>%
                attr("term.labels") %>%
                map(~ vegan::adonis(
                  formula = str_glue("dist ~ {.x}") %>% as.formula(),
                  data = meta,
                  parallel = 1 # avoid overhead
                ))

              bind_rows(
                univariate_adonis %>%
                  map(~ .x$aov.tab %>%
                    tidy() %>%
                    head(1)) %>%
                  bind_rows() %>%
                  mutate(type = "univariate"),

                step_selected_adonis$aov.tab %>% tidy() %>%
                  mutate(type = "multivariate") %>%
                  filter(!term %in% c("Residuals", "Total"))
              ) %>%
                rename(taxon = term) %>%
                mutate(kingdom = kingdom) %>%
                group_by(type) %>%
                mutate(q.value = p.value %>% p.adjust(method = "fdr"))
            }, NA))
          ) %>%
          select(-where(is.list), contribution)
      },
      dynamic = map(stepwise_constrained_ordinations),
      hpc = TRUE
    ),

    stepwise_constrained_ordinations_anova = target(
      {
        stepwise_constrained_ordinations %>%
          mutate(
            # avoid overhead
            anova_all = ord %>% map(~ .x$selected_dbrda %>% anova(permutations = 999, parallel = 1)),
            anova_terms = ord %>% map(~ .x$selected_dbrda %>% anova(by = "terms", permutations = 999, parallel = 1))
          ) %>%
          select(-ord)
      },
      hpc = FALSE # do not copy big target to all workers
    ),

    stepwise_constrained_ordinations_plots = target(
      {
        stepwise_constrained_ordinations %>%
          inner_join(stepwise_constrained_ordinations_anova) %>%
          mutate(
            plt = list(kingdom, ord, anova_all, anova_terms) %>% pmap(~ plot_dbrda(
              kingdom = ..1,
              samples = samples,
              lineages = lineages,
              constrained_ordination = ..2,
              anova_all = ..3,
              anova_terms = ..4,
              taxrank = "genus"
            ))
          ) %>%
          select(-c("ord", "anova_all", "anova_terms"))
      },
      trigger = trigger(change = c(habitat_colors, plot_dbrda)),
      hpc = FALSE
    ),

    stepwise_constrained_plot = target(
      {
        selected_stepwise_constrained_ordinations_contribution <-
          stepwise_constrained_ordinations_contribution %>%
          filter(subset_name == "kingdom") %>%
          select(-kingdom) %>%
          unnest(contribution)

        contribution_taxa_shown <-
          selected_stepwise_constrained_ordinations_contribution %>%
          filter(type == "multivariate") %>%
          group_by(kingdom) %>%
          arrange(-R2) %>%
          slice(1:10) %>%
          pull(taxon)

        #  stepwise_constrained_ordinations_contribution_plt
        tbl <-
          selected_stepwise_constrained_ordinations_contribution %>%
          filter(taxon %in% contribution_taxa_shown) %>%
          mutate(
            R2 = R2 * 100,
            label = p.value %>% map_chr(significance_label),
            taxon = {
              taxon %>%
                recode("`Burkholderia-Caballeronia-Paraburkholderia`" = "Burkholderia-Caballeronia") %>%
                str_remove_all("`")
            },
            taxon_short = taxon %>%
              str_remove_all("`") %>%
              map_chr(trim_long_names) %>%
              map_chr(~ .x %>% abbreviate(2)),
            kingdom = subset_value
          )

        stepwise_constrained_ordinations_contribution_plt <-
          tbl %>%
          left_join(selected_generalists_specialists) %>%
          mutate(
            taxon_color = prevalence_group %>% map_chr(possibly(~ prevalence_group_colors[[.x]], "black")),
            taxon = glue("<i style='color:{taxon_color}'>{taxon}</i>")
          ) %>%
          ggplot(aes(reorder(taxon, R2), R2, fill = type)) +
          geom_bar(stat = "identity", position = "dodge") +
          geom_text(aes(label = label), position = position_dodge(width = 1), hjust = -0.2, vjust = 0.5, size = 2.5) +
          coord_flip() +
          facet_wrap(~kingdom, nrow = 1, scales = "free") +
          scale_fill_modeltype() +
          scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
          theme_pub() +
          theme(axis.text.y = ggtext:::element_markdown()) +
          labs(
            tag = "B",
            x = "Explaining genus",
            y = "partial R² (%)",
            fill = "Model type"
          )

        wrap_plots(
          wrap_plots(
            stepwise_constrained_ordinations_plots %>%
              filter(kingdom == "Bacteria") %>%
              pull(plt) %>%
              first() +
              labs(tag = "A", color = "Explaining phylum") +
              guides() +
              theme_pub(),
            stepwise_constrained_ordinations_plots %>%
              filter(kingdom == "Fungi") %>%
              pull(plt) %>%
              first() +
              labs(tag = "", color = "Explaining phylum") +
              theme_pub()
          ),
          stepwise_constrained_ordinations_contribution_plt,
          ncol = 1,
          guides = "collect"
        ) &
          theme(legend.position = "right", legend.direction = "vertical")
      },
      hpc = FALSE
    ),

    stepwise_constrained_environment_group_plot = target(
      {
        data <-
          stepwise_constrained_ordinations_plots %>%
          rename(ord_plt = plt) %>%
          inner_join(stepwise_constrained_ordinations_contribution) %>%
          filter(subset_name == "environment_group,kingdom") %>%
          arrange(kingdom, subset_value) %>%
          mutate(
            taxa_shown = contribution %>% map(~ {
              .x %>%
                filter(type == "multivariate") %>%
                group_by(kingdom) %>%
                arrange(-R2) %>%
                slice(1:10) %>%
                pull(taxon) %>%
                trim_long_names()
            })
          )

        ord_plts <-
          data %>%
            mutate(
              ord_plt = ord_plt %>% map2(subset_value, ~ {
                .x +
                  facet_wrap(~ .y %>% str_replace_all(",", " ")) +
                  theme_pub() +
                  # remove environment color
                  guides(color_new = FALSE)
              })
            ) %>%
            pull(ord_plt) %>%
            wrap_plots(guides = "collect") &
            theme(legend.position = "right")

        contrib_plt <-
          data %>%
          select(subset_value, contribution) %>%
          separate(subset_value, into = c("environment_group", "kingdom")) %>%
          select(-kingdom) %>%
          unnest(contribution) %>%
          group_by(environment_group, kingdom) %>%
          mutate(
            R2 = R2 * 100,
            label = p.value %>% map_chr(significance_label),
            taxon = taxon %>% str_remove_all("`") %>% recode(`Burkholderia-Caballeronia-Paraburkholderia` = "Burkholderia-C-P")
          ) %>%
          {
            # Only keep tax among highest multivariate R2 per subset
            x <- .
            top_taxa <-
              x %>%
              distinct(taxon, type, R2) %>%
              pivot_wider(names_from = type, values_from = R2) %>%
              arrange(-multivariate) %>%
              slice(1:5) %>%
              pull(taxon) %>%
              unique()

            x %>% filter(taxon %in% top_taxa)
          } %>%
          left_join(selected_generalists_specialists) %>%

          # sort by explained variance
          group_by(kingdom) %>%
          arrange(R2) %>%
          mutate(taxon = taxon %>% fct_inorder()) %>%
          nest() %>%
          transmute(
            plt = data %>% map2(kingdom, ~ {
              .x %>%
                mutate(
                  taxon_color = prevalence_group %>% map_chr(possibly(~ prevalence_group_colors[[.x]], "black")),
                  taxon = glue("<i style='color:{taxon_color}'>{taxon}</i>")
                ) %>%
                ggplot(aes(taxon, R2, fill = type)) +
                geom_bar(stat = "identity", position = "dodge") +
                geom_text(aes(label = label), position = position_dodge(width = 1), hjust = -0.2, vjust = .5, size = 2) +
                scale_color_prevalence_group(drop = TRUE) +
                coord_flip(clip = "off") +
                guides(fill = guide_legend(reverse = TRUE)) +
                facet_wrap(~environment_group, ncol = 1, scales = "free_y", strip.position = "right") +
                scale_fill_modeltype() +
                scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
                theme_pub() +
                theme(
                  axis.text.y = ggtext::element_markdown(),
                  legend.direction = "vertical"
                ) +
                labs(
                  x = "Explanatory genus",
                  y = "partial R² (%)",
                  title = .y,
                  fill = "Model type",
                  color = "Prevalence group"
                )
            })
          ) %>%
          pull(plt) %>%
          wrap_plots(nrow = 1) +
          plot_layout(guides = "collect")

        list(
          ord_plts,
          contrib_plt
        ) %>%
          wrap_plots(ncol = 1) +
          plot_annotation(tag_levels = "A")
      },
      trigger = trigger(change = habitat_colors),
      hpc = FALSE
    ),

    ordinations = target(
      {
        sub_phys %>%
          filter(norm_method == "tss") %>%
          select(-norm_method) %>%
          expand_grid(
            ord_method = ord_methods,
            dist_method = ord_distances,
            dist_type = dist_types
          ) %>%
          mutate(
            ord = list(phy, ord_method, dist_method, dist_type) %>% pmap(possibly(~ {
              # disable internal parallelizaion to make it work with drake
              options(mc.cores = 1)

              switch(..2,
                "DCA" = phyloseq::ordinate(
                  physeq = ..1,
                  distance = "bray",
                  method = ..2
                ),
                phyloseq::ordinate(
                  physeq = ..1,
                  distance = ..3,
                  method = ..2,
                  type = ..4
                )
              )
            }, NA))
          ) %>%
          select(-where(is.list), ord)
      },
      dynamic = cross(sub_phys, ord_methods, ord_distances, dist_types),
      hpc = TRUE
    ),

    kingdom_env_habitat_ordination_plt = target(
      {
        adonis_tbl <-
          distances %>%
          filter(subset_name == "kingdom" & dist_type == "samples" & norm_method == "tss") %>%
          head(2) %>%
          expand_grid(term = c("habitat", "environment_group")) %>%
          mutate(
            adonis = dist %>% map2(term, function(dist, term) {
              data <-
                dist %>%
                as.matrix() %>%
                rownames() %>%
                tibble(sample_id = .) %>%
                left_join(samples)

              adonis_tbl <-
                adonis(
                  formula = str_glue("dist ~ {term}") %>% as.formula(),
                  strata = data$bioproject_id,
                  data = data
                ) %>%
                pluck("aov.tab") %>%
                as_tibble()

              list(
                p.value = adonis_tbl[[6]][[1]],
                r2 = adonis_tbl[[5]][[1]]
              )
            })
          ) %>%
          unnest_wider(adonis)

        plts <-
          ordinations %>%
          # filter(! is.na(ord)) %>%
          # left_join(sub_phys) %>%
          left_join(adonis_tbl) %>%
          filter(subset_name == "kingdom" & dist_type == "samples") %>%
          rename(kingdom = subset_value) %>%
          group_by(kingdom) %>%
          mutate(
            ord_tbl = ord %>% map2(phy, ~ {
              plot_ordination(
                .y,
                .x,
                justDF = TRUE
              ) %>%
                as_tibble()
            }),
            r2_axes = ord %>% map2(phy, ~ {
              plt <- plot_ordination(.y, .x)
              list(
                x = plt %>% pluck("labels", "x") %>% str_extract("[0-9.]+%"),
                y = plt %>% pluck("labels", "y") %>% str_extract("[0-9.]+%")
              )
            })
          ) %>%
          mutate(
            plt = list(ord_tbl, term, r2_axes, r2, p.value, kingdom) %>% pmap(function(ord_tbl, term, r2_axes, r2, p.value, kingdom) {
              err_tbl <-
                ord_tbl %>%
                group_by_at(term) %>%
                summarise(
                  x = mean(Axis.1),
                  y = mean(Axis.2),
                  xmin = x - sqrt(var(Axis.1)),
                  xmax = x + sqrt(var(Axis.1)),
                  ymin = y - sqrt(var(Axis.2)),
                  ymax = y + sqrt(var(Axis.2))
                )

              plt <- switch(term,
                "habitat" = {
                  ord_tbl %>%
                    mutate(
                      habitat = habitat %>% factor(levels = names(habitat_colors)),
                      facet_x = "",
                      facet_y = kingdom
                    ) %>%
                    ggplot(aes(Axis.1, Axis.2)) +
                    geom_point(color = "grey", alpha = 0.2) +
                    geom_errorbarh(
                      data = err_tbl,
                      mapping = aes(x = x, y = y, xmin = xmin, xmax = xmax, color = habitat)
                    ) +
                    geom_errorbar(
                      data = err_tbl,
                      mapping = aes(x = x, y = y, ymin = ymin, ymax = ymax, color = habitat)
                    ) +
                    scale_color_habitat() +
                    coord_fixed() +
                    annotate("text",
                      x = Inf, y = Inf, vjust = "inward", hjust = "inward",
                      label = str_glue("p={p.value}, R²={r2 %>% sprintf(fmt = '%.2f')}")
                    ) +
                    facet_grid(facet_x ~ facet_y) +
                    labs(
                      x = str_glue("PCo1 ({r2_axes$x})"),
                      y = str_glue("PCo1 ({r2_axes$y})"),
                      color = "Habitat"
                    )
                },
                "environment_group" = {
                  ord_tbl %>%
                    mutate(
                      facet_x = "",
                      facet_y = kingdom
                    ) %>%
                    ggplot(aes(Axis.1, Axis.2)) +
                    geom_point(color = "grey", alpha = 0.2) +
                    geom_errorbarh(
                      data = err_tbl,
                      mapping = aes(x = x, y = y, xmin = xmin, xmax = xmax, color = environment_group)
                    ) +
                    geom_errorbar(
                      data = err_tbl,
                      mapping = aes(x = x, y = y, ymin = ymin, ymax = ymax, color = environment_group)
                    ) +
                    scale_color_environment_group() +
                    coord_fixed() +
                    facet_grid(facet_x ~ facet_y) +
                    annotate("text",
                      x = Inf, y = Inf, vjust = "inward", hjust = "inward",
                      label = str_glue("p={p.value}, R²={r2 %>% sprintf(fmt = '%.2f')}")
                    ) +
                    labs(
                      x = str_glue("PCo1 ({r2_axes$x})"),
                      y = str_glue("PCo1 ({r2_axes$y})"),
                      color = "Environment"
                    )
                }
              )
            })
          ) %>%
          rename(Kingdom = kingdom, Group = term) %>%
          mutate(Group = Group %>% recode("habitat" = "Habitat", "environment_group" = "Environment"))

        plts %>%
          arrange(Group, Kingdom) %>%
          pull(plt) %>%
          wrap_plots(ncol = 2, guides = "collect")
      },
      # trigger on non change of non-targets
      trigger = trigger(change = habitat_colors, condition = TRUE)
    ),

    phate_data = {
      abundances %>%
        filter(norm_method == "tss") %>%
        pull(data) %>%
        first() %>%
        mutate(abundance = abundance %>% sqrt()) %>%
        pivot_wider(names_from = taxon, values_from = abundance, values_fill = list(abundance = 0)) %>%
        column_to_rownames("sample_id") %>%
        as.matrix() %>%
        phateR::phate() %>%
        pluck("embedding") %>%
        as_tibble(rownames = "sample_id") %>%
        inner_join(
          samples %>% select(sample_id, habitat, environment_group) %>% mutate(habitat = habitat %>% factor(levels = habitats$habitat %>% setdiff("unknown")))
        ) %>%
        inner_join(
          dysbalances %>% pivot_wider(names_from = kingdom, values_from = c(distance_to_medoid, cluster, rank, norm_rank, dysbalance))
        )
    },

    phate_plt = {
      list(
        phate_data %>%
          ggplot(aes(PHATE1, PHATE2)) +
          geom_point(aes(color = environment_group)) +
          scale_color_environment_group() +
          facet_wrap(~"Environment") +
          labs(color = "Environment", x = ""),

        phate_data %>%
          ggplot(aes(PHATE1, PHATE2)) +
          geom_point(aes(color = habitat)) +
          scale_color_habitat() +
          facet_wrap(~"Habitat") +
          guides(col = guide_legend(ncol = 4)) +
          labs(color = "Habitat", x = "", y = ""),

        phate_data %>%
          ggplot(aes(PHATE1, PHATE2)) +
          geom_point(aes(color = distance_to_medoid_Bacteria)) +
          scale_color_viridis_c(na.value = "transparent") +
          facet_wrap(~"Bacterial dysbalance") +
          labs(color = "Dysbalance (Bray-Curtis)"),

        phate_data %>%
          ggplot(aes(PHATE1, PHATE2)) +
          geom_point(aes(color = distance_to_medoid_Fungi)) +
          scale_color_viridis_c(na.value = "transparent") +
          facet_wrap(~"Fungal dysbalance") +
          labs(color = "Dysbalance (Bray-Curtis)", y = "")
      ) %>%
        wrap_plots(guides = "collect", ncol = 2)
    },

    environemnt_groups_canonical_correlations = {
      sub_phys %>%
        filter(subset_name == "environment_group,kingdom" & norm_method == "tss") %>%
        separate(subset_value, into = c("environment_group", "kingdom"), sep = ",") %>%
        transmute(
          environment_group, kingdom,
          otu_data = phy %>% map(~ {
            .x %>%
              otu_table() %>%
              as.data.frame()
          })
        ) %>%
        pivot_wider(names_from = kingdom, values_from = otu_data) %>%
        transmute(
          environment_group,
          cancor = Bacteria %>% map2(Fungi, possibly(~ {
            blocks <-
              .x %>%
              rownames() %>%
              tibble(sample_id = .) %>%
              left_join(samples %>% select(sample_id, bioproject_id)) %>%
              pull(bioproject_id) %>%
              factor()

            vegan::CCorA(
              Y = .x,
              X = .y,
              permutations = permute::how(blocks = blocks),
              stand.Y = TRUE,
              stand.X = TRUE
            )
          }, NA))
        ) %>%
        transmute(
          environment_group,
          mean_cor = cancor %>% map_dbl(~ .x$CanCorr %>% mean()),
          p.value = cancor %>% map_dbl(~ .x$p.perm),
          p.value.pillar = cancor %>% map_dbl(~ .x$p.Pillai),
          R2 = cancor %>% map(~ .x$RDA.adj.Rsq %>% set_names(c("pred_Fungi_from_Bacteria", "pred_Bacteria_from_Fungi")))
        ) %>%
        unnest_wider(R2)
    },

    within_habitat_dissimilarity_plt = {
      distances %>%
        filter(dist_method == "bray" & subset_name == "habitat,kingdom" & dist_type == "samples") %>%
        transmute(
          subset_value,
          dist = dist %>% map(~ .x %>% as.numeric())
        ) %>%
        unnest(dist) %>%
        separate(subset_value, into = c("habitat", "kingdom"), sep = ",") %>%
        mutate(habitat = habitat %>% factor() %>% fct_reorder(dist)) %>%
        ggplot(aes(habitat, dist, color = habitat)) +
        facet_wrap(~kingdom) +
        geom_boxplot() +
        scale_color_habitat() +
        coord_flip() +
        guides(color = FALSE) +
        labs(y = "Within habitat dissimilarity (Bray Curtis)")
    }
  )
}
