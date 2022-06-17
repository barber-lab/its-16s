get_graph_plan <- function() {
  drake_plan(
    graph_cor_method = "glasso",
    node_topology_metrics = c("closeness", "betweeness", "degree"),
    coabundance_max_pval = 0.05,
    coabundance_min_abs_estimate = 0.1,

    graph_params_description = paste(
      graph_cor_method,
      ", Prevalence >=",
      sprintf("%.i%%", coabundance_min_prevalence_perc),
      ", p < ",
      coabundance_max_pval,
      ", |r| >=",
      coabundance_min_abs_estimate,
      sep = ""
    ),

    cooccurrence_graphs = target(
      {
        cooccurrences %>%
          mutate(
            cor_method = "jaccard",
            cor_graph = list(cooccurence) %>% pmap(possibly(~ {
              .x %>%
                activate(nodes) %>%
                rename(taxon = name) %>%
                left_join(lineages) %>%
                annotate_node_attributes_in_edges() %>%
                as_tbl_graph(
                  cor_res = .,
                  nodes = lineages,
                  method = "jaccard",
                  max_pval = coabundance_max_pval,
                  min_abs_estimate = coabundance_min_abs_estimate
                )
            }, NA))
          ) %>%
          select(-where(is.list), cor_graph)
      },
      dynamic = map(cooccurrences),
      hpc = FALSE
    ),

    single_graphs = target(
      {
        coabundances %>%
          mutate(
            cor_graph = list(cor_res, cor_method) %>% pmap(possibly(~ {
              if ("tbl_graph" %in% class(.x)) {
                .x %>%
                  activate(edges) %>%
                  filter(
                    p.value <= coabundance_max_pval &
                      abs(estimate) >= coabundance_max_pval
                  ) %>%
                  activate(nodes) %>%
                  mutate(taxon = name) %>%
                  left_join(lineages) %>%
                  annotate_node_attributes_in_edges()
              } else {
                as_tbl_graph(
                  cor_res = .x,
                  nodes = lineages,
                  method = .y,
                  max_pval = coabundance_max_pval,
                  min_abs_estimate = coabundance_min_abs_estimate
                )
              }
            }, NA))
          ) %>%
          select(-where(is.list), cor_graph)
      },
      dynamic = map(coabundances),
      hpc = FALSE
    ),

    # ensemble_graphs = {
    #   single_graphs %>%
    #     pivot_wider(names_from = cor_method, values_from = cor_graph) %>%
    #     dplyr::filter(!is.na(sparcc) & !is.na(glasso)) %>%
    #     mutate(
    #       cor_method = "sparcc&glasso",
    #       cor_graph = sparcc %>% map2(glasso, ~ possibly(ensemble_coabundance, NA)(.x, .y, method = "intersection"))
    #     ) %>%
    #     filter(!is.na(cor_graph)) %>%
    #     select(subset_name, subset_value, norm_method, cor_method, cor_graph)
    # },

    graphs = target(
      {
        bind_rows(
          single_graphs,
          cooccurrence_graphs,
          # ensemble_graphs
        )
      },
      hpc = FALSE
    ),

    graphs_file = {
      graphs %>%
        write_rds(file_out("results/graphs.rds"), compress = "gz")
    },

    graph_summary = target(
      {
        summarise_graphs(graphs = graphs %>% filter(cor_method != "spearman"))
      },
      hpc = FALSE
    ),

    #
    # (Multi)graphs per environment and kingdom
    #

    subset_names = {
      subsets$subset_name %>% unique()
    },

    multi_graphs = target(
      {
        # clean
        graphs %>%
          filter(subset_name == subset_names & cor_method == cor_methods) %>%
          filter(!is.na(cor_graph)) %>%
          mutate(
            cor_graph = cor_graph %>% map2(subset_value, ~ {
              .x %>%
                activate(edges) %>%
                mutate(subset_value = .y) %>%
                select(from, to, estimate, subset_value) %>%

                # remove topology to unify nodes
                activate(nodes) %>%
                select(any_of("taxon"))
            })
          ) %>%
          # reduce
          group_by(subset_name, norm_method, cor_method) %>%
          summarise(cor_graph = .$cor_graph %>% reduce(graph_join) %>% list())
      },
      dynamic = cross(subset_names, cor_methods),
      hpc = TRUE
    ),

    multigraphs_by = {
      graphs %>%
        unite("group", subset_name, norm_method, cor_method) %>%
        pull(group)
    },

    multigraphs = target(
      {
        all_edges <-
          graphs %>%
          filter(!is.na(cor_graph)) %>%
          rename(data = cor_graph) %>%
          mutate(data = data %>% map(~ {
            .x %>%
              activate(edges) %>%
              as_tibble() %>%
              distinct(from_taxon, to_taxon, sign_estimate = sign(estimate), estimate)
          })) %>%
          unnest(data) %>%
          mutate(
            # ensure that cor(a,b) is identical to cor(b,a)
            from_taxon = map2_chr(from_taxon, to_taxon, ~ c(.x, .y) %>%
              sort() %>%
              first()),
            to_taxon = map2_chr(from_taxon, to_taxon, ~ c(.x, .y) %>%
              sort() %>%
              last())
          ) %>%
          group_by(from_taxon, to_taxon, sign_estimate) %>%
          arrange(from_taxon, to_taxon, sign_estimate)

        multigraphs <-
          all_edges %>%
          nest(-c(subset_name, norm_method, cor_method)) %>%
          mutate(
            data = data %>% map2(subset_name, function(data, subset_name) {
              edges <-
                data %>%
                left_join(lineages %>% rename_all(~ paste0("from_", .))) %>%
                left_join(lineages %>% rename_all(~ paste0("to_", .))) %>%
                mutate(
                  kingdom_group = case_when(
                    from_kingdom != to_kingdom ~ "inter kingdom",
                    from_kingdom == "Bacteria" & to_kingdom == "Bacteria" ~ "within bacteria",
                    from_kingdom == "Fungi" & to_kingdom == "Fungi" ~ "within fungi"
                  )
                ) %>%
                rename_at("subset_value", ~subset_name) %>%
                rename(from = from_taxon, to = to_taxon)

              nodes <-
                lineages %>%
                filter(taxon %in% edges$from | taxon %in% edges$to) %>%
                mutate(phylum_color = phyla_colors[phylum] %>% replace_na("") %>% as.character())

              tbl_graph(nodes, edges)
            })
          ) %>%
          rename(cor_graph = data)

        multigraphs
      },
      dynamic = group(graphs, .by = multigraphs_by),
      hpc = TRUE
    ),

    environment_group_multi_graph = {
      graphs %>%
        filter(
          norm_method == "raw" &
            cor_method == graph_cor_method &
            subset_name == "environment_group"
        ) %>%
        rename(environment_group = subset_value) %>%
        mutate(
          cor_graph = cor_graph %>% map2(environment_group, ~ {
            .x %>%
              activate(edges) %>%
              mutate(environment_group = .y) %>%

              # remove topology to unify nodes
              activate(nodes) %>%
              select(any_of(c("taxon", "genus", "phylum", "kingdom")))
          })
        ) %>%
        pull(cor_graph) %>%
        reduce(graph_join) %>%
        to_undirected() %>%
        tidygraph::activate(nodes) %>%
        tidygraph::left_join(lineages) %>%
        annotate_node_attributes_in_edges() %>%
        topologize_graph() %>%
        activate(edges) %>%
        arrange(from, to, environment_group)
    },

    environment_group_multi_graph_merged_plt = {
      multiedge_nodes <-
        environment_group_multi_graph %>%
        # filter multiedges
        activate(edges) %>%
        group_by(from, to) %>%
        mutate(n = n()) %>%
        ungroup() %>%
        filter(n > 1) %>%
        activate(nodes) %>%
        filter(!node_is_isolated()) %>%
        as_tibble()

      environment_group_multi_graph %>%
        activate(nodes) %>%
        mutate(
          Phylum = phylum %>% {
            .x <- .
            ifelse(.x %in% names(phyla_colors), .x, NA)
          },
          `Multi env node` = taxon %in% multiedge_nodes$taxon
        ) %>%
        activate(edges) %>%
        mutate(
          Environment = environment_group,
          `Effect size` = abs(estimate),
        ) %>%
        ggraph(layout = "igraph", algorithm = "fr", weights = abs(estimate)) +
        geom_edge_link(aes(color = Environment, alpha = `Effect size`)) +
        geom_node_point(aes(color = Phylum, shape = `Multi env node`), size = 3) +
        # geom_node_text(
        #   mapping = aes(label = ifelse(taxon %in% multiedge_nodes$taxon, taxon, NA)),
        #   nudge_y = -0.1
        # ) +
        scale_edge_color_manual(values = environment_group_colors) +
        scale_color_phyla() +
        theme_void() +
        labs(
          title = "Coabundance multigraph of environments",
          subtitle = graph_params_description
        )
    },

    environment_group_coabundance_venn_plt = {
      environment_group_multi_graph %>%
        activate(edges) %>%
        as_tibble() %>%
        mutate(edge = paste0(from, "-", to)) %>%
        group_by(environment_group) %>%
        nest() %>%
        mutate(
          edges = data %>% map(~ {
            .x %>% pull(edge)
          })
        ) %>%
        select(environment_group, edges) %>%
        # arrange to fix venn fill color codes
        ungroup() %>%
        mutate(environment_group = environment_group %>% factor(levels = environment_groups)) %>%
        deframe() %>%
        eulerr::euler(shape = "circle") %>%
        plot(fills = environment_group_colors[environment_groups], quantities = TRUE)
    },

    environment_group_graph_plt = {
      graphs %>%
        filter(subset_name == "environment_group" & cor_method == graph_cor_method) %>%
        mutate(
          cor_plt = cor_graph %>% map2(subset_value, ~ {
            .x %>%
              # component size
              activate(nodes) %>%
              mutate(component = group_components()) %>%
              group_by(component) %>%
              mutate(component_size = length(taxon)) %>%
              filter(component_size > 3) %>%
              ungroup() %>%
              morph(to_components) %>%
              keep(~ .x %>%
                activate(nodes) %>%
                as_tibble() %>%
                nrow() > 10) %>%
              # unmorph does not work
              reduce(graph_join) %>%
              activate(nodes) %>%
              mutate(phylum = ifelse(phylum %in% names(phyla_colors), phylum, "other")) %>%
              ggraph(layout = "igraph", algorithm = "nicely") +
              geom_edge_link(aes(color = estimate), size = 6) +
              geom_node_point(aes(color = phylum)) +
              scale_edge_color_gradient2(
                high = "#770000", low = "#000077", midpoint = 0,
                limits = c(-1, 1)
              ) +
              scale_color_phyla() +
              theme_void() +
              coord_fixed() +
              labs(subtitle = .y)
          })
        ) %>%
        pull(cor_plt) %>%
        wrap_plots() +
        plot_layout(guides = "collect") +
        plot_annotation(
          title = "Coabundance per environment",
          subtitle = graph_params_description
        )
    },

    partial_sub_mats = sub_mats %>% filter(subset_name == "environment_group"),

    partial_coabundances = target(
      {
        res <- tryCatch(
          error = function(x) NA,
          {
            subset <-
              sub_mats %>%
              left_join(graphs) %>%
              left_join(sub_abundances %>% rename(abundance = data)) %>%
              # pearson like estimates
              filter(norm_method == "raw" & cor_method == "sparcc") %>%
              slice(1) %>%
              as.list()

            if (subset$data[[1]] %>% length() > 1) {
              meta_tbl <-
                subset$data[[1]] %>%
                # merge all abundances together if list of matricies provided
                map(~ .x %>% as_tibble(rownames = "sample_id")) %>%
                bind_cols() %>%
                left_join(samples)
            } else {
              meta_tbl <-
                subset$data[[1]] %>%
                # merge all abundances together if list of matricies provided
                as_tibble(rownames = "sample_id") %>%
                left_join(samples)
            }

            subset$cor_graph[[1]] %>%
              activate(edges) %>%
              as_tibble() %>%
              select(from_taxon, to_taxon, estimate, p.value) %>%
              dplyr::mutate(
                pcor = list(from_taxon, to_taxon, estimate) %>%
                  pmap(~ possibly(partial_correlation, NA)(
                    from_taxon = ..1,
                    to_taxon = ..2,
                    estimate = ..3,
                    # covariates = c("lat", "lon", "collection_datetime", "habitat"),
                    covariates = "habitat",
                    data = meta_tbl
                  ))
              ) %>%
              select(-estimate) %>%
              unnest_wider(pcor)
          }
        )

        if (is.na(res)) {
          return(tibble())
        }

        sub_mats %>%
          select(subset_name, subset_value) %>%
          mutate(
            cor_method = "sparcc",
            cor_res = list(res)
          )
      },
      dynamic = map(partial_sub_mats),
      hpc = TRUE
    ),

    #
    # Analyze graphs -----
    #

    environment_group_multi_graph_edges = {
      environment_group_multi_graph %>%
        annotate_node_attributes_in_edges() %>%
        activate(edges) %>%
        as_tibble()
    },

    bioproject_specific_coabundances = {
      graphs %>%
        filter(subset_name == "environment_group,bioproject_id,kingdom" & cor_method == "sparcc") %>%
        filter(!is.na(cor_graph)) %>%
        separate(subset_value, into = c("environment_group", "bioproject_id", "kingdom"), sep = ",") %>%
        transmute(
          environment_group,
          bioproject_id,
          kingdom,
          data = cor_graph %>% map(~ .x %>%
            activate(edges) %>%
            as_tibble())
        ) %>%
        unnest(data) %>%
        inner_join(samples %>% distinct(bioproject_id, environment_group, habitat)) %>%
        group_by(habitat, from_taxon, to_taxon) %>%
        nest() %>%
        ungroup() %>%
        transmute(
          habitat, from_taxon, to_taxon,
          is_specific = data %>% map_lgl(~ .x$estimate %>% has_iqr_outlier())
        )
    },

    environment_group_multi_graph_edges_by = {
      environment_group_multi_graph_edges %>%
        transmute(batch = row_number() %% n_batches) %>%
        pull(batch)
    },

    environment_group_differential_coabundances = target(
      {
        environment_group_multi_graph_edges %>%
          rowwise() %>%
          mutate(env_differential_fdr = list(from_taxon, to_taxon, estimate) %>% pmap_dbl(~ {
            fdr_counts <-
              permuted_coabundances %>%
              mutate(cor_res = cor_res %>% map(~
              .x %>%
                possibly(annotate_name_in_edges, tbl_graph())(from_names = "name", to_names_suffix = "taxon") %>%
                activate(edges) %>%
                as_tibble())) %>%
              unnest(cor_res) %>%
              filter(from_taxon == ..1 & to_taxon == ..2) %>%
              # tow sided test
              mutate(false_discovery = abs(estimate) >= abs(..3)) %>%
              group_by(false_discovery) %>%
              count() %>%
              ungroup() %>%
              complete(false_discovery = c(TRUE, FALSE), fill = list(n = 0)) %>%
              deframe()

            case_when(
              # no false discovery found, so p value is zero
              sum(fdr_counts) == 0 ~ 0,

              # default case
              TRUE ~ fdr_counts[["TRUE"]] / sum(fdr_counts)
            )
          })) %>%
          ungroup()
      },
      dynamic = group(environment_group_multi_graph_edges, .by = environment_group_multi_graph_edges_by),
      hpc = TRUE
    ),

    # Does the coabundance strength varies across environments?
    environment_group_differential_coabundances_coefficient = {
      n_samples <-
        sub_mats %>%
        filter(subset_name == "environment_group") %>%
        transmute(
          environment_group = subset_value,
          n_samples = data %>% map_dbl(~ .x[[1]] %>% nrow())
        )

      # Correlation coeff varies between environments
      environment_group_multi_graph %>%
        activate(edges) %>%
        as_tibble() %>%

        # exclude edges unique to an environment
        group_by(from_taxon, to_taxon) %>%
        mutate(environment_group_occurance = n()) %>%
        filter(environment_group_occurance > 1) %>%
        left_join(n_samples, by = "environment_group") %>%
        mutate(
          # Assume Pearson standard error
          estimate_se = sqrt((1 - estimate^2) / (n_samples - 2))
        ) %>%
        group_by(from_taxon, to_taxon) %>%
        nest() %>%
        mutate(
          q_test = data %>% map(~ metagen(TE = .x$estimate, seTE = .x$estimate_se) %>% tidy()),
          max_group = data %>% map_chr(~ .x %>%
            arrange(-abs(estimate)) %>%
            pluck("environment_group", 1)),
          min_group = data %>% map_chr(~ .x %>%
            arrange(abs(estimate)) %>%
            pluck("environment_group", 1)),
        ) %>%
        unnest_wider(q_test) %>%
        select(-data)
    },

    environment_group_differential_coabundances_plt = {
      environment_group_differential_coabundances_counts <-
        environment_group_differential_coabundances %>%
        mutate(
          is_env_differential = env_differential_fdr <= 0.05,
          kingdom = case_when(
            from_kingdom != to_kingdom ~ "Interkingdom",
            from_kingdom == "Bacteria" & to_kingdom == "Bacteria" ~ "within Bacteria",
            from_kingdom == "Fungi" & to_kingdom == "Fungi" ~ "within Fungi"
          )
        )

      environment_group_differential_coabundances_counts %>%
        ggplot(aes(environment_group, alpha = is_env_differential, fill = kingdom)) +
        geom_bar() +
        scale_fill_manual(values = c(
          "Interkingdom" = "#5a189a",
          "within Bacteria" = kingdoms_colors[["Bacteria"]],
          "within Fungi" = kingdoms_colors[["Fungi"]]
        )) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_alpha_discrete(range = c(0.2, 1)) +
        labs(
          title = "Environment differential coabundances",
          subtitle = paste0(
            "Permutation test: SparCC on samples with shuffled environment lables, both kingdoms in abundance matrix\n",
            graph_params_description,
            ", Permutations:",
            coabundance_n_permutations
          ),
          x = "Environment",
          y = "Coabundances",
          alpha = "Environment differential",
          fill = "Kingdom"
        )
    },

    topology_plt = {
      graphs %>%
        filter(subset_name == "environment_group") %>%
        transmute(
          environment_group = subset_value,
          data = cor_graph %>% map(~ .x %>%
            activate(nodes) %>%
            as_tibble())
        ) %>%
        unnest(data) %>%
        pivot_longer(cols = node_topology_metrics) %>%
        ggplot(aes(environment_group, value, fill = environment_group)) +
        geom_boxplot() +
        stat_compare_means_environment_group() +
        scale_fill_environment_group() +
        scale_y_continuous(labels = scales::scientific) +
        facet_wrap(~name, scales = "free_y", nrow = 1) +
        guides(fill = FALSE) +
        labs(y = "Node topology", x = "")
    },

    node_topology_formulas = "degree ~ 0 + phylum + habitat + (1|bioproject_id)",

    node_topology_models = target(
      {
        data <-
          graphs %>%
          filter(subset_name == "environment_group,bioproject_id,kingdom" & cor_method == "sparcc") %>%
          distinct(subset_value, cor_graph) %>%
          separate(subset_value, sep = ",", into = c("environment_group", "bioproject_id", "kingdom")) %>%
          transmute(
            bioproject_id,
            data = cor_graph %>% map(~ .x %>%
              activate(nodes) %>%
              as_tibble())
          ) %>%
          unnest(data) %>%
          inner_join(
            samples %>% distinct(bioproject_id, environment_group, habitat)
          ) %>%
          inner_join(lineages) %>%
          filter(kingdom == kingdoms_groups) %>%

          # subset
          {
            .x <- .

            selected_phyla <- .x %>%
              count(phylum) %>%
              pull(phylum) %>%
              head(20)

            .x %>% filter(phylum %in% selected_phyla)
          }

        options(mc.cores = 20)

        res <- tryCatch(
          error = function(x) NA,
          {
            glmer.nb(as.formula(node_topology_formulas), data = data)
          }
        )

        if (is.na(res)) {
          return(tibble())
        }

        tibble(
          kingdom = kingdoms_groups,
          formula = node_topology_formulas,
          model = list(res)
        )
      },
      dynamic = cross(kingdoms_groups, node_topology_formulas),
      hpc = TRUE
    ),

    common_envs_graphs = target(
      {
        multigraphs %>%
          filter(subset_name == "environment_group") %>%
          mutate(
            cor_graph = list(cor_graph) %>% pmap(possibly(function(cor_graph) {
              # convert multigraph to one edge per taxon and sign_estimate pair
              edges <-
                cor_graph %>%
                activate(edges) %>%
                annotate_node_attributes_in_edges() %>%
                as_tibble() %>%
                select(from_taxon, to_taxon, sign_estimate, environment_group, estimate) %>%
                # ensure there is only one edge per group
                group_by(from_taxon, to_taxon, sign_estimate, environment_group) %>%
                summarise(estimate = mean(estimate)) %>%
                group_by(from_taxon, to_taxon, sign_estimate) %>%
                # edge must be found in most environments
                filter(n() >= 2) %>%
                pivot_wider(names_from = environment_group, values_from = estimate) %>%
                mutate(
                  # to decide which edge to keep in MST
                  weight = -sum(abs(aquatic), abs(host), abs(soil), na.rm = TRUE),
                  environment_group = case_when(
                    is.na(host) & !is.na(aquatic) & !is.na(soil) ~ "aquatic-soil",
                    !is.na(host) & !is.na(aquatic) & is.na(soil) ~ "host-aquatic",
                    !is.na(host) & is.na(aquatic) & !is.na(soil) ~ "host-soil",
                    !is.na(host) & !is.na(aquatic) & !is.na(soil) ~ "host-aquatic-soil"
                  )
                ) %>%
                left_join(lineages %>% rename_all(~ paste0("from_", .))) %>%
                left_join(lineages %>% rename_all(~ paste0("to_", .))) %>%
                mutate(
                  kingdom_group = case_when(
                    from_kingdom == "Fungi" & to_kingdom == "Fungi" ~ "within Fungi",
                    from_kingdom == "Bacteria" & to_kingdom == "Bacteria" ~ "within Bacteria",
                    from_kingdom != to_kingdom ~ "inter kingdom"
                  )
                ) %>%
                filter(from_genus != to_genus)

              nodes <-
                lineages %>%
                filter(taxon %in% edges$from_taxon | taxon %in% edges$to_taxon) %>%
                mutate(phylum_color = phyla_colors[phylum] %>% replace_na("") %>% as.character())

              graph <-
                tbl_graph(nodes, edges) %>%
                activate(edges) %>%
                filter(from_taxon != to_taxon) %>%
                activate(nodes) %>%
                filter(!node_is_isolated()) %>%
                annotate_node_attributes_in_edges()

              graph
            }, NA))
          )
      },
      hpc = FALSE
    ),

    intersect_inter_kingdoms_weightings = c("weighted", "unweighted"),
    intersect_inter_kingdoms_hops = c(1, 2),

    #' Interactions found in most environments
    #' Keep inter kingdom edges + strongest edges of taxa
    #' Keep taxa having maximal 2 hops to a fungus involved in inter kingdom interaction
    intersect_env_inter_kingdoms_hops_graphs = target(
      {
        common_envs_graphs %>%
          filter(subset_name == "environment_group") %>%
          expand_grid(hops = intersect_inter_kingdoms_hops, weighting = intersect_inter_kingdoms_weightings) %>%
          mutate(
            cor_graph = list(cor_graph, hops, weighting) %>% pmap(possibly(function(cor_graph, hops, weighting) {
              #
              # Filter edges: max 2 hops to inter kingdom interaction
              #

              # Fungi involved in inter kingdom interactions
              # This is faster than taking all of the bacteria
              inter_kingdom_fungi_names <-
                graph %>%
                activate(edges) %>%
                filter(kingdom_group == "inter kingdom") %>%
                activate(nodes) %>%
                mutate(name = row_number()) %>%
                filter(kingdom == "Fungi") %>%
                as_tibble() %>%
                pull(name)

              graph %>%
                as.igraph() %>%
                # get all nodes <=.y hops to an inter kingdom interaction
                make_ego_graph(nodes = inter_kingdom_fungi_names, order = hops) %>%
                # keep only shortest paths
                map(~ {
                  switch(weighting,
                    weighted = igraph::mst(.x, weights = E(.x)$weight),
                    unweighted = igraph::mst(.x, algorithm = "unweighted"),
                    nofilter = .x
                  )
                }) %>%
                # join all graphs together
                map(as_tbl_graph) %>%
                reduce(graph_join) %>%
                activate(edges) %>%
                distinct(from, to, .keep_all = TRUE) %>%
                activate(nodes) %>%
                filter(!node_is_isolated())
            }, NA))
          )
      },
      hpc = TRUE,
      dynamic = cross(intersect_inter_kingdoms_weightings, intersect_inter_kingdoms_hops)
    ),

    #
    # Analyze taxa importance ----
    #

    keystones = {
      graphs %>%
        filter(!is.na(cor_graph)) %>%
        mutate(nodes_tbl = cor_graph %>% map(~ .x %>%
          activate(nodes) %>%
          as_tibble())) %>%
        select(-cor_graph) %>%
        unnest(nodes_tbl) %>%
        group_by(samples_grouping, samples_group, taxa_grouping, taxa_group, norm_method, cor_method) %>%
        nest() %>%
        mutate(data = data %>% map(~ .x %>% select(taxon, keystoneness)))
    },

    foundations = {
      graphs %>%
        filter(!is.na(cor_graph)) %>%
        mutate(nodes_tbl = cor_graph %>% map(~ .x %>%
          activate(nodes) %>%
          as_tibble())) %>%
        select(-cor_graph) %>%
        unnest(nodes_tbl) %>%
        group_by(samples_grouping, samples_group, taxa_grouping, taxa_group, norm_method, cor_method) %>%
        nest() %>%
        mutate(data = data %>% map(~ .x %>% select(taxon, foundationess)))
    },

    abundance_topology = {
      get_abundance_topology(sub_abundances, graphs, node_topology_metrics, keystones, foundations)
    },

    abundance_topology_summary = {
      summarise_abundance_topology(abundance_topology)
    },

    get_signed_node_degrees = target(
      {
        plt <- get_signed_node_degrees(graphs$cor_graph[[1]])

        graphs %>%
          select(-cor_graph) %>%
          mutate(plt = list(plt))
      },
      dynamic = map(graphs),
      hpc = FALSE
    ),

    environment_group_inter_intra_kingdom_euler_plts = target(
      {
        environment_group_correlations <-
          graphs %>%
          filter(subset_name == "environment_group" & cor_method == "sparcc") %>%
          rename(environment_group = subset_value) %>%
          distinct(environment_group, .keep_all = TRUE) %>%
          transmute(
            environment_group,
            data = cor_graph %>% map2(environment_group, ~ {
              selected_taxa <-
                majority_prevalent_taxa %>%
                filter(environment_group == .y) %>%
                pull(taxon)

              .x %>%
                activate(nodes) %>%
                mutate(taxon = name) %>%
                annotate_node_attributes_in_edges() %>%
                activate(edges) %>%
                as_tibble() %>%
                select(from_taxon, to_taxon, p.value, estimate) %>%
                left_join(lineages %>% select(from_taxon = taxon, from_kingdom = kingdom)) %>%
                left_join(lineages %>% select(to_taxon = taxon, to_kingdom = kingdom)) %>%
                mutate(
                  kingdom_group = case_when(
                    from_kingdom == "Bacteria" & to_kingdom == "Bacteria" ~ "within bacteria",
                    from_kingdom == "Fungi" & to_kingdom == "Fungi" ~ "within fungi",
                    from_kingdom != to_kingdom ~ "inter kingdom"
                  )
                )
              # filter(from_taxon %in% selected_taxa & to_taxon %in% selected_taxa)
            })
          ) %>%
          unnest(data) %>%
          mutate(estimate_sign = sign(estimate)) %>%
          filter(p.value < 0.05)


        taxa_pairs <-
          environment_group_correlations$from_taxon %>%
          union(environment_group_correlations$to_taxon) %>%
          combn(2) %>%
          t() %>%
          as_tibble() %>%
          rename(from_taxon = V1, to_taxon = V2)

        plt <-
          environment_group_correlations %>%
          nest(-kingdom_group) %>%
          mutate(
            plt = data %>% map2(kingdom_group, possibly(~ {
              res <-
                .x %>%
                  mutate(
                    correlation = list(from_taxon, to_taxon, estimate_sign) %>% pmap_chr(~ {
                      c(..1, ..2) %>%
                        sort() %>%
                        paste(collapse = "~") %>%
                        paste0("~", ..3)
                    })
                  ) %>%
                  distinct(environment_group, correlation) %>%
                  nest(-environment_group) %>%
                  deframe() %>%
                  map(~ .x %>% pluck("correlation")) -> res

              # Note: Ensure distinct abbreviations
              colors <- environment_group_colors
              colors <- colors[names(res)]

              res %>%
                eulerr::euler(shape = "circle") %>%
                plot(quantities = list(cex = 0.7), labels = FALSE, fills = colors, alpha = hmsc_euler_alpha, main = .y)
            }, NA))
          ) %>%
          filter(!is.na(kingdom_group) & !is.na(plt)) %>%
          pull(plt) %>%
          wrap_plots(nrow = 1)

        plt
      },
      hpc = FALSE
    ),

    habitat_inter_intra_kingdom_euler_plts = {
      majority_prevalent_taxa <-
        prevalences %>%
        filter(subset_name == "habitat") %>%
        rename(habitat = subset_value) %>%
        left_join(habitats) %>%
        unnest(data) %>%
        nest(-taxon, -environment_group) %>%
        mutate(
          is_majority_prevalent = data %>% map2_lgl(environment_group, ~ {
            n_habitats <- habitats %>%
              filter(environment_group == .y) %>%
              nrow()
            nrow(.x) / n_habitats > 0.5
          })
        ) %>%
        filter(is_majority_prevalent) %>%
        select(environment_group, taxon)

      habitat_correlations <-
        graphs %>%
        filter(subset_name == "habitat" & cor_method == "sparcc") %>%
        rename(habitat = subset_value) %>%
        distinct(habitat, .keep_all = TRUE) %>%
        left_join(habitats) %>%
        transmute(
          habitat,
          data = cor_graph %>% map2(environment_group, ~ {
            selected_taxa <-
              majority_prevalent_taxa %>%
              filter(environment_group == .y) %>%
              pull(taxon)

            .x %>%
              activate(nodes) %>%
              mutate(taxon = name) %>%
              annotate_node_attributes_in_edges() %>%
              activate(edges) %>%
              as_tibble() %>%
              select(from_taxon, to_taxon, p.value, estimate) %>%
              left_join(lineages %>% select(from_taxon = taxon, from_kingdom = kingdom)) %>%
              left_join(lineages %>% select(to_taxon = taxon, to_kingdom = kingdom)) %>%
              mutate(
                kingdom_group = case_when(
                  from_kingdom == "Bacteria" & to_kingdom == "Bacteria" ~ "within bacteria",
                  from_kingdom == "Fungi" & to_kingdom == "Fungi" ~ "within fungi",
                  from_kingdom != to_kingdom ~ "inter kingdom"
                )
              )
            # filter(from_taxon %in% selected_taxa & to_taxon %in% selected_taxa)
          })
        ) %>%
        unnest(data) %>%
        mutate(estimate_sign = sign(estimate)) %>%
        left_join(habitats) %>%
        filter(p.value < 0.05)

      taxa_pairs <-
        habitat_correlations$from_taxon %>%
        union(habitat_correlations$to_taxon) %>%
        combn(2) %>%
        t() %>%
        as_tibble() %>%
        rename(from_taxon = V1, to_taxon = V2)

      plt <-
        habitat_correlations %>%
          nest(-c(environment_group, kingdom_group)) %>%
          mutate(
            plt = data %>% map(~ {
              res <-
                .x %>%
                mutate(
                  correlation = list(from_taxon, to_taxon, estimate_sign) %>% pmap_chr(~ {
                    c(..1, ..2) %>%
                      sort() %>%
                      paste(collapse = "~") %>%
                      paste0("~", ..3)
                  })
                ) %>%
                distinct(habitat, correlation) %>%
                nest(-habitat) %>%
                mutate(habitat = habitat %>% factor(levels = habitats$habitat) %>% map_chr(abbreviate)) %>%
                deframe() %>%
                map(~ .x %>% pluck("correlation"))

              # Note: Ensure distinct abbreviations
              colors <- habitat_colors
              names(colors) <- names(colors) %>% map_chr(abbreviate)
              colors <- colors[names(res)]

              res <- res %>% ggeulerr() + scale_fill_manual(values = colors)
            })
          ) %>%
          filter(!is.na(kingdom_group) & !is.na(plt)) %>%
          rename(Environment = environment_group, `Kingdom group` = kingdom_group) %>%
          wrap_plots_grid(formula = formula(Environment ~ `Kingdom group`), plot_column = "plt") &
          theme(legend.position = "none")

      (wrap_elements(plt) / hmsc_euler_legend_plt) +
        plot_layout(heights = c(5, 1)) +
        plot_annotation()
    }
  )
}
