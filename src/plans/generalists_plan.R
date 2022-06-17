get_generalists_plan <- function() {
  drake_plan(
    generalists_min_prevalence_abundance = 0.01e-2,

    bioproject_prevalences = target(
      sub_abundances %>%
        filter(subset_name == "bioproject_id" & norm_method == "tss") %>%
        select(bioproject_id = subset_value, data) %>%
        mutate(
          data = data %>% map(~ {
            n_samples <- .x$sample_id %>%
              unique() %>%
              length()

            .x %>%
              # ensure only one abudnance entry per sample and taxon
              group_by(taxon, sample_id) %>%
              slice(1) %>%
              group_by(taxon) %>%
              filter(abundance >= generalists_min_prevalence_abundance) %>%
              summarise(
                prevalence_n = n(),
                prevalence_perc = prevalence_n / n_samples * 100
              )
          })
        ) %>%
        unnest(data),
      hpc = FALSE
    ),

    generalists = target(
      expand_grid(
        min_prevalence_perc = c(10, 20, 30, 40, 50, 60, 70),
        min_n_bioprojects = c(30, 25, 20, 15)
        # min_prevalence_perc = seq(0, 100, length.out = 15),
        # min_n_bioprojects = seq(0, 30, length.out = 15)
      ) %>%
        mutate(
          data = list(min_prevalence_perc, min_n_bioprojects) %>% pmap(~ {
            # high prevalence in most bioprojects
            bioproject_prevalences %>%
              group_by(taxon) %>%
              mutate(is_prevalent = prevalence_perc > .x) %>%
              filter(is_prevalent) %>%
              count(taxon, name = "n_bioprojects") %>%
              transmute(taxon, is_generalist = n_bioprojects > .y)
          })
        ) %>%
        unnest(data) %>%
        complete(min_prevalence_perc, min_n_bioprojects, taxon, fill = list(is_generalist = FALSE)),
      hpc = FALSE
    ),

    generalists_summary = target(
      generalists %>%
        filter(is_generalist) %>%
        inner_join(lineages) %>%
        count(min_prevalence_perc, min_n_bioprojects, kingdom, name = "n_general_taxa") %>%
        complete(min_prevalence_perc, min_n_bioprojects, kingdom, fill = list(n_general_taxa = 0)),
      hpc = FALSE
    ),

    bioproject_environment_groups = {
      samples %>%
        count(bioproject_id, environment_group) %>%
        group_by(bioproject_id) %>%
        arrange(-n) %>%
        slice(1) %>%
        select(environment_group) %>%
        ungroup()
    },

    bioproject_habitats = {
      samples %>%
        count(bioproject_id, habitat) %>%
        group_by(bioproject_id) %>%
        arrange(-n) %>%
        slice(1) %>%
        select(habitat) %>%
        ungroup()
    },

    # high prevalence in at least one project of all environments
    common_generalists = target(
      {
        tibble(min_prevalence_perc = seq(0, 100, by = 5)) %>%
          mutate(
            taxa = min_prevalence_perc %>% map(~ {
              bioproject_prevalences %>%
                group_by(taxon) %>%
                filter(prevalence_perc > .x) %>%
                inner_join(bioproject_environment_groups, by = "bioproject_id") %>%
                count(taxon, environment_group, name = "n_bioprojects") %>%
                group_by(taxon) %>%
                mutate(is_generalist = n_bioprojects > 0 & n() == 3) %>%
                pivot_wider(
                  names_from = environment_group, values_from = n_bioprojects,
                  values_fill = list(n_bioprojects = 0), names_prefix = "n_bioprojects_"
                )
            })
          ) %>%
          unnest(taxa) %>%
          select(-n_bioprojects)
      },
      hpc = FALSE
    ),

    specialists = target(
      expand_grid(
        bioproject_id = bioproject_groups,
        min_prevalence_perc_this_bioproject = c(10, 20, 30, 40, 50, 60, 70, 80),
        max_prevalence_perc_other_bioprojects = c(1, 5, 10)
        # min_prevalence_perc_this_bioproject = seq(0, 100, length.out = 13),
        # max_prevalence_perc_other_bioprojects = seq(0, 25, length.out = 13),
      ) %>%
        mutate(
          data = list(bioproject_id, min_prevalence_perc_this_bioproject, max_prevalence_perc_other_bioprojects) %>% pmap(~ {
            # high prevalence in the current bioproject
            current_bioproject_data <-
              bioproject_prevalences %>%
              filter(bioproject_id == .x) %>%
              transmute(
                bioproject_id, taxon,
                is_specialist_current = prevalence_perc > .y
              )

            # low prevalence in all other projects
            other_bioprojects_data <-
              bioproject_prevalences %>%
              filter(bioproject_id != .x) %>%
              transmute(
                bioproject_id, taxon,
                is_specialist_other = prevalence_perc < ..3
              ) %>%
              group_by(taxon) %>%
              summarise(is_specialist_other = all(is_specialist_other))

            current_bioproject_data %>%
              inner_join(other_bioprojects_data, by = "taxon") %>%
              transmute(
                bioproject_id, taxon,
                is_specialist = is_specialist_current & is_specialist_other
              )
          })
        ) %>%
        select(-bioproject_id) %>%
        unnest(data),
      hpc = FALSE
    ),

    specialists_summary = target(
      specialists %>%
        group_by(min_prevalence_perc_this_bioproject, max_prevalence_perc_other_bioprojects, taxon) %>%
        summarise(is_specialist = any(is_specialist)) %>%
        left_join(lineages) %>%
        filter(is_specialist) %>%
        count(min_prevalence_perc_this_bioproject, max_prevalence_perc_other_bioprojects, kingdom),
      hpc = FALSE
    ),

    generalists_and_specialists = target(
      common_generalists %>%
        select(-matches("n_bioprojects_")) %>%
        full_join(specialists) %>%
        filter(is_generalist | is_specialist) %>%
        # relevant parameter space
        # discard trivial cases
        filter(
          min_prevalence_perc %in% seq(40, 60, by = 10) &
            min_prevalence_perc_this_bioproject %in% seq(40, 60, by = 10) &
            max_prevalence_perc_other_bioprojects == 5
        ),
      hpc = FALSE
    ),

    taxa_prevalence_summary = target(
      generalists_and_specialists %>%
        group_by(taxon) %>%
        summarise(
          prevalence_group = case_when(
            is_generalist & is_specialist ~ "both",
            is_generalist ~ "generalist",
            is_specialist ~ "specialist"
          ),
          prevalence_strictness = case_when(
            prevalence_group == "generalist" ~ min_prevalence_perc,
            prevalence_group == "specialist" ~ min_prevalence_perc_this_bioproject
          ) %>% max()
        ) %>%
        slice(1) %>%
        ungroup(),
      hpc = FALSE
    ),

    selected_generalists_specialists = target(
      bind_rows(
        common_generalists %>%
          filter(is_generalist & min_prevalence_perc == 40) %>%
          distinct(taxon) %>%
          mutate(prevalence_group = "generalist"),
        specialists %>%
          filter(is_specialist & min_prevalence_perc_this_bioproject == 40 & max_prevalence_perc_other_bioprojects == 5) %>%
          distinct(taxon) %>%
          mutate(prevalence_group = "specialist")
      ),
      hpc = FALSE
    ),

    selected_generalists = {
      common_generalists %>%
        filter(is_generalist & min_prevalence_perc == 40) %>%
        distinct(taxon) %>%
        mutate(prevalence_group = "generalist")
    },

    selected_specialists = {
      specialists %>%
        filter(is_specialist & min_prevalence_perc_this_bioproject == 40 & max_prevalence_perc_other_bioprojects == 5) %>%
        distinct(taxon, bioproject_id) %>%
        mutate(prevalence_group = "specialist")
    },

    selected_generalists_specialists_table = {
      bind_rows(
        selected_generalists,
        selected_specialists
      ) %>%
        left_join(lineages) %>%
        left_join(bioproject_habitats) %>%
        arrange(prevalence_group, kingdom, phylum, order, class, family, genus)
    },

    # continious metric based on gini
    specificities = target(
      bioproject_prevalences %>%
        complete(bioproject_id, taxon, fill = list(prevalence_n = 0, prevalence_perc = 0)) %>%
        nest(-taxon) %>%
        transmute(
          taxon,
          specificity_gini = data %>% map_dbl(~ {
            .x$prevalence_perc %>% ineq::Gini(corr = TRUE)
          })
        ) %>%
        left_join(lineages),
      hpc = FALSE
    ),

    levins_sample_groupings = c("bioproject_id", "habitat", "environment_group"),

    # calculate the gold standard of generalists
    levins_indicies = target(
      {
        library(MicroNiche)

        taxa_data <-
          sub_abundances %>%
          filter(norm_method == "tss" & subset_value == "all") %>%
          unnest(data) %>%
          # pseudo counts to normalize for seqdepth
          # use a max size according to resolution of the data
          mutate(abundance = round(abundance * 10e3)) %>%
          select(sample_id, taxon, abundance) %>%
          pivot_wider(names_from = sample_id, values_from = abundance, values_fill = list(abundance = 0)) %>%
          as.data.frame()

        tibble(sample_grouping = levins_sample_groupings) %>%
          mutate(
            levins_index = sample_grouping %>% map(possibly(~ {
              sample_data <-
                taxa_data %>%
                colnames() %>%
                tibble(sample_id = .) %>%
                left_join(samples) %>%
                pluck(.x)

              n_sample_groups <- sample_data %>%
                unique() %>%
                length()

              levins.Bn(taxa_data, n_sample_groups, sample_data) %>%
                as_tibble(rownames = "taxon")
            }, NA))
          )
      },
      dynamic = map(levins_sample_groupings)
    ),

    levins_plot = target(
      {
        levins_indicies %>%
          unnest(levins_index) %>%
          left_join(selected_generalists_specialists) %>%
          filter(!is.na(prevalence_group)) %>%
          mutate(
            sample_grouping = sample_grouping %>%
              recode("bioproject_id" = "Project", "habitat" = "Habitat", "environment_group" = "Environment") %>%
              factor(levels = c("Project", "Habitat", "Environment"))
          ) %>%
          ggplot(aes(prevalence_group, Bn, color = prevalence_group)) +
          geom_boxplot() +
          scale_x_discrete(labels = str_to_sentence) +
          scale_color_prevalence_group() +
          stat_compare_means(
            method = "wilcox",
            comparisons = list(c("generalist", "specialist"))
          ) +
          facet_wrap(~sample_grouping) +
          guides(color = "none") +
          theme(panel.grid.major.x = element_blank()) +
          labs(
            x = "",
            y = TeX("Levins' niche breadth index B_n")
          )
      },
      hpc = FALSE
    ),

    # What is the most special specialist?
    selected_specialists_specialness = target(
      {
        prevalence_ranks <-
          prevalences %>%
          # specialists are defined by bioproject_id
          filter(subset_name == "bioproject_id") %>%
          select(-subset_name) %>%
          rename(bioproject_id = subset_value) %>%
          unnest(data) %>%
          # ensure ranks are comparable
          complete(bioproject_id, taxon, fill = list(prevalence_n = 0, prevalence_perc = 0)) %>%
          left_join(lineages %>% select(taxon, kingdom)) %>%
          group_by(kingdom) %>%
          mutate(prevalence_rank = rank(prevalence_perc)) %>%
          mutate(prevalence_rank = prevalence_rank / max(prevalence_rank))

        selected_specialists %>%
          mutate(
            specialness = list(bioproject_id, taxon) %>% pmap_dbl(~ {
              # special iff high prev in this project and low prev in an average other project
              this_project_score <-
                prevalence_ranks %>%
                filter(taxon == .y & bioproject_id == .x) %>%
                pull(prevalence_rank) %>%
                first()

              other_projects_score <-
                prevalence_ranks %>%
                filter(taxon == .y & bioproject_id != .x) %>%
                pull(prevalence_rank) %>%
                mean()

              this_project_score - other_projects_score
            })
          ) %>%
          left_join(lineages %>% select(taxon, kingdom)) %>%
          arrange(kingdom, -specialness)
      },
      hpc = FALSE
    ),

    generalists_min_prevalence_perc = target(
      {
        unique(common_generalists$min_prevalence_perc)
      },
      hpc = FALSE
    ),

    generalists_alphadiv_test = target(
      {
        # Is median alpha diversity in samples without n generalists lower than
        # in samples without n random taxa ?

        lacking_taxa <-
          common_generalists %>%
          filter(is_generalist) %>%
          distinct(min_prevalence_perc, taxon) %>%
          nest(taxon) %>%
          transmute(min_prevalence_perc, taxa = data %>% map(~ .x$taxon))

        random <-
          tibble(trail = seq(10e3)) %>%
          expand_grid(min_prevalence_perc = generalists_min_prevalence_perc) %>%
          left_join(lacking_taxa) %>%
          mutate(
            data = list(trail, min_prevalence_perc, taxa) %>% pmap(~ {
              get_median_alpha_div_of_samples_lacking_taxa(
                permute = TRUE,
                seed = ..1,
                taxa = ..3,
                sub_abundances = sub_abundances,
                alphadiv = alphadiv
              )
            })
          ) %>%
          unnest(data) %>%
          select(-taxa) %>%
          rename(random = median_alpha_div)

        truth <-
          random %>%
          distinct(min_prevalence_perc) %>%
          left_join(lacking_taxa) %>%
          mutate(
            data = taxa %>% map(~ {
              get_median_alpha_div_of_samples_lacking_taxa(
                permute = FALSE,
                taxa = ..1,
                sub_abundances = sub_abundances,
                alphadiv = alphadiv
              )
            })
          ) %>%
          select(-taxa) %>%
          unnest(data) %>%
          rename(truth = median_alpha_div)

        random %>%
          inner_join(truth) %>%
          group_by(min_prevalence_perc, kingdom, alphadiv_metric) %>%
          # only need to look at samples without the taxa
          # because sum of alpha diversity of both groups are const.
          filter(!has_generalist) %>%
          mutate(false_positive = random <= truth) %>%
          summarise(p.value = sum(false_positive) / n())
      },
      dynamic = map(generalists_min_prevalence_perc)
    ),

    # have generalists higher abundance than specialists using different strictness criteria?
    generalists_specialists_abundance_test = target(
      {
        full_join(
          common_generalists %>%
            filter(is_generalist) %>%
            select(min_prevalence_perc, taxon) %>%
            nest(taxon) %>%
            rename(generalists = data),
          specialists %>%
            filter(is_specialist) %>%
            select(min_prevalence_perc_this_bioproject, max_prevalence_perc_other_bioprojects, taxon, bioproject_id) %>%
            nest(taxon, bioproject_id) %>%
            rename(specialists = data),
          # cross join i. e. all combinations of generalists and specialists thresholds
          by = character()
        ) %>%
          filter(
            # get relevant search space
            min_prevalence_perc %>% between(20, 60) &
              min_prevalence_perc_this_bioproject %>% between(20, 60) &
              max_prevalence_perc_other_bioprojects %>% between(0, 20)
          ) %>%
          expand_grid(kingdom = kingdoms_groups) %>%
          transmute(
            min_prevalence_perc, min_prevalence_perc_this_bioproject, max_prevalence_perc_other_bioprojects,
            test = list(generalists, specialists, kingdom) %>% pmap(~ {
              abundances <-
                sub_abundances %>%
                filter(subset_value == "all" & norm_method == "tss") %>%
                pull(data) %>%
                first() %>%
                filter(kingdom == ..3)

              bind_rows(
                abundances %>% semi_join(.x, by = "taxon") %>% transmute(sample_id, taxon, abundance, group = "generalist"),
                abundances %>% semi_join(.y, by = c("taxon", "bioproject_id")) %>% transmute(sample_id, taxon, abundance, group = "specialist")
              ) %>%
                wilcox.test(abundance ~ group, data = .) %>%
                tidy()
            })
          ) %>%
          unnest(test) %>%
          count(is_significant = p.value < 0.05, name = "n_thresholds_sets") %>%
          complete(is_significant = c(TRUE, FALSE), fill = list(n_thresholds_sets = 0))
      },
      hpc = FALSE
    ),

    generalists_common_coabundance_graph = target(
      {
        node_majority_edge_group <- function(graph, cur_taxon) {
          graph %>%
            activate(nodes) %>%
            filter(node_is_adjacent(taxon == cur_taxon)) %>%
            activate(edges) %>%
            as_tibble() %>%
            count(environment_group) %>%
            arrange(-n) %>%
            pull(environment_group) %>%
            first()
        }

        taxa_prevalent_all_environment_groups <-
          prevalences %>%
          filter(subset_name == "environment_group") %>%
          unnest(data) %>%
          filter(prevalence_n > 0) %>%
          select(subset_value, taxon, prevalence_n) %>%
          pivot_wider(names_from = subset_value, values_from = prevalence_n) %>%
          filter(!is.na(host) & !is.na(aquatic) & !is.na(soil)) %>%
          pull(taxon)

        graph <-
          graphs %>%
          filter(subset_name == "environment_group") %>%
          rename(environment_group = subset_value) %>%
          filter(cor_method == "sparcc" & norm_method == "raw") %>%
          mutate(
            cor_graph = cor_graph %>% map(~ .x %>%
              activate(edges) %>%
              as_tibble())
          ) %>%
          unnest(cor_graph) %>%
          select(from_taxon, to_taxon, environment_group, estimate, q.value) %>%
          filter(
            q.value < 0.05 &
              abs(estimate) > generalists_min_abs_correlation &
              from_taxon %in% taxa_prevalent_all_environment_groups &
              to_taxon %in% taxa_prevalent_all_environment_groups
          ) %>%
          mutate(estimate_direction = estimate %>% map_chr(~ ifelse(sign(.x) == 1, "positive", "negative"))) %>%
          group_by(from_taxon, to_taxon, estimate_direction) %>%
          arrange(from_taxon, to_taxon) %>%
          filter(n() >= 2) %>% # present in most environments
          pivot_wider(names_from = environment_group, values_from = c(estimate, q.value)) %>%
          tbl_graph(edges = .) %>%
          activate(nodes) %>%
          mutate(taxon = name) %>%
          left_join(lineages) %>%
          left_join(selected_generalists) %>%
          annotate_node_attributes_in_edges() %>%
          activate(edges) %>%
          mutate(
            environment_group = case_when(
              !is.na(estimate_host) & !is.na(estimate_aquatic) & !is.na(estimate_soil) ~ "host-aquatic-soil",
              !is.na(estimate_host) & !is.na(estimate_aquatic) ~ "host-aquatic",
              !is.na(estimate_host) & !is.na(estimate_soil) ~ "host-soil",
              !is.na(estimate_aquatic) & !is.na(estimate_soil) ~ "aquatic-soil",
              !is.na(estimate_host) & !is.na(estimate_soil) ~ "host-soil"
            ),
            kingdoms_group = case_when(
              from_kingdom == "Bacteria" & to_kingdom == "Bacteria" ~ "intra Bacteria",
              from_kingdom == "Fungi" & to_kingdom == "Fungi" ~ "intra Fungi",
              from_kingdom != to_kingdom ~ "inter kingdom"
            ),
            involves_prevalence_group_taxon = !is.na(from_prevalence_group) | !is.na(to_prevalence_group)
          ) %>%
          activate(nodes)

        graph %>%
          mutate(
            is_generalist = (prevalence_group == "generalist") %>% replace_na(FALSE),
            majority_edge_group = taxon %>% map_chr(~ node_majority_edge_group(graph, .x)),
            taxon_group = paste0(majority_edge_group, kingdom)
          ) %>%
          topologize_graph() %>%
          annotate_node_attributes_in_edges()
      },
      hpc = FALSE
    ),

    generalists_common_coabundance_graph_plt = target(
      {
        generalists_common_coabundance_graph %>%
          activate(edges) %>%
          mutate(estimate_direction = estimate_direction %>% recode("positive" = "Positive", "negative" = "Negative")) %>%
          ggraph(layout = layout_in_circle(generalists_common_coabundance_graph, order = order(V(generalists_common_coabundance_graph)$taxon_group))) +
          geom_edge_link(aes(color = environment_group, alpha = involves_prevalence_group_taxon)) +
          geom_node_point(aes(color = kingdom, size = is_generalist)) +
          scale_edge_color_manual(values = environment_group_sets_colors, labels = str_to_sentence) +
          scale_color_kingdom(na.value = "black") +
          scale_size_manual(values = c(`TRUE` = 5, `FALSE` = 2), labels = c("No", "Yes")) +
          scale_edge_alpha_manual(values = c(`TRUE` = 0.5, `FALSE` = 0.05), guide = "none") +
          coord_fixed() +
          facet_edges(~estimate_direction) +
          guides(color = guide_legend(override.aes = list(size = 5))) +
          theme_void() +
          labs(
            color = "Kingdom",
            edge_color = "Common in environments",
            size = "Generalist"
          )
      },
      hpc = FALSE
    ),

    generalists_specialists_overall_prevalences = target(
      {
        raw_top_prevalent_generalists_specialists <-
          prevalences %>%
          filter(subset_name == "kingdom") %>%
          select(-subset_name) %>%
          separate(subset_value, into = "kingdom", sep = ",") %>%
          unnest(data) %>%
          left_join(selected_generalists_specialists) %>%
          group_by(kingdom) %>%
          # rank 1 iff most prevalent
          mutate(prevalence_rank = rank(100 - prevalence_perc)) %>%
          arrange(prevalence_rank) %>%
          filter(!is.na(prevalence_group)) %>%
          group_by(kingdom, prevalence_group) %>%
          select(kingdom, taxon, prevalence_group, raw_prevalence_rank = prevalence_rank)


        adjusted_top_prevalent_generalists_specialists <-
          prevalences %>%
          filter(subset_name == "environment_group") %>%
          select(-subset_name) %>%
          separate(subset_value, into = "environment_group", sep = ",") %>%
          unnest(data) %>%
          left_join(lineages %>% select(taxon, kingdom)) %>%

          # control for environment_group sample sizes
          group_by(kingdom, taxon) %>%
          summarise(prevalence_perc = mean(prevalence_perc)) %>%

          # rank 1 iff most prevalent
          mutate(prevalence_rank = rank(100 - prevalence_perc)) %>%
          arrange(prevalence_rank) %>%
          left_join(selected_generalists_specialists) %>%
          filter(!is.na(prevalence_group)) %>%
          group_by(kingdom, prevalence_group) %>%
          select(kingdom, taxon, prevalence_group, adjusted_prevalence_rank = prevalence_rank)


        list(
          raw_top_prevalent_generalists_specialists,
          adjusted_top_prevalent_generalists_specialists
        ) %>%
          reduce(full_join) %>%
          arrange(raw_prevalence_rank) %>%
          slice(1:5)
      },
      hpc = FALSE
    ),

    generalists_alphadiv_test_plt = {
      generalists_alphadiv_test %>%
        filter(min_prevalence_perc < 60 & min_prevalence_perc > 20) %>%
        ggplot(aes(min_prevalence_perc, p.value)) +
        geom_point(size = 3) +
        geom_line() +
        geom_hline(aes(color = "significance", yintercept = 0.05)) +
        geom_vline(aes(color = "prevalence", xintercept = 40)) +
        scale_color_manual(values = c(prevalence = "blue", significance = "red")) +
        facet_grid(alphadiv_metric ~ kingdom) +
        labs(
          x = "Generalists prevalence threshold (%)",
          y = "Permutation significance (FDR)",
          color = "Threshold"
        )
    }
  )
}
