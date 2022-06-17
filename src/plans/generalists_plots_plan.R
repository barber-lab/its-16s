get_generalists_plots_plan <- function() {
  drake_plan(
    fig1_overview = target(
      {
        prevalences_venn_plts <-
          prevalence_summary %>%
          filter(
            subset_name == "kingdom"
          ) %>%
          mutate(venn_plt = taxa %>% map2(subset_value, ~ {
            .x %>%
              enframe() %>%
              unnest(value) %>%
              mutate(is_present = TRUE) %>%
              pivot_wider(names_from = name, values_from = is_present, values_fill = list(is_present = FALSE)) %>%
              select_at(environment_groups) %>%
              ggeulerr(show_labels = FALSE) +
              scale_fill_environment_group()
          })) %>%
          ungroup() %>%
          select(subset_value, venn_plt) %>%
          deframe()

        collection_times_plt <-
          samples %>%
          semi_join(sub_abundances$data[[1]]) %>%
          # keep same plot order as habitat pictograms
          mutate(environment_group = environment_group %>% factor(levels = c("host", "soil", "aquatic"))) %>%
          ggplot(aes(collection_datetime, fill = environment_group)) +
          geom_histogram(position = "identity", binwidth = 60 * 60 * 24 * 60) +
          facet_wrap(~environment_group, ncol = 1, strip.position = "right", labeller = labeller(environment_group = str_to_sentence)) +
          scale_y_log10() +
          scale_x_datetime(
            labels = scales::date_format("%Y"),
            date_breaks = "3 years",
            minor_breaks = scales::date_breaks("1 year")
          ) +
          scale_fill_environment_group() +
          annotation_logticks(sides = "l") +
          theme_pub() +
          theme(
            panel.grid.minor.x = element_line(),
            panel.grid.major.x = element_line(),
            panel.grid.major.y = element_line()
          ) +
          theme(panel.margin.y = unit(0, "lines")) +
          labs(x = "Collection time", y = "Samples", fill = "Environment")

        geocoded_bioprojects <- read_csv("raw/meta/geocoded_bioprojects.csv")

        worldmap_tbl <-
          samples %>%
          filter(!is.na(lat)) %>%
          bind_rows(
            samples %>% filter(is.na(lat)) %>% select(-lat, -lon) %>% left_join(geocoded_bioprojects)
          ) %>%
          filter(
            !is.na(environment_group) &
              sample_id %in% sub_abundances$data[[1]]$sample_id
          ) %>%
          mutate(
            lon = lon %>% round(-0.5),
            lat = lat %>% round(-0.5)
          ) %>%
          group_by(environment_group, lon, lat) %>%
          count()

        worldmap_plt <-
          rnaturalearth::ne_countries(scale = "small", returnclass = "sf") %>%
          ggplot() +
          geom_sf(color = "lightgrey") +
          geom_point(
            data = worldmap_tbl,
            mapping = aes(lon, lat, color = environment_group, size = n),
          ) +
          scale_x_continuous(
            expand = c(0.03, 0.03),
            limits = c(worldmap_tbl$lon %>% min(na.rm = TRUE), worldmap_tbl$lon %>% max(na.rm = TRUE))
          ) +
          scale_y_continuous(
            expand = c(0.03, 0.03),
            limits = c(worldmap_tbl$lat %>% min(na.rm = TRUE), worldmap_tbl$lat %>% max(na.rm = TRUE))
          ) +
          scale_color_environment_group() +
          guides(color = guide_legend(override.aes = list(size = 3))) +
          labs(
            x = "",
            y = "",
            color = "Environment",
            size = "Samples"
          ) +
          theme_pub()

        prevalent_taxa <-
          sub_abundances$data[[1]] %>%
          left_join(lineages) %>%
          left_join(samples) %>%
          filter(abundance > 0 & kingdom %in% kingdoms_groups) %>%
          distinct(environment_group, kingdom, phylum, taxon) %>%
          mutate(
            is_prevalent = TRUE,
            phylum = case_when(
              phylum %in% names(phyla_colors) ~ phylum,
              kingdom == "Bacteria" ~ "other Bacteria",
              kingdom == "Fungi" ~ "other Fungi"
            ) %>%
              # Keep all levels in ggplot allowing guide collect
              factor(levels = names(phyla_colors))
          ) %>%
          pivot_wider(
            names_from = environment_group,
            values_from = is_prevalent,
            values_fill = list(is_prevalent = FALSE)
          ) %>%
          distinct(taxon, .keep_all = TRUE)

        common_unique_taxa_counts_per_sample_subset_plts <-
          kingdoms_groups %>%
          map(~ {
            data <-
              prevalent_taxa %>% filter(kingdom == .x)

            tibble() %>%
              ggplot(aes(y = ..count.. / nrow(data) * 100, fill = phylum)) +
              geom_bar(
                data = data %>% filter(aquatic & !soil & !host),
                mapping = aes(x = "Unique to\naquatic")
              ) +
              geom_bar(
                data = data %>% filter(!aquatic & !soil & host),
                mapping = aes(x = "Unique to\nhost")
              ) +
              geom_bar(
                data = data %>% filter(!aquatic & soil & !host),
                mapping = aes(x = "Unique to\nsoil")
              ) +
              geom_bar(
                data = data %>% filter(aquatic & host & soil),
                mapping = aes(x = "Common\nin all")
              ) +
              labs(
                x = "",
                y = str_glue("Genera (% of all {str_to_lower(.x)})")
              ) +
              theme(
                panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank()
              ) +
              scale_fill_phyla() +
              scale_y_continuous(expand = c(0, 0), limits = c(0, 45)) +
              facet_wrap(~kingdom, ncol = 1, scales = "free") +
              coord_flip() +
              labs(fill = "Phylum")
          }) %>%
          set_names(kingdoms_groups)

        common_unique_abundances <-
          prevalent_taxa %>%
          nest(-kingdom) %>%
          mutate(
            data = data %>% map(~ {
              .x %>%
                mutate(
                  group = case_when(
                    soil & host & aquatic ~ "Common\nin all",
                    soil & !host & !aquatic ~ "Unique to\nsoil",
                    !soil & !host & aquatic ~ "Unique to\naquatic",
                    !soil & host & !aquatic ~ "Unique to\nhost"
                  ),
                  environment_group = group %>% recode(
                    "Common\nin all" = "all", "Unique to\nsoil" = "soil", "Unique to\naquatic" = "aquatic",
                    "Unique to\nhost" = "host"
                  )
                ) %>%
                left_join(
                  sub_abundances$data[[2]] %>% select(sample_id, taxon, abundance)
                ) %>%
                filter(abundance > 0.01e-2) %>%
                filter(!is.na(group)) %>%
                left_join(sub_abundances$data[[2]]) %>%
                ggplot(aes(group, abundance, color = environment_group)) +
                geom_boxplot(width = 0.5) +
                stat_compare_means(
                  method = "wilcox",
                  comparisons = list(
                    c("Common\nin all", "Unique to\naquatic"),
                    c("Common\nin all", "Unique to\nhost"),
                    c("Common\nin all", "Unique to\nsoil")
                  )
                ) +
                scale_y_log10(expand = c(0, 0.5)) +
                scale_color_environment_group() +
                annotation_logticks(sides = "l") +
                guides(color = FALSE) +
                theme(
                  panel.grid.minor.x = element_blank(),
                  panel.grid.major.x = element_blank()
                ) +
                labs(x = "", y = "Abundance (TSS)")
            })
          ) %>%
          deframe()

        wrap_plots(
          collection_times_plt + guides(fill = FALSE) + labs(tag = "A"),
          rarefactions_bioproject_id_plt + labs(tag = "B"),
          worldmap_plt + labs(tag = "C"),
          wrap_plots(
            prevalences_venn_plts$Bacteria + facet_wrap(~"Bacteria") + guides(fill = FALSE) + labs(tag = "D"),
            prevalences_venn_plts$Fungi + facet_wrap(~"Fungi") + guides(fill = FALSE)
          ),
          common_unique_taxa_counts_per_sample_subset_plts$Bacteria + labs(tag = "E"),
          common_unique_taxa_counts_per_sample_subset_plts$Fungi + labs(tag = "F"),
          common_unique_abundances$Bacteria + facet_wrap(~"Bacteria") + labs(tag = "G"),
          common_unique_abundances$Fungi + facet_wrap(~"Fungi") + labs(tag = "H"),
          ncol = 2,
          widths = c(1, 1),
          heights = c(2.5, 2, 2, 2.5),
          guides = "collect"
        )
      },
      hpc = FALSE
    ),
    generalists_and_specialists_abundances_plt = target(
      {
        kingdoms_groups %>%
          set_names(kingdoms_groups) %>%
          purrr::map(function(current_kingdom) {
            sub_abundances %>%
              filter(subset_value == "all" & norm_method == "tss") %>%
              pull(data) %>%
              first() %>%
              inner_join(selected_generalists_specialists) %>%
              filter(kingdom == current_kingdom & !is.na(prevalence_group)) %>%
              filter(abundance > 0.01e-2) %>%
              ggplot(aes(x = prevalence_group, y = abundance, color = prevalence_group)) +
              geom_boxplot(width = 0.5) +
              facet_grid(~kingdom) +
              scale_x_discrete(labels = str_to_sentence) +
              scale_y_log10(expand = c(0, 1)) +
              stat_compare_means(
                method = "wilcox",
                comparisons = list(c("generalist", "specialist"))
              ) +
              annotation_logticks(sides = "l") +
              scale_color_prevalence_group() +
              theme(
                legend.position = "right",
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                panel.grid.minor.y = element_blank()
              ) +
              guides(color = FALSE) +
              labs(
                x = "",
                y = "Abundance (TSS)"
              )
          })
      },
      hpc = FALSE
    ),
    fig3_alphadiv = target(
      {
        # taxa detected in a sample given the generalists abundance threshold
        identifications <-
          sub_abundances$data[[2]] %>%
          transmute(
            bioproject_id,
            sample_id,
            taxon,
            is_identified =  abundance > generalists_min_prevalence_abundance
          )

        samples_having_any_generalists <-
          sub_abundances %>%
          filter(subset_name == "kingdom" & norm_method == "tss") %>%
          transmute(kingdom = subset_value, data) %>%
          select(-kingdom) %>%
          unnest(data) %>%
          left_join(
            selected_generalists_specialists %>% filter(prevalence_group == "generalist")
          ) %>%
          filter(abundance > 0) %>%
          group_by(sample_id, kingdom) %>%
          summarise(has_generalist = any(prevalence_group == "generalist")) %>%
          mutate(has_generalist = has_generalist %>% replace_na(FALSE)) %>%
          transmute(sample_id, kingdom, group = ifelse(has_generalist, "with generalists", "without generalists"))

        samples_having_anything <-
          sub_abundances$data[[2]] %>%
          distinct(sample_id, kingdom) %>%
          transmute(sample_id, kingdom, group = "all")

        # see file src/lab/77-perm-test-aphadiv-generalists.R
        alphadiv_generalists_tests <- tribble(
          ~kingdom, ~alphadiv_metric, ~group1, ~group2, ~label, ~y.position,
          "Bacteria", "Chao1", "with generalists", "without generalists", "\u2605\u2605\u2605", 2000,
          "Bacteria", "Shannon", "with generalists", "without generalists", "\u2605\u2605\u2605", 5.5,
          "Fungi", "Chao1", "with generalists", "without generalists", "\u2605", 1000,
          "Fungi", "Shannon", "with generalists", "without generalists", "\u2605", 4.5,
        )

        generalists_alphadiv_plts <-
          samples_having_any_generalists %>%
          inner_join(
            alphadiv %>%
              filter(subset_name == "kingdom") %>%
              select(kingdom = subset_value, data) %>%
              unnest(data)
          ) %>%
          pivot_longer(cols = c(Chao1, Shannon), names_to = "alphadiv_metric") %>%
          mutate(group = group %>% factor(
            levels = c("without generalists", "with generalists")
          )) %>%
          nest(-c(kingdom, alphadiv_metric)) %>%
          mutate(
            plt = list(kingdom, alphadiv_metric, data) %>% pmap(~ {
              plt <-
                ..3 %>%
                mutate(
                  facet_x = .x,
                  facet_y = .y,
                  kingdom = .x,
                  has_any_generalist = group %>% str_detect("without") %>% ifelse(FALSE, TRUE)
                ) %>%
                ggplot(aes(group, value)) +
                geom_boxplot(aes(alpha = group), width = 0.3, fill = "#2d51a5", color = "#2d51a5") +
                scale_alpha_manual(values = c("with generalists" = 0.5, "without generalists" = 0)) +
                scale_x_discrete(label = function(x) x %>% str_to_sentence() %>% str_replace(" ", "\n")) +
                stat_pvalue_manual(
                  data = {
                    data <-
                      alphadiv_generalists_tests %>%
                      filter(kingdom == .x & alphadiv_metric == .y)

                    if (.y == "Chao1") {
                      data %>% mutate(y.position = log10(y.position))
                    } else {
                      data
                    }
                  }
                ) +
                facet_wrap(~facet_x, scales = "free", nrow = 1) +
                scale_y_continuous(expand = expand_scale(mult = c(0, 0.15))) +
                guides(color = FALSE, alpha = FALSE) +
                theme(
                  panel.grid.minor.x = element_blank(),
                  panel.grid.major.x = element_blank()
                ) +
                labs(x = "", y = .y, alpha = "sample group")

              if (.y == "Chao1") {
                plt +
                  scale_y_log10(expand = expansion(mult = c(0, 0.15))) +
                  expand_limits(x = 1) +
                  annotation_logticks(sides = "l")
              } else {
                plt
              }
            })
          ) %>%
          arrange(alphadiv_metric, kingdom)

        selected_specialists <-
          specialists %>%
          filter(is_specialist & min_prevalence_perc_this_bioproject == 40 & max_prevalence_perc_other_bioprojects == 5) %>%
          select(bioproject_id, taxon) %>%
          mutate(is_specialist = TRUE)

        samples_detected_taxa <-
          sub_abundances %>%
          filter(subset_name == "all" & norm_method == "tss") %>%
          unnest(data) %>%
          transmute(
            sample_id, bioproject_id, taxon,
            detected = abundance > generalists_min_prevalence_abundance
          )

        samples_having_any_specialists <-
          samples_detected_taxa %>%
          left_join(selected_specialists) %>%
          mutate(is_specialist = is_specialist %>% replace_na(FALSE)) %>%
          left_join(lineages) %>%
          group_by(sample_id, kingdom) %>%
          summarise(
            has_any_specialist = any(detected & is_specialist) %>% replace_na(FALSE)
          )

        # need samples with and without any specialist for any kingdom
        # discard other habitats
        selected_habitats <-
          samples_having_any_specialists %>%
          filter(has_any_specialist) %>%
          left_join(samples %>% select(sample_id, bioproject_id)) %>%
          left_join(bioproject_habitats) %>%
          ungroup() %>%
          distinct(habitat) %>%
          pull(habitat)

        sample_groups <-
          habitats %>%
          filter(habitat %in% selected_habitats) %>%
          expand_grid(has_any_specialist = c(TRUE, FALSE)) %>%
          arrange(environment_group, habitat, has_any_specialist) %>%
          mutate(
            sample_group = has_any_specialist %>% ifelse("with specialists", "without specialists") %>% paste0(habitat, " ", .) %>% fct_inorder()
          )

        alphadiv_specialists_habitats <-
          alphadiv %>%
          filter(subset_name == "kingdom") %>%
          rename(kingdom = subset_value) %>%
          unnest(data) %>%
          pivot_longer(c(Chao1, Shannon), names_to = "alphadiv_metric", values_to = "alphadiv_value") %>%
          left_join(samples %>% select(sample_id, bioproject_id)) %>%
          inner_join(samples_having_any_specialists) %>%
          # alternative definition: left_join(samples %>% select(sample_id, habitat))
          left_join(bioproject_habitats) %>%
          full_join(sample_groups) %>%
          filter(!is.na(sample_group) & !is.na(kingdom))

        alphadiv_specialists_habitats_tests <-
          alphadiv_specialists_habitats %>%
          nest(-c(kingdom, habitat, alphadiv_metric)) %>%
          mutate(
            data = data %>% map(possibly(~ {
              .x %>%
                wilcox.test(alphadiv_value ~ sample_group, data = .) %>%
                tidy()
            }, NA))
          ) %>%
          unnest(data) %>%
          mutate(q.value = p.value %>% p.adjust(method = "fdr")) %>%
          filter(q.value < 0.05) %>%
          transmute(
            kingdom,
            habitat,
            group1 = paste0(habitat, " with specialists"),
            group2 = paste0(habitat, " without specialists"),
            p.value,
            q.value,
            label = q.value %>% significance_label()
          ) %>%
          left_join(
            alphadiv_specialists_habitats %>%
              group_by(habitat, alphadiv_metric, kingdom) %>%
              summarise(y.position = 1.1 * max(alphadiv_value))
          )

        specialists_alphadiv_plts <-
          alphadiv_specialists_habitats %>%
          # Remove sample group all from legend
          mutate(environment_group = environment_group %>% factor(levels = c("host", "aquatic", "soil"))) %>%
          nest(-c(kingdom, alphadiv_metric)) %>%
          transmute(
            kingdom, alphadiv_metric,
            plt = list(data, kingdom, alphadiv_metric) %>% pmap(~ {
              plt <-
                .x %>%
                mutate(kingdom = .y) %>%
                ggplot(aes(sample_group, alphadiv_value)) +
                geom_boxplot(aes(color = environment_group, fill = environment_group, alpha = has_any_specialist)) +
                scale_x_discrete(
                  drop = FALSE, expand = expansion(0.1, 0),
                  labels = function(x) x %>% str_remove("with(out)? specialists") %>% str_to_sentence()
                ) +
                scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
                scale_color_environment_group(drop = TRUE) +
                scale_fill_environment_group(drop = TRUE) +
                scale_alpha_manual(values = c("TRUE" = 0.5, "FALSE" = 0)) +
                stat_pvalue_manual(
                  data = {
                    data <-
                      alphadiv_specialists_habitats_tests %>%
                      filter(kingdom == .y & alphadiv_metric == ..3)

                    if (..3 == "Chao1") {
                      data %>% mutate(y.position = log10(y.position))
                    } else {
                      data
                    }
                  },
                  label = "label",
                  label.size = 2
                ) +
                facet_wrap(~kingdom) +
                theme(
                  axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                  panel.grid.major.x = element_blank()
                ) +
                guides(alpha = FALSE, color = FALSE) +
                labs(x = "", y = ..3, fill = "Environment")

              if (..3 == "Chao1") {
                plt +
                  scale_y_log10(expand = expansion(mult = c(0, 0.1))) +
                  expand_limits(x = 1) +
                  annotation_logticks(sides = "l")
              } else {
                plt
              }
            })
          ) %>%
          arrange(alphadiv_metric, kingdom)


        layout <- "
        ABCD
        EEFF
        GGHH
        "
        generalists_alphadiv_plts %>%
          arrange(kingdom, alphadiv_metric) %>%
          pull(plt) %>%
          c(specialists_alphadiv_plts$plt) %>%
          wrap_plots(design = layout, guides = "collect", heights = c(1, 1.5, 1.5)) +
          plot_annotation(tag_levels = "A")
      },
      hpc = FALSE
    ),

    levins_plt = target(
      {
        levins_indicies %>%
          unnest(levins_index) %>%
          left_join(selected_generalists_specialists) %>%
          filter(!is.na(prevalence_group)) %>%
          mutate(sample_grouping = sample_grouping %>% factor(levels = c("bioproject_id", "habitat", "environment_group"))) %>%
          ggplot(aes(prevalence_group, Bn, color = prevalence_group)) +
          geom_boxplot() +
          scale_color_prevalence_group() +
          stat_compare_means(
            method = "wilcox",
            comparisons = list(c("generalist", "specialist"))
          ) +
          facet_wrap(~sample_grouping) +
          guides(color = "none") +
          theme(panel.grid.major.x = element_blank()) +
          labs(y = TeX("Levins' niche breadth index B_n"))
      },
      hpc = FALSE
    ),

    generalists_min_abs_correlation = 0.2,

    correlation_pos_neg_plts = target(
      {
        edges <-
          graphs %>%
          filter(subset_name %in% c("all", "environment_group") & cor_method == "sparcc") %>%
          transmute(
            subset_value,
            data = cor_graph %>% map(~ {
              .x %>%
                activate(edges) %>%
                # treat small correlation e.g. +0.001 and -0.001 the same
                # to reduce noise
                filter(q.value < 0.05 & abs(estimate) > generalists_min_abs_correlation) %>%
                as_tibble() %>%
                left_join(selected_generalists_specialists %>% rename_all(~ paste0("from_", .))) %>%
                left_join(selected_generalists_specialists %>% rename_all(~ paste0("to_", .))) %>%
                select(contains("kingdom"), contains("taxon"), estimate, contains("prevalence_group"))
            })
          ) %>%
          unnest() %>%
          mutate(
            kingdom = case_when(
              from_kingdom == "Bacteria" & to_kingdom == "Bacteria" ~ "intra bacteria",
              from_kingdom == "Fungi" & to_kingdom == "Fungi" ~ "intra fungi",
              from_kingdom != to_kingdom ~ "inter kingdom"
            ),
            prevalence_group = case_when(
              !is.na(from_prevalence_group) ~ from_prevalence_group,
              !is.na(to_prevalence_group) ~ to_prevalence_group,
              TRUE ~ "other",
            ) %>%
              recode(
                "generalist" = "with generalists",
                "specialist" = "with specialists"
              ),
            edge_group = paste(kingdom, prevalence_group, sep = "\n"),
            edge_sub_group = paste0(sign(estimate), prevalence_group, kingdom),
            special_taxon = case_when(
              !is.na(from_prevalence_group) ~ from_taxon,
              !is.na(to_prevalence_group) ~ to_taxon,
              TRUE ~ from_taxon
            ),
            other_taxon = case_when(
              to_taxon == special_taxon ~ from_taxon,
              TRUE ~ to_taxon
            )
          )

        box_plt <-
          edges %>%
          ggplot(aes(edge_group, estimate)) +
          geom_rect(xmin = -Inf, xmax = Inf, ymin = -min_abs_correlation, ymax = min_abs_correlation, fill = "lightgrey") +
          geom_hline(yintercept = 0) +
          geom_tile(
            # show all combinations
            data = edges %>% distinct(edge_group, edge_sub_group, kingdom),
            mapping = aes(fill = kingdom, y = -1),
            height = 0.2
          ) +
          scale_fill_aaas() +
          scale_fill_manual(values = c(
            "intra bacteria" =  as.character(kingdoms_colors["Bacteria"]),
            "inter kingdom" = "black",
            "intra fungi" =   as.character(kingdoms_colors["Fungi"])
          )) +
          geom_boxplot(
            mapping = aes(color = prevalence_group, group = edge_sub_group),
            position = position_dodge(0),
            varwidth = TRUE
          ) +
          scale_color_manual(values = c(
            "with generalists" = as.character(prevalence_group_colors["generalist"]),
            "with specialists" = as.character(prevalence_group_colors["specialist"]),
            "other" = "black"
          )) +
          stat_compare_means(
            method = "wilcox",
            comparisons = list(
              c("intra bacteria\nwith generalists", "intra bacteria\nother"),
              c("intra fungi\nwith generalists", "intra fungi\nother"),
              c("inter kingdom\nwith generalists", "inter kingdom\nother"),

              c("intra bacteria\nwith generalists", "intra fungi\nwith generalists"),
              c("intra bacteria\nwith generalists", "inter kingdom\nwith generalists")
            )
          ) +
          scale_y_continuous(expand = expansion(c(0, 0.3))) +
          scale_x_discrete(expand = c(0, 0)) +
          facet_wrap(~subset_value, ncol = 1) +
          theme(
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.text.x = element_blank()
          ) +
          labs(
            x = "Edge group",
            y = str_glue("SparCC correlation ($|r| > {min_abs_correlation}$)") %>% latex2exp::TeX(),
            color = "",
            fill = ""
          )

        bar_plt <-
          edges %>%
          mutate(
            edge_group = edge_group %>% str_replace_all("[\n]", " "),
            edge_direction = ifelse(estimate > 0, "positive", "negative")
          ) %>%
          ggplot(aes(edge_group, fill = edge_direction)) +
          geom_bar(position = position_dodge()) +
          geom_text(stat = "count", aes(label = ..count.., color = edge_direction), position = position_dodge(1), hjust = -0.1) +
          scale_fill_aaas() +
          scale_color_aaas() +
          coord_flip() +
          facet_wrap(~subset_value, ncol = 1, scales = "free_y") +
          scale_y_continuous(expand = expansion(c(0, 0.3))) +
          labs(
            x = "Edge group",
            y = "Edges"
          )

        list(
          box_plt = box_plt,
          bar_plt = bar_plt
        )
      },
      hpc = FALSE
    ),

    correlation_common_taxa_across_environment_group_topology_plt = target(
      {
        # Does coabundance differ across environments?
        # Do node topology instead of drawing networks due to > 30k edges

        taxa_prevalent_all_environment_groups <-
          prevalences %>%
          filter(subset_name == "environment_group") %>%
          unnest(data) %>%
          filter(prevalence_n > 0) %>%
          select(subset_value, taxon, prevalence_n) %>%
          pivot_wider(names_from = subset_value, values_from = prevalence_n) %>%
          filter(!is.na(host) & !is.na(aquatic) & !is.na(soil)) %>%
          pull(taxon)

        nodes <-
          graphs %>%
          filter(subset_name %in% c("environment_group") & cor_method == "sparcc") %>%
          transmute(
            subset_value,
            data = cor_graph %>% map(~ {
              .x %>%
                activate(edges) %>%
                filter(
                  q.value < 0.05 &
                    # treat small correlation e.g. +0.001 and -0.001 the same to reduce noise
                    abs(estimate) > generalists_min_abs_correlation &
                    # only look at taxa prevalent in all environmentrs to make them comparable
                    # i.e. remove trivial case of diff coabundance just because of diff abundance
                    from_taxon %in% taxa_prevalent_all_environment_groups &
                    to_taxon %in% taxa_prevalent_all_environment_groups
                ) %>%
                activate(nodes) %>%
                # filtering edges requires to re-calculate topology
                topologize_graph() %>%
                left_join(selected_generalists_specialists) %>%
                as_tibble()
            })
          ) %>%
          unnest() %>%
          mutate(
            prevalence_group = prevalence_group %>% replace_na("other"),
            taxon_group = paste(prevalence_group, kingdom, sep = "\n")
          ) %>%
          pivot_longer(c(degree, closeness, betweeness, hub), names_to = "node_topology_name", values_to = "node_topology_value")

        nodes %>%
          filter(!node_topology_name %in% c("closeness", "degree")) %>%
          nest(-node_topology_name) %>%
          mutate(
            y_log = node_topology_name == "betweeness",
            node_topology_name = node_topology_name %>% recode("hub" = "Kleinberg's hub centrality")
          ) %>%
          mutate(
            plt = list(data, y_log, node_topology_name) %>% pmap(~ {
              plt <-
                .x %>%
                ggplot(aes(subset_value, node_topology_value, color = subset_value)) +
                geom_boxplot() +
                scale_x_discrete(labels = str_to_sentence) +
                scale_color_environment_group() +
                stat_compare_means_environment_group(method = "wilcox") +
                theme(panel.grid.major.x = element_blank()) +
                labs(
                  x = "",
                  y = ..3 %>% str_to_sentence(),
                  color = "Environment"
                ) +
                guides(color = FALSE)

              if (.y) {
                plt + scale_y_log10() + annotation_logticks(sides = "l")
              } else {
                plt
              }
            })
          ) %>%
          pull(plt) %>%
          wrap_plots(nrow = 1, guides = "collect")
      },
      hpc = FALSE
    ),

    correlation_nodes_box_plt = target(
      {
        nodes <-
          graphs %>%
          filter(subset_name %in% c("all", "environment_group") & cor_method == "sparcc") %>%
          transmute(
            subset_value,
            data = cor_graph %>% map(~ {
              .x %>%
                activate(edges) %>%
                # treat small correlation e.g. +0.001 and -0.001 the same
                # to reduce noise
                filter(q.value < 0.05 & abs(estimate) > generalists_min_abs_correlation) %>%
                activate(nodes) %>%
                # filtering edges requires to re-calculate topology
                topologize_graph() %>%
                left_join(selected_generalists_specialists) %>%
                as_tibble()
            })
          ) %>%
          unnest() %>%
          mutate(
            prevalence_group = prevalence_group %>% replace_na("other"),
            taxon_group = paste(prevalence_group, kingdom, sep = "\n")
          ) %>%
          pivot_longer(c(degree, closeness, betweeness, hub), names_to = "node_topology_name", values_to = "node_topology_value")

        nodes %>%
          nest(-c(subset_value, node_topology_name)) %>%
          mutate(
            log_transform = node_topology_name %in% c("betweeness", "degree"),
            plt = list(subset_value, node_topology_name, log_transform, data) %>% pmap(~ {
              plt <-
                ..4 %>%
                mutate(
                  facet_x = ..1,
                  facet_y = ..2,
                  # keep all levels to combine guide
                  prevalence_group = prevalence_group %>% factor(levels = prevalence_group_colors %>% names() %>% c("other"))
                ) %>%
                ggplot(aes(taxon_group, node_topology_value)) +
                geom_boxplot(aes(color = prevalence_group)) +
                scale_color_manual(values = prevalence_group_colors %>% append(c("other" = "black")), drop = FALSE) +
                stat_compare_means(
                  method = "wilcox",
                  comparisons = {
                    # compare only those who are available
                    # e.g. specialists might fail in some plots
                    ..4$taxon_group %>%
                      unique() %>%
                      combn(2) %>%
                      t() %>%
                      as_tibble() %>%
                      filter(V1 %>% str_detect("generalist") | V2 %>% str_detect("generalist")) %>%
                      transmute(comp = list(V1, V2) %>% pmap(~ c(.x, .y))) %>%
                      pull(comp)
                  }
                ) +
                theme(
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank()
                ) +
                facet_wrap(~facet_x, scale = "free") +
                labs(x = "", color = "Taxon group", y = ..2)

              if (..3) {
                plt <-
                  plt +
                  scale_y_log10(expand = expansion(c(0, 0.1))) +
                  annotation_logticks(sides = "l")
              }

              plt
            })
          ) %>%
          pull(plt) %>%
          wrap_plots(guides = "collect", ncol = 4)
      },
      hpc = FALSE
    ),

    fig6_heatmap = target({
      mat <-
        sub_abundances %>%
        filter(subset_name == "kingdom" & norm_method == "tss") %>%
        unnest(data) %>%
        filter(abundance > generalists_min_prevalence_abundance) %>%
        filter(!is.na(environment_group)) %>%
        left_join(selected_generalists_specialists) %>%
        semi_join(
          prevalences %>%
            filter(subset_value == "all") %>%
            unnest(data) %>%
            full_join(selected_generalists_specialists) %>%
            filter(prevalence_perc > 5) %>%
            select(taxon)
        ) %>%
        filter(taxon != "Blattella germanica (German cockroach)") %>%
        arrange(environment_group, habitat, bioproject_id, kingdom, prevalence_group) %>%
        select(sample_id, taxon, abundance) %>%
        # normalization
        group_by(taxon) %>%
        mutate(abundance = scale(log10(abundance))) %>%
        pivot_wider(names_from = sample_id, values_from = abundance, values_fill = list(abundance = -10)) %>%
        column_to_rownames("taxon") %>%
        as.matrix()

      samples_annot_tbl <-
        mat %>%
        colnames() %>%
        tibble(sample_id = .) %>%
        left_join(samples)

      taxa_annot_tbl <-
        mat %>%
        rownames() %>%
        tibble(taxon = .) %>%
        left_join(lineages) %>%
        left_join(
          prevalences %>%
            filter(subset_value == "all") %>%
            unnest(data)
        ) %>%
        left_join(selected_generalists_specialists) %>%
        mutate(
          position = row_number(),
          fontcolor = prevalence_group %>% map_chr(~ prevalence_group_colors[.x]),
          colname = ifelse(is.na(prevalence_group), "", taxon),
          short_taxon = taxon %>% shorten_taxon(),
          phylum = case_when(
            phylum %in% names(phyla_colors) ~ phylum,
            kingdom == "Bacteria" ~ "other Bacteria",
            kingdom == "Fungi" ~ "other Fungi"
          ) %>% factor(levels = names(phyla_colors))
        )

      label_font_size <- 10

      Heatmap(
        matrix = mat,
        name = "Abundance [z(log10(TSS))]",
        col = circlize::colorRamp2(c(-10, 0, 10), viridis::viridis(3)),
        na_col = "white",
        use_raster = TRUE,

        # UPGMA clustering
        clustering_method_rows = "average",
        clustering_method_columns = "average",
        clustering_distance_rows = "manhattan",
        clustering_distance_columns = "manhattan",

        row_split = taxa_annot_tbl$kingdom,
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        top_annotation = HeatmapAnnotation(
          which = "column",
          annotation_name_gp = gpar(fontsize = label_font_size),
          Environment = samples_annot_tbl$environment_group,
          col = list(
            Environment = environment_group_colors %>% simplify()
          ),
          annotation_legend_param = list(
            Environment = list(labels_gp = gpar(fontsize = label_font_size))
          )
        ),
        right_annotation = HeatmapAnnotation(
          which = "row",
          annotation_name_gp = gpar(fontsize = label_font_size),
          Phylum = taxa_annot_tbl$phylum,
          `Prevalence (%)` = taxa_annot_tbl$prevalence_perc,
          Group = taxa_annot_tbl$prevalence_group,
          `Prevalence taxon` = anno_mark(
            at = taxa_annot_tbl %>% filter(!is.na(prevalence_group)) %>% pull(position),
            labels = taxa_annot_tbl %>% filter(!is.na(prevalence_group)) %>% pull(short_taxon),
            lines_gp = gpar(lwd = 0.2),
            labels_gp = gpar(
              fontsize = 6,
              fontface = 3,
              col = {
                taxa_annot_tbl %>%
                  filter(!is.na(prevalence_group)) %>%
                  transmute(taxon, prevalence_group_colors[prevalence_group]) %>%
                  deframe()
              }
            )
          ),
          col = list(
            Phylum = phyla_colors %>% simplify(),
            Kingdom = kingdoms_colors %>% simplify(),
            Group = prevalence_group_colors %>% simplify(),
            `Prevalence (%)` = circlize::colorRamp2(c(0, 100), c("white", "black"))
          ),
          annotation_legend_param = list(
            Phylum = list(nrow = 2, labels_gp = gpar(fontsize = label_font_size)),
            `Prevalence (%)` = list(labels_gp = gpar(fontsize = label_font_size)),
            Group = list(labels_gp = gpar(fontsize = label_font_size))
          )
        )
      )
    })
  )
}
