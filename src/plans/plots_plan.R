get_plots_plan <- function() {
  drake_plan(
    #
    # Parameters ----
    #
    min_heatmap_prevalence_perc = 10,

    fig_prevalance_abundance = target(
      {
        prevalence_abundance_plt <- plot_prevalences_abundances(prevalences, abundances, samples, lineages)
        prevalence_abundance_plt$patches[[1]][[1]] <- prevalence_abundance_plt$patches[[1]][[1]] + labs(tag = "C")

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

        common_unique_taxa_counts_plt <- plot_common_unique_taxa_counts_per_sample_subset(
          abundances = abundances, lineages = lineages, samples = samples, kingdoms_groups = kingdoms_groups
        )

        design <- "ABD
                   CCE"

        venn_bar_plts <-
          wrap_plots(
            prevalences_venn_plts$Bacteria + facet_wrap(~"Bacteria") + theme(legend.position = "none") + labs(tag = "A"),
            prevalences_venn_plts$Fungi + facet_wrap(~"Fungi") + theme(legend.position = "none"),
            guide_area(),
            common_unique_taxa_counts_plt$Bacteria +
              guides(fill = guide_legend(
                ncol = 2, title.position = "top",
                label.theme = element_text(size = 8),
                title.theme = element_text(size = 8),
                keywidth = 0.5,
                keyheight = 0.5
              )) +
              labs(tag = "B"),
            common_unique_taxa_counts_plt$Fungi + guides(fill = guide_none()),
            widths = c(1, 1, 1)
          ) +
          plot_layout(design = design, guides = "collect")

        cowplot::plot_grid(venn_bar_plts, prevalence_abundance_plt,
          ncol = 1, rel_heights = c(1, 1.2),
          label_fontface = "plain"
        )
      },
      hpc = FALSE,
      trigger = trigger(condition = TRUE)
    ),

    #
    # Figure: Prevalence ----http://amantia:50639/graphics/plot_zoom_png?width=530&height=681
    #

    fig1_prevalence = target(
      {
        top_prevalences_plt <-
          plot_top_prevalent(prevalence_summary, lineages) +
          # Fix phylum and enviornment guides to appear twice
          guides(fill = FALSE, color = FALSE)

        habitat_plt <-
          "raw/media/envs.png" %>%
          magick::image_read() %>%
          magick::image_rotate(180) %>%
          magick::image_ggplot()

        collection_times_plt <-
          samples %>%
          # keep same plot order as habitat pictograms
          mutate(environment_group = environment_group %>% factor(levels = c("host", "soil", "aquatic"))) %>%
          ggplot(aes(collection_datetime, fill = environment_group)) +
          geom_histogram(position = "identity", binwidth = 60 * 60 * 24 * 60) +
          facet_wrap(~environment_group, ncol = 1, strip.position = "right") +
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
              eulerr::euler(shape = "circle") %>%
              plot(fills = environment_group_colors[environment_groups], main = list(label = .y, cex = 0.7), quantities = list(cex = 0.7), labels = FALSE)
            # ggplotify::as.ggplot()
          })) %>%
          ungroup() %>%
          select(subset_value, venn_plt) %>%
          deframe()

        prevalence_intersection_plt <- plot_prevalence_comparision(
          prevalence_summary,
          lineages,
          environment_groups,
          type = "Intersection"
        ) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
          guides(fill = guide_legend(ncol = 3))

        layout <- "
      AF
      BF
      CF
      DF
      EF
      "

        list(
          collection_times_plt + labs(tag = "A"),
          prevalence_intersection_plt + labs(tag = "B"),

          prevalences_venn_plts$Bacteria,
          prevalences_venn_plts$Fungi,
          plot_spacer(), # spacer to block floating
          top_prevalences_plt + theme(axis.text = element_text(size = 7)) + labs(tag = "C")
        ) %>%

          # patchwork
          wrap_plots() +
          plot_layout(
            design = layout,
            guides = "collect",
            heights = c(2, 1, 2, 2, 0.1, 1),
            widths = c(1, 3)
          )
      },
      trigger = trigger(change = habitat_colors)
    ),

    #
    # Figure: Diversity ----
    #

    fig_alpha_diversity = target(
      {
        collection_times_plt <-
          samples %>%
          # keep same plot order as habitat pictograms
          mutate(environment_group = environment_group %>% factor(levels = c("host", "soil", "aquatic"))) %>%
          ggplot(aes(collection_datetime, fill = environment_group)) +
          geom_histogram(position = "identity", binwidth = 60 * 60 * 24 * 60) +
          facet_wrap(~environment_group, ncol = 1, strip.position = "right") +
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

        geocoded_bioprojects <- read_csv("/analysis/raw/meta/geocoded_bioprojects.csv")

        worldmap_tbl <-
          samples %>%
          filter(!is.na(lat)) %>%
          bind_rows(
            samples %>% filter(is.na(lat)) %>% select(-lat, -lon) %>% left_join(geocoded_bioprojects)
          ) %>%
          filter(
            !is.na(environment_group) &
              sample_id %in% abundant_samples
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

        sampling_plt <-
          wrap_plots(
            collection_times_plt + theme(legend.position = "none") + labs(tag = "A"), worldmap_plt + labs(tag = "B"),
            nrow = 1, widths = c(1, 1.5)
          )

        alphadiv_plt <-
          plot_alphadiv(alphadiv, alphadiv_metrics, samples = samples) &
            guides(color = FALSE, fill = FALSE)
        alphadiv_plt$patches[[1]][[1]] <- alphadiv_plt$patches[[1]][[1]] + labs(tag = "C", y = "Alpha Diversity")

        alphadiv_cor_plt <-
          plot_alphadiv_cor(alphadiv, alphadiv_metrics, samples = samples) +
            guides(color = FALSE) &
            theme_pub()
        alphadiv_cor_plt$patches[[1]][[1]] <- alphadiv_cor_plt$patches[[1]][[1]] + labs(tag = "D")

        wrap_plots(sampling_plt, alphadiv_plt, alphadiv_cor_plt, ncol = 1, heights = c(0.6, 1, 0.8), guides = "collect")
      },
      hpc = FALSE
    ),

    fig_beta_diversity = target(
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
          ggplot(aes(reorder(taxon, R2), R2, fill = type)) +
          geom_bar(stat = "identity", position = "dodge") +
          geom_text(aes(label = label), position = position_dodge(width = 1), hjust = -0.2, vjust = 0.5, size = 2.5) +

          # # dummy geom to create legend for short axis labels
          # geom_text(aes(label = taxon_short), size = 0) + # dummy to create legend
          # scale_discrete_identity(
          #     aesthetics = "label",
          #     name = "Explaining genus",
          #     breaks = tbl %>% distinct(taxon, taxon_short) %>% pull(taxon_short),
          #     labels = tbl %>% distinct(taxon, taxon_short) %>% pull(taxon),
          #     guide = "legend"
          # ) +
          # guides(label = guide_legend(nrow = 7, override.aes = list(size = 2.5))) +

          coord_flip() +
          guides(fill = guide_legend(reverse = TRUE)) +
          facet_wrap(~kingdom, nrow = 1, scales = "free") +
          scale_fill_modeltype() +
          scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
          theme_pub() +
          theme(axis.text.y = element_text(face = "italic")) +
          labs(
            tag = "B",
            x = "Explaining genus",
            y = "partial R² (%)",
            fill = "Model type"
          )

        stepwise_constrained_ordinations_plots$plt[[1]] <-
          stepwise_constrained_ordinations_plots$plt[[1]] +
          labs(tag = "A")

        list(
          wrap_plots(
            stepwise_constrained_ordinations_plots %>%
              filter(kingdom == "Bacteria") %>%
              pull(plt) %>%
              first() +
              labs(color = "Explaining phylum") +
              guides() +
              theme_pub(),
            stepwise_constrained_ordinations_plots %>%
              filter(kingdom == "Fungi") %>%
              pull(plt) %>%
              first() +
              labs(color = "Explaining phylum") +
              theme_pub()
          ),
          stepwise_constrained_ordinations_contribution_plt
        ) %>%
          wrap_plots(ncol = 1) +
          plot_layout(guides = "collect", heights = c(1, 1)) +
          plot_annotation() &
          theme(
            legend.position = "right",
            legend.box = "vertical",

            # squeeze tag labels
            plot.tag = element_text(hjust = 1, vjust = -3),
            plot.margin = margin(t = -10)
          )
      },
      hpc = FALSE
    ),

    fig2_diversity = target(
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
          ggplot(aes(reorder(taxon, R2), R2, fill = type)) +
          geom_bar(stat = "identity", position = "dodge") +
          geom_text(aes(label = label), position = position_dodge(width = 1), hjust = -0.2, vjust = 0.5, size = 2.5) +

          # # dummy geom to create legend for short axis labels
          # geom_text(aes(label = taxon_short), size = 0) + # dummy to create legend
          # scale_discrete_identity(
          #     aesthetics = "label",
          #     name = "Explaining genus",
          #     breaks = tbl %>% distinct(taxon, taxon_short) %>% pull(taxon_short),
          #     labels = tbl %>% distinct(taxon, taxon_short) %>% pull(taxon),
          #     guide = "legend"
          # ) +
          # guides(label = guide_legend(nrow = 7, override.aes = list(size = 2.5))) +

          coord_flip() +
          facet_wrap(~kingdom, nrow = 1, scales = "free") +
          scale_fill_modeltype() +
          scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
          theme_pub() +
          theme(axis.text.y = element_text(face = "italic")) +
          labs(
            tag = "D",
            x = "Explaining genus",
            y = "partial R² (%)",
            fill = "Model type"
          )

        stepwise_constrained_ordinations_plots$plt[[1]] <-
          stepwise_constrained_ordinations_plots$plt[[1]] +
          labs(tag = "C")

        list(
          alphadiv_plt,
          alphadiv_cor_plt,
          wrap_plots(
            stepwise_constrained_ordinations_plots %>%
              filter(kingdom == "Bacteria") %>%
              pull(plt) %>%
              first() +
              labs(color = "Explaining phylum") +
              guides() +
              theme_pub(),
            stepwise_constrained_ordinations_plots %>%
              filter(kingdom == "Fungi") %>%
              pull(plt) %>%
              first() +
              labs(color = "Explaining phylum") +
              theme_pub()
          ),
          stepwise_constrained_ordinations_contribution_plt
        ) %>%
          wrap_plots(ncol = 1) +
          plot_layout(guides = "collect", heights = c(2, 1.2, 1, 1)) +
          plot_annotation() &
          theme(
            legend.position = "right",
            legend.box = "vertical",

            # squeeze tag labels
            plot.tag = element_text(hjust = 1, vjust = -3),
            plot.margin = margin(t = -10)
          )
      },
      trigger = trigger(change = habitat_colors),
      hpc = FALSE
    ),


    #
    # Figure: Abundance ----
    #

    fig3_abundance = target(
      {
        tbl <-
          sub_abundances %>%
          filter(subset_name == "environment_group,kingdom" & norm_method == "tss") %>%
          discard(~ .x %>%
            unique() %>%
            length() < 2) %>%
          unnest(data) %>%
          inner_join(prevalences %>% unnest(data), by = c("subset_value", "taxon"))

        prevalent_taxa <-
          tbl %>%
          filter(prevalence_perc > min_heatmap_prevalence_perc) %>%
          pull(taxon) %>%
          unique()

        filtered_tbl <-
          tbl %>%
          filter(
            !is.na(dominance) & !is.na(environment_group) & taxon %in% prevalent_taxa
          ) %>%
          group_by(taxon) %>%
          mutate(
            abundance = abundance %>% log10() %>% scale() %>% pluck(1)
          ) %>%
          ungroup() %>%
          arrange(environment_group, habitat, bioproject_id)

        mat <-
          filtered_tbl %>%
          select(sample_id, taxon, abundance) %>%
          pivot_wider(names_from = sample_id, values_from = abundance, values_fill = list(abundance = -100)) %>%
          as_matrix("taxon")

        samples_tbl <-
          mat %>%
          colnames() %>%
          tibble(sample_id = .) %>%
          left_join(samples) %>%
          mutate(
            habitat = habitat %>% factor(levels = names(habitat_colors))
          )

        taxa_tbl <-
          mat %>%
          rownames() %>%
          tibble(taxon = .) %>%
          left_join(lineages) %>%
          left_join(
            dominances %>%
              select(taxon, environment_group, dominance) %>%
              pivot_wider(names_from = environment_group, values_from = dominance, names_prefix = "dominance_")
          ) %>%
          mutate(
            # set other phyla to NA
            phylum = ifelse(phylum %in% names(phyla_colors), phylum, NA)
          ) %>%
          mutate(
            phylum = phylum %>% factor(levels = names(phyla_colors))
          )

        abundance_span <- filtered_tbl$abundance %>%
          map_dbl(abs) %>%
          max() %>%
          round()

        abundance_plt <-
          Heatmap(
            matrix = mat,
            name = "z(log10(Abundance))",

            row_split = taxa_tbl %>% select(kingdom),
            # column_split = samples_tbl %>% select(environment_group),

            cluster_rows = TRUE,
            cluster_columns = FALSE,

            # UPGMA clustering
            clustering_method_rows = "average",
            clustering_method_columns = "average",

            # clustering_distance_rows = "manhattan",
            # clustering_distance_columns = "manhattan",

            top_annotation = HeatmapAnnotation(
              which = "column",

              Environment = samples_tbl$environment_group,
              Habitat = samples_tbl$habitat,
              # `Dysbalance Bacteria` = samples_tbl$dysbalance_Bacteria,
              # `Dysbalance Fungi` = samples_tbl$dysbalance_Fungi,


              col = list(
                Environment = environment_group_colors %>% simplify(),
                Habitat = habitat_colors %>% simplify()
                # `Dysbalance Bacteria` = dysbalance_colors %>% simplify(),
                # `Dysbalance Fungi` = dysbalance_colors %>% simplify()
              ),

              annotation_legend_param = list(
                # `Dysbalance Bacteria` = list(title = "Dysbalance")
                Habitat = list(ncol = 2)
              ),

              show_legend = c(TRUE, TRUE)
            ),

            left_annotation = HeatmapAnnotation(
              which = "row",

              Kingdom = taxa_tbl$kingdom,
              Phylum = taxa_tbl$phylum,
              # `Dominance aquatic` = taxa_tbl$dominance_aquatic,
              # `Dominance host` = taxa_tbl$dominance_host,
              # `Dominance soil` = taxa_tbl$dominance_soil,

              col = list(
                Kingdom = kingdoms_colors %>% simplify(),
                Phylum = phyla_colors %>% simplify()

                # `Dominance aquatic` = dominance_colors %>% simplify(),
                # `Dominance host` = dominance_colors %>% simplify(),
                # `Dominance soil` = dominance_colors %>% simplify()
              ),

              # annotation_legend_param = list(
              #   `Dominance aquatic` = list(title = "Dominance")
              # ),

              show_legend = c(TRUE, TRUE)
            ),

            heatmap_legend_param = list(ncol = 1),

            col = circlize::colorRamp2(c(-abundance_span, 0, abundance_span), viridis::viridis(3)),

            show_column_names = FALSE,
            show_row_names = FALSE,
            show_row_dend = TRUE,
            show_column_dend = TRUE
          ) %>%
          draw(heatmap_legend_side = "bottom", annotation_legend_side = "bottom")

        abundance_plt
      },
      # trigger on change of non-targets
      trigger = trigger(change = habitat_colors)
    ),

    #
    # Figure: Coabundance ----
    #

    fig4_coabundance = {
      graph_plt <-
        environment_group_multi_graph_merged_plt +
        labs(title = NULL, subtitle = NULL)

      specific_coabundances_plt <-
        environment_group_specific_coabundances_plt + labs(title = NULL, subtitle = NULL)

      layout <- "
      AB
      CB
      DB
      "

      list(
        environment_group_coabundance_venn_plt + theme_pub(show_axes = FALSE),
        graph_plt + theme_pub(show_axes = FALSE),
        specific_coabundances_plt + theme_pub(),
        topology_plt + theme_pub() + theme(axis.text.x = element_blank())
      ) %>%
        wrap_plots(tag_level = "new") +
        plot_layout(design = layout, guides = "collect", widths = c(1, 2), heights = c(1, 2, 2)) +
        plot_annotation(tag_levels = "A")
    }
  )
}
