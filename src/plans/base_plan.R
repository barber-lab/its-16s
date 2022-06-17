get_base_plan <- function() {
  drake_plan(
    #
    # Params ----
    #

    # Must be flat and separated to not always trigger cache invalidation
    # on every target
    n_batches = 32,
    norm_methods = c("raw", "tss"),
    pooling_col = "genus",
    min_prevalence_abundance = 0,
    taxranks = get_taxranks(pooling_col, direction = "with_upstream"),

    taxa_groupings = c("kingdom", "dominance"),

    dominance_quantile = 0.1,
    dominance_samples_grouping = "environment_group",
    dominance_taxa_grouping = "kingdom",

    alphadiv_metrics = c("Chao1", "Shannon"),
    dist_types = c("samples", "taxa"),
    dist_methods = "bray",
    dist_norm_methods = "tss",

    dysbalance_groupings = c("bioproject_id", "kingdom"),
    dysbalance_dist_method = "bray",
    max_balanced_quantile = 0.8,

    adonis_dist_method = "bray",

    groupings = colnames(subsets),

    focus_genera = {
      prev_abundance_taxa <- c(
        "Flavobacterium", "hgcI clade", "Pseudomonas", "CL500-29 marine group",
        "Rhodoferax", "Limnohabitans", "Polynucleobacter", "OM43 clade",
        "Arenimonas", "Polaromonas", "Malassezia", "Cortinarius", "Aureobasidium",
        "Diplodia", "Crassiclypeus", "Peniophora", "Acidea", "Phaeoacremonium",
        "Phylloporus", "Cytospora", "Streptococcus", "Escherichia-Shigella",
        "Lactobacillus", "Lachnospiraceae UCG-008", "Roseburia", "Veillonella",
        "Blautia", "Bacteroides", "Ruminococcaceae UCG-014", "Candida",
        "Cladosporium", "Penicillium", "Rhodotorula", "Aspergillus",
        "Trichosporon", "Mycosphaerella", "Alternaria", "Saccharomyces",
        "Fusarium", "Bradyrhizobium", "Bacillus", "Sphingomonas", "Bryobacter",
        "Streptomyces", "Haliangium", "Acidibacter", "Rhodoplanes", "Pseudolabrys",
        "Mortierella", "Trichoderma", "Chaetomium", "Cladophialophora",
        "Pseudogymnoascus", "Chloridium"
      ) %>% unique()

      dbrda_taxa <- c(
        # fungi explaining bacterial dissimilarity
        "Mortierella", "Trichocladium", "Candida", "Trichoderma",
        "Russula", "Oidiodendron", "Penicillium",
        "Malassezia", "Chaetomium", "Podospora",

        # bacteria explaining fungal dissimilarity
        "Conexibacter", "Bacillus", "Lysobacter", "Burkholderia-Caballeronia-Paraburkholderia",
        "Escherichia-Shigella", "Candidatus Solibacter", "Acidovorax", "Pseudolabrys",
        "Legionella", "Gemmatimonas",

        # bacteria explaining fungal dissimilarity in one env
        "Comamonas", "Legionella", "Sphingorhabdus", "Aquabacterium",
        "Flavobacterium", "Escherichia-Shigella", "Staphylococcus", "Bacteroides",
        "Roseburia", "Rothia", "Jatrophihabitans", "Lysobacter", "Pseudolabrys",
        "Candidatus Koribacter", "Gemmatimonas",

        # fungal explaining bacterial dissimilarity in one env
        "Cortinarius", "Hebeloma", "Acidea", "Peniophora", "Oidiodendron", "Aspergillus",
        "Chaetomium", "Candida", "Saccharomyces", "Penicillium", "Alternaria", "Trichocladium",
        "Fusarium"
      ) %>% unique()


      text_genera <- c(
        "Mortierella", "Trichocladium", "Candida", "Conexibacter",
        "Bacillus", "Escherichia-Shigella", "Lysobacter", "Penicillium",
        "Aspergillus", "Chaetomium", "Alternaria", "Hebeloma"
      ) %>% unique()

      focus_genera <-
        list(prev_abundance_taxa, dbrda_taxa, text_genera) %>%
        reduce(union) %>%
        unique()

      focus_genera
    },

    #
    # Subsets -----
    #

    environment_groups = c("host", "aquatic", "soil"),
    kingdoms_groups = c("Bacteria", "Fungi"),
    dominance_groups = c("dominant", "rare"),
    dysbalance_groups = c("dysbalanced", "balanced"),
    # manual curation for prevalent habitats
    habitat_groups = c(
      "temperate forests",
      "conifer forests",
      "mediterranean forests",
      "other hosts",
      "insects",
      "other freshwater",
      "boreal forests",
      "lung",
      "tropical forests",
      "gut",
      "other aquatic",
      "large lakes",
      "plants",
      "marine",
      "skin",
      "mouth"
    ),
    bioproject_groups = bioprojects$bioproject_id,

    all_subsets = tibble(subset_name = "all", subset_value = "all", subset_query = "TRUE"),

    kingdom_subsets = {
      tibble(kingdom = kingdoms_groups) %>%
        transmute(
          subset_name = "kingdom",
          subset_value = str_glue("{kingdom}"),
          subset_query = str_glue("kingdom == '{kingdom}'")
        )
    },

    habitat_subsets = {
      tibble(habitat = habitat_groups) %>%
        transmute(
          subset_name = "habitat",
          subset_value = str_glue("{habitat}"),
          subset_query = str_glue("habitat == '{habitat}'")
        )
    },

    bioproject_subsets = {
      tibble(bioproject_id = bioproject_groups) %>%
        transmute(
          subset_name = "bioproject_id",
          subset_value = str_glue("{bioproject_id}"),
          subset_query = str_glue("bioproject_id == '{bioproject_id}'")
        )
    },

    environment_group_subsets = {
      tibble(environment_group = environment_groups) %>%
        transmute(
          subset_name = "environment_group",
          subset_value = str_glue("{environment_group}"),
          subset_query = str_glue("environment_group == '{environment_group}'")
        )
    },

    environment_group_kingdom_subsets = {
      tibble(environment_group = environment_groups) %>%
        expand_grid(kingdom = kingdoms_groups) %>%
        transmute(
          subset_name = "environment_group,kingdom",
          subset_value = str_glue("{environment_group},{kingdom}"),
          subset_query = str_glue("environment_group == '{environment_group}' & kingdom == '{kingdom}'")
        )
    },

    environment_group_bioproject_kingdom_subsets = {
      samples %>%
        distinct(environment_group, bioproject_id) %>%
        expand_grid(kingdom = kingdoms_groups) %>%
        transmute(
          subset_name = "environment_group,bioproject_id,kingdom",
          subset_value = str_glue("{environment_group},{bioproject_id},{kingdom}"),
          subset_query = str_glue("environment_group == '{environment_group}' & bioproject_id == '{bioproject_id}' & kingdom == '{kingdom}'")
        )
    },

    habitat_kingdom_subsets = {
      tibble(habitat = habitat_groups) %>%
        expand_grid(kingdom = kingdoms_groups) %>%
        transmute(
          subset_name = "habitat,kingdom",
          subset_value = str_glue("{habitat},{kingdom}"),
          subset_query = str_glue("habitat == '{habitat}' & kingdom == '{kingdom}'")
        )
    },

    dysbalanced_kingdom_subsets = {
      expand_grid(dysbalance = dysbalance_groups) %>%
        expand_grid(kingdom = kingdoms_groups) %>%
        transmute(
          subset_name = "dysbalance,kingdom",
          subset_value = str_glue("{dysbalance},{kingdom}"),
          subset_query = str_glue("kingdom == '{kingdom}' & dysbalance == '{dysbalance}'")
        )
    },

    dysbalanced_environment_group_kingdom_subsets = {
      tibble(environment_group = environment_groups) %>%
        expand_grid(dysbalance = dysbalance_groups) %>%
        expand_grid(kingdom = kingdoms_groups) %>%
        transmute(
          subset_name = "dysbalance,environment_group,kingdom",
          subset_value = str_glue("{dysbalance},{environment_group},{kingdom}"),
          subset_query = str_glue("environment_group == '{environment_group}' & kingdom == '{kingdom}' & dysbalance == '{dysbalance}'")
        )
    },

    dominance_environment_group_kingdom_subsets = {
      tibble(environment_group = environment_groups) %>%
        expand_grid(dominance = dominance_groups) %>%
        expand_grid(kingdom = kingdoms_groups) %>%
        transmute(
          subset_name = "dominance,environment_group,kingdom",
          subset_value = str_glue("{dominance},{environment_group},{kingdom}"),
          subset_query = str_glue("dominance == '{dominance}' & environment_group == '{environment_group}' & kingdom == '{kingdom}'")
        )
    },

    subsets = {
      list(
        all_subsets,
        kingdom_subsets,
        environment_group_subsets,
        environment_group_kingdom_subsets,
        # environment_group_bioproject_kingdom_subsets,
        habitat_subsets,
        habitat_kingdom_subsets,
        bioproject_subsets
        # TODO: Add more
      ) %>%
        # ignore type glue
        map(~ .x %>% mutate_all(as.character)) %>%
        bind_rows() %>%
        view()
    },

    #
    # Load raw data ----
    #
    bioprojects = file_in("raw/meta/bioprojects.csv") %>% read_csv(),
    bioproject_annotations = file_in("raw/meta/bioproject_annotations.csv") %>% read_csv(),
    biosamples = file_in("raw/meta/biosamples.csv") %>% read_csv(),
    samples = {
      file_in("raw/meta/samples.csv") %>%
        read_csv() %>%
        left_join(biosamples) %>%
        mutate(
          biosample_id,
          location = attr_lat_lon %>% map(parse_lon_lat),
          date = attr_collection_date %>% lubridate::as_date(),
          collection_datetime = collection_datetime %>% as_datetime()
        ) %>%
        unnest_wider(location) %>%
        add_habitat() %>%
        filter(sample_id %in% abundant_samples)
    },

    lineages = {
      file_in("raw/meta/lineages.csv") %>%
        read_csv() %>%
        dplyr::group_by_at(pooling_col) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup() %>%
        dplyr::select_at(get_taxranks(pooling_col, direction = "with_upstream")) %>%
        dplyr::distinct() %>%
        dplyr::mutate(taxon = .[[pooling_col]])
    },

    traits = {
      traits_funguild <-
        file_in("raw/databases/FUNGuild.csv") %>%
        read_csv() %>%
        transmute(
          taxon = genus,
          trophicMode,
          guild
        )

      traits_faprotax <-
        file_in("raw/databases/faprotax.csv") %>%
        read_csv() %>%
        transmute(
          taxon = phylogeny %>% str_extract("g__[^;]+;") %>% str_remove("^g__") %>% str_remove(";$"),
          faprotax_entry = entry
        ) %>%
        group_by(taxon) %>%
        slice(1) %>%
        ungroup()

      traits_funguild %>%
        full_join(traits_faprotax)
    },

    unknown_habitat = "unknown",
    habitats = {
      samples %>%
        distinct(habitat, environment_group) %>%
        arrange(environment_group, habitat) %>%
        add_row(habitat = unknown_habitat, environment_group = unknown_habitat)
    },

    taxa_grouping_pairs = {
      taxa_groupings %>%
        map(~ .x %>%
          discard(~ .x == "all") %>%
          head(2)) %>%
        enframe() %>%
        unnest_wider(value) %>%
        set_colnames(c("taxa_grouping", "taxa_group_x", "taxa_group_y"))
    },

    abundances = {
      norm_methods %>%
        map(
          ~ get_abundance_tbl(
            taxrank = pooling_col,
            normalization_method = .x,
            filter_zeros = TRUE
          )
        ) %>%
        enframe(name = "norm_method", value = "data") %>%
        mutate(norm_method = norm_methods) %>%
        arrange_columns()
    },

    #
    # Misc ------
    #

    mapped_samples = {
      samples %>%
        filter(!is.na(environment_group))
    },

    abundant_samples = {
      abundances %>%
        filter(norm_method == "tss") %>%
        pluck("data", 1) %>%
        filter(abundance >= min_prevalence_abundance) %>%
        pull(sample_id) %>%
        unique()
    },

    failed_traces = get_failed_traces(),

    sub_dysbalances = target(
      {
        if (
          # This is of length 0 sometimes ...
          length(distances$dist_type) == 0 ||
            distances$dist_type[[1]] != "samples"
        ) {
          return(tibble())
        }

        res <- possibly(get_dysbalance, NULL)(
          distance = distances$dist[[1]],
          max_balanced_quantile = max_balanced_quantile,
          max_clusters = max_clusters
        )

        if (is.null(res)) {
          return(tibble())
        }

        distances %>%
          discard(is.list) %>%
          bind_cols(
            res %>%
              enframe() %>%
              pivot_wider(names_from = name, values_from = value)
          ) %>%
          arrange_columns()
      },
      dynamic = map(distances),
      hpc = TRUE
    ),

    #
    # Augment ----
    #

    dominances = get_dominances(
      abundances = abundances,
      samples = samples,
      lineages = lineages,
      pooling_col = pooling_col,
      samples_grouping = dominance_samples_grouping,
      taxa_grouping = dominance_taxa_grouping,
      dominance_quantile = dominance_quantile
    ),

    dysbalances = get_dysbalances(
      abundances = abundances,
      samples = samples,
      lineages = lineages,
      pooling_col = pooling_col,
      groupings = dysbalance_groupings,
      max_balanced_quantile = max_balanced_quantile,
      max_clusters = max_clusters,
      dist_method = dysbalance_dist_method
    ),

    augmented_abundances = target(
      {
        raw_abugmented_abundances <-
          abundances %>%
          mutate(
            data = data %>% map(~ {
              .x %>%
                left_join(
                  samples %>% select(sample_id, environment_group, habitat),
                  by = "sample_id"
                ) %>%
                left_join(lineages) %>%
                left_join(dysbalances) %>%
                left_join(dominances) %>%
                arrange_columns()
            })
          )

        #
        # sanity checks
        #

        both_kingdoms_present_samples <-
          raw_abugmented_abundances$data[[1]] %>%
          filter(abundance > 0) %>%
          distinct(sample_id, kingdom) %>%
          mutate(is_present = TRUE) %>%
          pivot_wider(names_from = kingdom, values_from = is_present) %>%
          filter(Bacteria & Fungi) %>%
          pull(sample_id) %>%
          unique()

        fungal_and_bacterial_taxa <-
          raw_abugmented_abundances$data[[1]] %>%
          filter(kingdom %in% c("Bacteria", "Fungi")) %>%
          pull(taxon) %>%
          unique()

        # intersect all criteria
        sanity_passed_samples <- both_kingdoms_present_samples
        sanity_passed_taxa <- fungal_and_bacterial_taxa

        raw_abugmented_abundances %>%
          mutate(
            data = data %>% map(~ filter(
              .x,
              sample_id %in% sanity_passed_samples & taxon %in% sanity_passed_taxa
            ))
          )
      },
      dynamic = map(abundances),
      hpc = TRUE
    ),

    #
    # Subset abundances ----
    #
    sub_abundances = target(
      {
        subsets %>%
          expand_grid(augmented_abundances) %>%
          mutate(data = data %>% map2(subset_query, ~ .x %>% filter(eval(parse(text = .y))))) %>%
          select(-subset_query) %>%
          # sanity checks
          filter(data %>% map_lgl(~ nrow(.x) > 0))
      },
      dynamic = map(subsets),
      hpc = TRUE
    ),

    #
    # Analyze abundance ----
    #
    prevalences = {
      sub_abundances %>%
        filter(norm_method == "tss") %>%
        select(-norm_method) %>%
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
              filter(abundance >= min_prevalence_abundance) %>%
              summarise(
                prevalence_n = n(),
                prevalence_perc = prevalence_n / n_samples * 100
              )
          })
        )
    },

    prevalence_summary = target(
      summarise_prevalences(prevalences = prevalences),
    ),

    # Influence of filter criteria for abundance and prevalence
    # on taxa counts
    prevalence_abundance_thresholds_plt = {
      unnested_sub_abundances <-
        sub_abundances %>%
        filter(subset_name == "kingdom" & norm_method == "tss") %>%
        unnest(data)

      # zero must be included for reference normalization
      min_abundances <- seq(0, 100, length.out = 11)
      min_prevalences <- seq(0, 50, length.out = 11)

      abundance_prevalence_tbl <-
        expand_grid(
          min_abundance_perc = min_abundances,
          min_prevalence_perc = min_prevalences
        ) %>%
        mutate(
          data = min_abundance_perc %>% map2(min_prevalence_perc, ~ {
            n_samples <- unnested_sub_abundances$sample_id %>%
              unique() %>%
              length()

            unnested_sub_abundances %>%
              group_by(environment_group, kingdom, taxon) %>%
              mutate(
                prevalence_n = n(),
                prevalence_perc = prevalence_n / n_samples * 100
              ) %>%
              filter(abundance * 100 > .x & prevalence_perc > .y) %>%
              group_by(environment_group, kingdom) %>%
              distinct(taxon) %>%
              count()
          })
        ) %>%
        unnest(data) %>%
        # sanity check
        filter(!is.na(environment_group)) %>%
        complete(
          min_abundance_perc = min_abundances,
          min_prevalence_perc = min_prevalences,
          environment_group = environment_groups,
          kingdom = kingdoms_groups,
          fill = list(n = 0)
        )

      # label in which soil ist
      abundance_prevalence_tbl %>%
        # normalise by percentag of taxa kept
        group_by(environment_group, kingdom) %>%
        mutate(n = n / max(n) * 100) %>%
        ungroup() %>%
        ggplot(aes(min_prevalence_perc, min_abundance_perc)) +
        geom_tile(aes(fill = log(n))) +
        geom_text(
          data = {
            # Area in parameter space in which there are moe soil taxa than host ones
            abundance_prevalence_tbl %>%
              group_by(min_abundance_perc, min_prevalence_perc, kingdom) %>%
              nest() %>%
              transmute(most_taxa_environment_group = data %>% map_chr(~ {
                .x %>%
                  arrange(-n) %>%
                  pull(environment_group) %>%
                  pluck(1)
              }, )) %>%
              mutate(label = most_taxa_environment_group %>% str_sub(1, 1))
          },
          mapping = aes(label = label, color = most_taxa_environment_group),
          size = 7
        ) +
        facet_grid(kingdom ~ environment_group, scales = "free") +
        scale_color_environment_group() +
        # log(0) = NA color
        scale_fill_viridis_c(na.value = viridis::viridis(1)) +
        theme(panel.grid.major = element_line()) +
        labs(
          title = "Taxa counts are highly influenced by filter critera",
          x = "Min. prevalence (%)",
          y = "Min. abundance (%)",
          fill = "log10(% genera kept)",
          color = "Environment with most taxa"
        )
    },

    fig_abundance_prevalence_cor = {
      bind_rows(
        prevalences %>%
          filter(subset_name == "environment_group,kingdom") %>%
          unnest(data) %>%
          separate(col = subset_value, into = c("environment_group", "kingdom"), sep = ",") %>%
          select(environment_group, kingdom, taxon, value = prevalence_perc) %>%
          mutate(value = value %>% rank() %>% scale_min_max() * 100) %>%
          complete(taxon, environment_group, fill = list(value = 0)) %>%
          mutate(type = "Prevalence (percentile)"),

        abundances %>%
          filter(norm_method == "raw") %>%
          pull(data) %>%
          first() %>%
          group_by(sample_id) %>%
          mutate(abundance = abundance %>% rank() %>% scale_min_max() * 100) %>%
          left_join(samples %>% select(sample_id, environment_group)) %>%
          left_join(lineages %>% select(taxon, kingdom)) %>%
          filter(!is.na(environment_group)) %>%
          group_by(kingdom, taxon, environment_group) %>%
          summarise(abundance = mean(abundance)) %>%
          ungroup() %>%
          select(kingdom, taxon, environment_group, value = abundance) %>%
          complete(taxon, environment_group, fill = list(value = 0)) %>%
          mutate(type = "Abundance (percentile)")
      ) %>%
        pivot_wider(names_from = type, values_from = value) %>%
        left_join(lineages) %>%
        filter(kingdom %in% kingdoms_groups) %>%
        ggplot(aes(`Prevalence (percentile)`, `Abundance (percentile)`)) +
        geom_point(aes(color = environment_group), alpha = 0.3) +
        facet_grid(environment_group ~ kingdom) +
        stat_cor(method = "pearson") +
        coord_fixed() +
        scale_color_environment_group()
    },

    fig_rel_core_metagenome = target({
      prevalences %>%
        filter(subset_name == "environment_group") %>%
        unnest(data) %>%
        transmute(
          environment_group = subset_value,
          taxon,
          is_core = prevalence_perc > 50
        ) %>%
        left_join(lineages) %>%
        mutate(
          phylum = phylum %>% map2_chr(kingdom, phylum_to_color_name) %>% factor(level = names(phyla_colors)),
          core = is_core %>% ifelse("core", "periphery")
        ) %>%
        count(environment_group, kingdom, phylum, core) %>%
        group_by(environment_group, kingdom) %>%
        mutate(n = n / sum(n) * 100) %>%
        mutate(x = paste0(kingdom, " ", core)) %>%
        ggplot(aes(x, n, fill = phylum)) +
        geom_bar(stat = "identity") +
        facet_wrap(~environment_group, ncol = 1, strip.position = "right") +
        scale_fill_phyla() +
        coord_flip() +
        scale_y_continuous(expand = c(0, 0)) +
        theme(panel.grid.major = element_blank()) +
        labs(x = "", y = "Genera (% per kingdom and environment)", fill = "Phylum")
    }),

    fig_venn_core_metagenome = {
      plt <-
        prevalences %>%
        filter(subset_name == "environment_group") %>%
        unnest(data) %>%
        transmute(
          environment_group = subset_value,
          taxon,
          core_group = ifelse(prevalence_perc > 50, "core", "periphery")
        ) %>%
        left_join(lineages %>% select(taxon, kingdom)) %>%
        nest(-c(core_group, environment_group, kingdom)) %>%
        mutate(data = data %>% map(~ .x$taxon)) %>%
        nest(-c(kingdom, core_group)) %>%
        mutate(
          plt = list(data, kingdom, core_group) %>% pmap(~ {
            .x %>%
              deframe() %>%
              ggeulerr() +
              facet_grid(~ str_glue("{.y}  {..3}")) +
              scale_fill_environment_group(drop = FALSE)
          })
        ) %>%
        arrange(core_group, kingdom) %>%
        pull(plt) %>%
        wrap_plots(nrow = 2, guides = "collect")
      # Manually remove redundant but different legend
      plt$patches$plots[[2]] <- plt$patches$plots[[2]] + theme(legend.position = "none")
      plt
    },

    core_metagenome_plt = target(
      {
        prevalences %>%
          filter(subset_name == "environment_group") %>%
          unnest(data) %>%
          transmute(
            environment_group = subset_value,
            taxon,
            is_core = prevalence_perc > 50
          ) %>%
          left_join(lineages) %>%
          mutate(
            phylum = phylum %>% map2_chr(kingdom, phylum_to_color_name) %>% factor(level = names(phyla_colors)),
            core = is_core %>% ifelse("core", "periphery")
          ) %>%
          count(environment_group, phylum, core) %>%
          ggplot(aes(core, n, fill = phylum)) +
          geom_bar(stat = "identity") +
          facet_wrap(~environment_group, ncol = 1, strip.position = "right") +
          scale_fill_phyla() +
          coord_flip() +
          scale_y_continuous(expand = c(0, 0)) +
          theme(panel.grid.major = element_blank()) +
          labs(x = "", y = "Genera", fill = "Phylum")
      },
      trigger = trigger(condition = TRUE)
    ),

    prevalence_abundance_heatmap_plt = {
      sub_abundances %>%
        rename(abundance = data) %>%
        inner_join(prevalences %>% rename(prevalence = data)) %>%
        filter(subset_name == "environment_group,kingdom" & norm_method == "tss") %>%
        head() %>%
        transmute(data = abundance %>% map2(prevalence, ~ {
          .x %>% inner_join(.y, by = "taxon")
        })) %>%
        unnest(data) %>%
        group_by(environment_group, kingdom, taxon) %>%
        # get one value per taxon and environment group
        summarise(
          # prevalence is the same for all samples
          prevalence_perc = first(prevalence_perc),
          abundance = median(abundance)
        ) %>%
        ggplot(aes(prevalence_perc, abundance)) +
        stat_density_2d(
          geom = "raster",
          mapping = aes(fill = after_stat(density)),
          contour = FALSE,
          n = 500
        ) +
        scale_y_log10(expand = c(0, 0)) +
        scale_x_log10(expand = c(0, 0)) +
        facet_grid(environment_group ~ kingdom) +
        labs(
          x = "Prevalence (% samples)",
          y = "Median abundance (%)",
          fill = "Taxa frequency"
        )
    },

    #
    # Analyze alpha diversity ----
    #
    alphadiv = target(
      {
        data <-
          subsets %>%
          inner_join(
            sub_abundances %>% filter(norm_method == "raw")
          ) %>%
          head(1)

        res <- possibly(get_alpha_div_tbl, NA)(data$data[[1]], pooling_col = "taxon")

        data %>%
          select(-subset_query, -norm_method) %>%
          select(-where(is.list)) %>%
          mutate(data = list(res)) %>%
          arrange_columns()
      },
      dynamic = map(subsets),
      hpc = TRUE
    ),

    alphadiv_summary = target(
      summarise_alphadiv(alphadiv = alphadiv, alphadiv_metrics = alphadiv_metrics),
      dynamic = group(alphadiv),
      hpc = FALSE
    ),

    alphadiv_by_taxa_group = target(
      {
        get_alphadiv_plot_by_taxa_group(
          alphadiv = alphadiv,
          taxa_grouping = taxa_grouping_pairs$taxa_grouping[[1]],
          taxa_group_x = taxa_grouping_pairs$taxa_group_x[[1]],
          taxa_group_y = taxa_grouping_pairs$taxa_group_y[[1]],
          metric = alphadiv_metrics,
          cor_method = "spearman"
        ) %>%
          mutate(metric = alphadiv_metrics) %>%
          arrange_columns()
      },
      dynamic = cross(taxa_grouping_pairs, alphadiv_metrics),
      format = "rds", # qs leads to
      hpc = FALSE
    ),

    alphadiv_lmm_tests = {
      alphadiv %>%
        filter(subset_name == "kingdom") %>%
        expand_grid(
          alphadiv_metric = alphadiv_metrics,
          sample_group = c("environment_group", "habitat")
        ) %>%
        mutate(
          test = list(data, sample_group, alphadiv_metric) %>% pmap(~ {
            # kruskal.test(
            #   str_glue("{..3} ~ {..2}") %>% as.formula(),
            #   data = ..1 %>% left_join(samples)
            # ) %>%
            #   tidy()

            lmerTest::lmer(
              formula = str_glue("{..3} ~ {..2} + (1|bioproject_id)") %>% as.formula(),
              data = ..1 %>% left_join(samples)
            ) %>%
              anova() %>%
              as_tibble()
          })
        ) %>%
        unnest(test)
    },

    alphadiv_lmm_plt = {
      alphadiv_pre_tbl <-
        alphadiv %>%
        filter(subset_name == "kingdom") %>%
        group_by(subset_value) %>%
        slice(1) %>%
        ungroup()

      alphadiv_model_tbl <-
        alphadiv_pre_tbl %>%
        expand_grid(
          samples_grouping = c("habitat", "environment_group"),
          alphadiv_metric = c("Chao1", "Shannon")
        ) %>%
        mutate(
          data = data %>% map(~ .x %>% inner_join(samples)),

          model = list(data, alphadiv_metric, samples_grouping) %>% pmap(~ {
            str_glue("{..2} ~ 0 + {..3} + (1|bioproject_id)") %>%
              as.formula() %>%
              lme4::lmer(data = ..1) %>%
              list()
          }),

          parameters = model %>% map2(samples_grouping, ~ {
            .x %>%
              first() %>%
              # Adjust for multiple testing and also standardize factors
              parameters::parameters(p_adjust = "fdr", standardize = "basic") %>%
              as_tibble()
          }),
        ) %>%
        transmute(
          alphadiv_metric,
          subset_value,
          samples_grouping,
          parameters
        ) %>%
        unnest(parameters) %>%
        mutate(
          predicted = CI_low + (CI_high - CI_low) / 2,
          samples_grouping = Parameter %>% str_extract("^(habitat|environment_group)") %>% dplyr::recode("environment_group" = "Environment"),
          samples_group = Parameter %>%
            str_remove("^(habitat|environment_group)") %>%
            factor(levels = names(environment_group_colors) %>%
              c(names(habitat_colors)))
        )

      alphadiv_model_tbl %>%
        filter(alphadiv_metric == "Shannon") %>%
        ggplot(aes(
          x = samples_group,
          y = predicted,
          ymin = CI_low,
          ymax = CI_high,
          color = subset_value,
          alpha = p < 0.05
        )) +
        geom_errorbar(position = position_dodge(width = 0.6)) +
        geom_point(position = position_dodge(width = 0.6)) +
        facet_grid(samples_grouping ~ "", scales = "free", space = "free_y") +
        scale_color_kingdom() +
        scale_alpha_discrete(range = c(0.3, 1)) +
        coord_flip() +
        theme(panel.grid.major.x = element_line()) +
        labs(
          x = "Samples group",
          y = "Alpha diversity (z-score Shannon)",
          alpha = "is significant",
          color = "Kingdom"
        )
    },

    alphadiv_lmm_level_plt = {
      alphadiv_pre_tbl <-
        alphadiv %>%
        filter(subset_name == "kingdom") %>%
        group_by(subset_value) %>%
        slice(1) %>%
        ungroup()

      alphadiv_model_tbl <-
        alphadiv_pre_tbl %>%
        expand_grid(
          samples_grouping = c("habitat"),
          alphadiv_metric = c("Shannon")
        ) %>%
        mutate(
          data = data %>% map(~ .x %>% inner_join(samples)),

          model = list(data, alphadiv_metric, samples_grouping) %>% pmap(~ {
            str_glue("{..2} ~ 0 + {..3} + (1|bioproject_id)") %>%
              as.formula() %>%
              lme4::lmer(data = ..1) %>%
              list()
          }),

          parameters = model %>% map2(samples_grouping, ~ {
            .x %>%
              first() %>%
              parameters::parameters(p_adjust = "fdr") %>%
              as_tibble()
          }),
        ) %>%
        transmute(
          alphadiv_metric,
          subset_value,
          samples_grouping,
          model,
          parameters
        )


      alphadiv_model_tbl %>%
        mutate(posthoc = model %>% map2(samples_grouping, ~ {
          str_glue("pairwise ~ 0 + {.y}") %>%
            as.formula() %>%
            emmeans::emmeans(object = .x[[1]], adjust = "tukey") %>%
            pluck("contrasts") %>%
            as_tibble() %>%
            separate(contrast, into = c("from", "to"), sep = " - ")
        })) %>%
        unnest(posthoc) %>%
        ggplot(aes(from, to, fill = estimate, alpha = p.value < 0.05)) +
        geom_tile() +
        facet_grid(~subset_value) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        labs(fill = "Shannon")
    },

    #
    # Analyze beta diversity ----
    #
    sub_phys = target(
      {
        sub_abundances %>%
          mutate(
            phy = data %>% map(possibly(~ {
              tbl_to_phy(
                abundance_tbl = .x %>% select(sample_id, taxon, abundance),
                samples_tbl = samples,
                lineages_tbl = lineages,
                pooling_col = pooling_col
              )
            }, NA))
          ) %>%
          select(-data)
      },
      dynamic = map(sub_abundances),
      hpc = TRUE
    ),

    distances = target(
      {
        subsets %>%
          inner_join(sub_phys) %>%
          filter(norm_method == dist_norm_methods) %>%
          expand_grid(
            dist_method = dist_methods,
            dist_type = dist_types
          ) %>%
          mutate(
            dist = list(phy, dist_method, dist_type) %>% pmap(possibly(
              ~ phyloseq::distance(..1, method = ..2, type = ..3), NA
            ))
          ) %>%
          arrange_columns()
      },
      dynamic = cross(subsets, dist_types, dist_methods, dist_norm_methods),
      hpc = TRUE
    ),

    distances_adonis = target(
      {
        if (
          sub_abundances$samples_group[[1]] != "all" |
            sub_abundances$norm_method[[1]] != "raw"
        ) {
          return(tibble())
        }

        res <- possibly(test_adonis, NA)(
          formula = abundance ~ samples_group,
          data = sub_abundances$data[[1]],
          method = adonis_dist_method,
          parallel = getOption("mc.cores")
        )

        sub_abundances %>%
          discard(is.list) %>%
          mutate(
            adonis = list(res),
            dist_method = adonis_dist_method
          ) %>%
          arrange_columns()
      },
      dynamic = map(sub_abundances),
      hpc = FALSE
    ),

    distances_by_taxa_group = target(
      {
        get_distances_plot_by_taxa_group(
          distances = distances,
          taxa_grouping = taxa_grouping_pairs$taxa_grouping[[1]],
          taxa_group_x = taxa_grouping_pairs$taxa_group_x[[1]],
          taxa_group_y = taxa_grouping_pairs$taxa_group_y[[1]]
        )
      },
      dynamic = map(taxa_grouping_pairs),
      format = "rds", # qs leads to RAM overflow
      hpc = FALSE
    ),

    dominances_summary = {
      sub_abundances %>%
        filter(
          taxa_grouping == "dominance" &
            norm_method == "tss"
        ) %>%
        group_by(samples_grouping, samples_group) %>%
        nest() %>%
        mutate(
          data = data %>% map(~ {
            .x %>%
              unnest(data) %>%
              group_by(dominance, kingdom) %>%
              dplyr::summarize(
                abundance = sum(abundance)
              )
          })
        ) %>%
        unnest(data) %>%
        group_by(samples_grouping, samples_group, kingdom) %>%
        nest() %>%
        mutate(pie_plt = data %>% map2(samples_group, ~ {
          .x %>%
            ggpie("abundance", fill = "dominance", label = "dominance") +
            scale_fill_manual(values = dominance_colors) +
            theme(axis.text = element_blank()) +
            labs(title = .y)
        }))
    },

    differential_abundances_sub_phys = {
      sub_phys %>%
        expand_grid(samples_covariate = c("environment_group", "habitat")) %>%
        filter(norm_method == "raw") %>%
        select(-norm_method) %>%
        # filter combos with insufficient subsets
        filter(!(subset_name %>% str_detect("habitat") & samples_covariate == "habitat")) %>%
        filter(!(subset_name %>% str_detect("(environment_group|habitat)") & samples_covariate == "environment_group")) %>%
        # arrange columns
        select(where(~ !.x %>% is.list()), where(is.list))
    },

    ancom_differential_abundances = target(
      {
        differential_abundances_sub_phys %>%
          mutate(
            ancombc = phy %>% map2(samples_covariate, possibly(~ ANCOMBC::ancombc(phyloseq = .x, formula = str_glue("{.y} - 1"), NA)))
          ) %>%
          select(-where(is.list), ancombc)
      },
      dynamic = map(differential_abundances_sub_phys),
      hpc = TRUE
    ),

    differential_abundances = {
      ancom_differential_abundances %>%
        mutate(
          diff = ancombc %>% map(~ {
            list(
              .x %>% parse_ancombc_result("p_val"),
              .x %>% parse_ancombc_result("q_val"),
              .x %>% parse_ancombc_result("W"),
              .x %>% parse_ancombc_result("diff_abn"),
              .x %>% parse_ancombc_result("beta")
            ) %>%
              reduce(full_join) %>%
              mutate(
                samples_covariate_level = samples_covariate_level %>% map_chr(~ {
                  .x %>%
                    str_replace("environment_group", "environment_group,") %>%
                    str_replace("habitat", "habitat,")
                })
              ) %>%
              separate(samples_covariate_level, into = c("samples_grouping", "samples_group"), sep = ",")
          })
        ) %>%
        select(-ancombc) %>%
        unnest(diff) %>%
        group_by(taxon) %>%
        # do not count diff taxa twice (e.g. in sample groupings kingdom and kingdom,environment_group)
        filter(subset_name == "kingdom") %>%
        mutate(
          result = case_when(
            beta > 0 & diff_abn ~ "higher",
            beta < 0 & diff_abn ~ "lower",
            TRUE ~ "n.s."
          ) %>% factor(levels = c("lower", "n.s.", "higher")),
          samples_group = samples_group %>% factor(levels = c(environment_groups, habitats$habitat))
        ) %>%
        # is identical wuth samples_grouping
        select(-samples_covariate) %>%
        group_by(subset_name, subset_value, samples_grouping) %>%
        nest() %>%
        rename(differential_abundance = data)
    },

    differential_abundances_box_plot = {
      tibble() %>%
        ggplot(aes(samples_group, beta, color = subset_value)) +
        geom_hline(yintercept = 0) +
        # draw up and down regulated taxa separateley (bimodal distribution)
        geom_boxplot(data = differential_abundances %>% filter(diff_abn & beta > 0)) +
        geom_boxplot(data = differential_abundances %>% filter(diff_abn & beta < 0)) +
        coord_flip() +
        scale_color_kingdom() +
        facet_grid(
          # name empty. Will be NA otherwise
          samples_grouping ~ "empty",
          space = "free", scales = "free",
          labeller = as_labeller(c("environment_group" = "Environ.", "habitat" = "Habitat", "empty" = ""))
        ) +
        scale_y_continuous(expand = c(-0, 0), limits = c(-8, 8)) +
        theme(
          strip.text.x = element_blank()
        ) +
        labs(
          title = "Differentially abundant genera",
          x = "Sample group",
          color = "Kingdom",
          y = TeX("Effect size $\\beta$")
        )
    },

    differential_abundances_bar_plot = {
      differential_abundances %>%
        ggplot(aes(samples_group, fill = result)) +
        geom_bar() +
        scale_fill_manual(values = c("lower" = "#457b9d", "n.s." = "grey", "higher" = "#e63946")) +
        coord_flip() +
        facet_grid(
          # name empty. Will be NA otherwise
          samples_grouping ~ subset_value,
          space = "free", scales = "free",
          labeller = as_labeller(c("environment_group" = "Environ.", "habitat" = "Habitat", "Bacteria" = "Bacteria", "Fungi" = "Fungi"))
        ) +
        scale_y_continuous(expand = c(-0, 0)) +
        scale_x_discrete(expand = c(-0, 0)) +
        labs(
          title = "Differentially abundant genera",
          x = "Sample group",
          fill = "Result",
          y = "Taxa count"
        )
    },

    environment_group_specific_taxa_summary = {
      taxa_deck <-
        prevalences %>%
        unnest(data) %>%
        filter(subset_name == "environment_group,kingdom") %>%
        separate(subset_value, into = c("environment_group", "kingdom"), sep = ",") %>%
        distinct(environment_group, taxon, kingdom)

      expand_grid(
        kingdom = kingdoms_groups,
        environment_group = environment_groups
      ) %>%
        mutate(
          specific_taxa = list(environment_group, kingdom) %>% pmap(~ {
            cur_subset_taxa <-
              taxa_deck %>%
              filter(environment_group == .x & kingdom == .y) %>%
              pull(taxon) %>%
              unique()
            other_subsets_taxa <-
              taxa_deck %>%
              filter(environment_group != .x & kingdom == .y) %>%
              pull(taxon) %>%
              unique()
            cur_specific_subset_taxa <-
              cur_subset_taxa %>%
              setdiff(other_subsets_taxa) %>%
              unique()

            list(
              n_specific = length(cur_specific_subset_taxa),
              n_all = length(cur_subset_taxa)
            )
          })
        ) %>%
        unnest_wider(specific_taxa) %>%
        mutate(fraction = n_specific / n_all)
    },
    
    environment_group_shared_genera_vs_sample_inbalance = {
      # Are there really more shared bacteria between soil and aquatic genera than
      # there are genera shared between host and aquatic? Is this influenced by sample size inbalance?
      
      # down sampling
      n_samples <-
        sub_abundances$data[[1]] %>%
        filter(abundance > 0 & kingdom %in% kingdoms_groups & environment_group %in% environment_groups) %>%
        distinct(sample_id, environment_group) %>%
        count(environment_group) %>%
        pull(n) %>%
        min() %>%
        `*`(0.5) %>% # allow randomnes even in the smallest subset
        as.integer()
      
      tibble(trail = seq(1000)) %>%
        mutate(
          data = trail %>% map(~ {
            set.seed(.x)
            
            sub_abundances$data[[1]] %>%
              filter(
                abundance > 0 &
                  kingdom %in% kingdoms_groups &
                  environment_group %in% environment_groups
              ) %>%
              # sub sampling  
              nest(-c(environment_group, sample_id)) %>%
              group_by(environment_group) %>%
              sample_n(n_samples) %>%
              unnest(data) %>%
              distinct(environment_group, kingdom, taxon) %>%
              mutate(is_prevalent = TRUE) %>%
              pivot_wider(
                names_from = environment_group,
                values_from = is_prevalent,
                values_fill = list(is_prevalent = FALSE)
              ) %>%
              distinct(taxon, .keep_all = TRUE) %>%
              group_by(kingdom) %>%
              summarise(
                soil_aquatic = sum(soil & aquatic),
                soil_host = sum(soil & host),
                host_aquatic = sum(host & aquatic)
              ) %>%
              pivot_longer(c(soil_aquatic, soil_host, host_aquatic))
          })
        ) %>%
        unnest() %>%
        group_by(trail, kingdom) %>%
        arrange(-value) %>%
        summarise(sample_subset_order = paste0(name, collapse = ",")) %>%
        group_by(kingdom, sample_subset_order) %>%
        count() %>%
        group_by(kingdom) %>%
        mutate(freq = n / sum(n))
    }
  )
}
