#' Sub plan for FEAST source tracking
get_feast_plan <- function() {
  drake_plan(
    feast_source_group = "environment_group",

    feast_samples = {
      set.seed(1337)
      
      # Random sampling stratified by habitat
      feast_source_samples <-
        samples %>%
        group_by(habitat) %>%
        sample_frac(0.5) %>%
        pull(sample_id) %>%
        unique()
      
      samples %>%
        mutate(SourceSink = ifelse(sample_id %in% feast_source_samples, "Source", "Sink")) %>%
        group_by(SourceSink) %>%
        mutate(
          id = (SourceSink == "Sink") %>% ifelse(row_number(), NA)
        ) %>%
        select(
          SourceSink,
          sample_id,
          Env = feast_source_group,
          id
        )
    },

    feast_sinks = {
      feast_samples %>%
        filter(SourceSink == "Sink") %>%
        pull(sample_id) %>%
        unique()
    },

    feast_abundances = target(
      {
        res <-
          sub_abundances %>%
          filter(norm_method == "raw" & subset_value == kingdoms_groups) %>%
          head(1) %>%
          pull(data) %>%
          first() %>%
          select(sample_id, taxon, abundance)

        tibble(kingdom = kingdoms_groups, data = list(res))
      },
      dynamic = map(kingdoms_groups)
    ),

    feast_outputs = target(
      {
        metadata <-
          feast_samples %>%
          filter(sample_id == feast_sinks | SourceSink == "Source") %>%
          # each sink against all sources
          group_by(SourceSink) %>%
          mutate(id = row_number() %>% map2_int(SourceSink, ~ ifelse(.y == "Sink", .x, NA))) %>%
          as.data.frame() %>%
          set_rownames(.$sample_id) %>%
          select(-sample_id)

        otus <-
          feast_abundances %>%
          filter(kingdom == kingdoms_groups) %>%
          pull(data) %>%
          first() %>%
          filter(sample_id %in% rownames(metadata)) %>%
          pivot_wider(
            names_from = taxon,
            values_from = abundance,
            values_fill = list(abundance = 0)
          ) %>%
          as.data.frame() %>%
          set_rownames(.$sample_id) %>%
          select(-sample_id) %>%
          as.matrix()

        # ensure thread safety
        set.seed(1337)

        feast_dir <- "/analysis/feast"
        dir.create(feast_dir, showWarnings = FALSE, recursive = TRUE)

        res <- possibly(FEAST, NA)(
          C = otus,
          metadata = metadata,
          different_sources_flag = 0,
          dir_path = feast_dir,
          outfile = tempfile(tmpdir = feast_dir),
        )

        # FEAST changes wd
        setwd("/analysis/")

        if (is.na(res)) {
          return(tibble())
        }

        tibble(
          kingdom = kingdoms_groups,
          sink_id = feast_sinks,
          output = list(res)
        )
      },
      dynamic = cross(feast_sinks, kingdoms_groups)
    ),
    
    feast_outputs_file = {
      permuted_coabundances %>%
        write_rds(file_out("results/feast_outputs.rds"), compress = "gz")
    },

    feast_outputs_wide = {
      feast_outputs %>%
        mutate(output = output %>% map(~ {
          .x %>%
            as_tibble(rownames = "sink") %>%
            pivot_longer(-sink, names_to = "source", values_to = "contribution")
        })) %>%
        unnest(output) %>%
        left_join(feast_samples, by = c("sink_id" = "sample_id")) %>%
        select(-sink_id) %>%
        mutate(
          sink = sink %>% str_remove("_[A-z]+$"),
          source = source %>% str_remove("_[A-z]+$")
        ) %>%
        left_join(
          samples %>% set_colnames(map_chr(samples %>% colnames(), ~ .x %>% paste0("sink_", .))),
          by = c("sink" = "sink_sample_id")
        ) %>%
        left_join(
          samples %>% set_colnames(map_chr(samples %>% colnames(), ~ .x %>% paste0("source_", .))),
          by = c("source" = "source_sample_id")
        ) %>%
        left_join(
          samples %>% select(sample_id, source_habitat = habitat),
          by =  c("source" = "source_habitat")
        ) %>%

        # Normalizing and pooling
        group_by(kingdom, source_environment_group, source_habitat, sink_environment_group, sink_habitat) %>%
        summarise(contribution = sum(contribution)) %>%
        group_by(kingdom, sink_environment_group, sink_habitat) %>%
        mutate(contribution = contribution / sum(contribution) * 100) %>%
        ungroup() %>%
        mutate(contribution = contribution / sum(contribution)) %>%
        mutate_at(c("sink_habitat", "source_habitat"), ~ .x %>% factor(levels = habitats$habitat %>% unique())) %>%
        replace_na(
          list(
            sink_habitat = unknown_habitat,
            source_habitat = unknown_habitat,
            sink_environment_group = unknown_habitat,
            source_environment_group = unknown_habitat
          )
        )
    },

    feast_plt_sankey = {
      feast_outputs_wide %>%
        # arrange strata
        mutate(source_environment_group = source_environment_group %>%
          factor(levels = c("unknown", "aquatic", "host", "soil"))) %>%
        # 100% is bacteria and fungi combined
        mutate(contribution = contribution * 2) %>%
        ggplot(aes(
          y = contribution,
          axis1 = source_environment_group,
          axis2 = source_habitat,
          axis3 = sink_habitat,
          axis4 = sink_environment_group
        )) +
        geom_alluvium(aes(fill = sink_environment_group)) +
        geom_stratum() +
        facet_wrap(~kingdom, nrow = 1) +
        ggfittext::geom_fit_text(stat = "stratum", aes(label = after_stat(stratum)), reflow = TRUE) +
        scale_x_discrete(limits = c("Source environment", "Source habitat", "Sink habitat", "Sink environment"), expand = c(0, 0)) +
        scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
        scale_fill_environment_group() +
        labs(
          y = "Source contribution",
          fill = "Environment"
        )
    },

    feast_plt_bar = {
      feast_outputs %>%
        mutate(output = output %>% map(~ {
          .x %>%
            as_tibble(rownames = "sink") %>%
            pivot_longer(-sink, names_to = "source", values_to = "contribution")
        })) %>%
        unnest(output) %>%
        left_join(feast_samples, by = c("sink_id" = "sample_id")) %>%
        select(-sink_id) %>%
        mutate(
          sink = sink %>% str_remove("_[A-z]+$"),
          source = source %>% str_remove("_[A-z]+$")
        ) %>%
        left_join(
          samples %>% set_colnames(map_chr(samples %>% colnames(), ~ .x %>% paste0("sink_", .))),
          by = c("sink" = "sink_sample_id")
        ) %>%
        left_join(
          samples %>% set_colnames(map_chr(samples %>% colnames(), ~ .x %>% paste0("source_", .))),
          by = c("source" = "source_sample_id")
        ) %>%

        # Normalizing and pooling
        group_by(kingdom, source_environment_group, sink_environment_group, sink_habitat) %>%
        summarise(contribution = sum(contribution)) %>%
        group_by(kingdom, sink_environment_group, sink_habitat) %>%
        mutate(contribution = contribution / sum(contribution)) %>%
        mutate(
          contribution = contribution * 100
        ) %>%
        ggplot(aes(
          x = sink_habitat,
          y = contribution,
          fill = source_environment_group
        )) +
        geom_bar(stat = "identity") +
        scale_fill_environment_group() +
        facet_grid(sink_environment_group ~ kingdom, scales = "free_y", space = "free") +
        coord_flip() +
        scale_y_continuous(expand = c(0, 0)) +
        theme(panel.spacing.x = unit(0.8, "cm")) +
        labs(
          x = "Sink habitat",
          y = "Source contribution (FEAST, %)",
          fill = "Source environment"
        )
    }
  )
}
