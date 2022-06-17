get_enrichment_plan <- function() {
  drake::drake_plan(
    enrichment_taxa = {
      enrichment_taxa_faprotax <-
        read_csv("raw/databases/faprotax.csv") %>%
        transmute(
          db = "FAPROTAX",
          kingdom = "Bacteria",
          taxon = phylogeny %>% str_extract("g__[^;]+;") %>% str_remove("^g__") %>% str_remove(";$"),
          entry
        )
      
      enrichment_taxa_fusiondb_oxygen <-
        read_csv("raw/databases/FusionDB.csv") %>%
          transmute(
            db = "FusionDB-oxygen",
            kingdom = "Bacteria",
            taxon = `Organism Name` %>% str_extract("^[A-z]+"),
            entry = OxygenPreference %>% str_extract("^[A-z]+")
          ) %>%
        filter(! is.na(entry))
      
      enrichment_taxa_fusiondb_temperature <-
        read_csv("raw/databases/FusionDB.csv") %>%
        transmute(
          db = "FusionDB-temperature",
          kingdom = "Bacteria",
          taxon = `Organism Name` %>% str_extract("^[A-z]+"),
          entry = TemperaturePreference %>% str_extract("^[A-z]+")
        ) %>%
        filter(! is.na(entry))

      enrichment_taxa_funguild_guild <-
        read_csv("raw/databases/FUNGuild.csv") %>%
        transmute(
          db = "FUNGuild-guild",
          kingdom = "Fungi",
          taxon = genus,
          entry = guild
        )
      
      enrichment_taxa_funguild_trophicMode <-
        read_csv("raw/databases/FUNGuild.csv") %>%
        transmute(
          db = "FUNGuild-trophicMode",
          kingdom = "Fungi",
          taxon = genus,
          entry = trophicMode
        )

      enrichment_taxa_phyla <-
        lineages %>%
        transmute(
          db = "phylum",
          kingdom,
          taxon,
          entry = phylum
        )

      enrichment_taxa <-
        list(
          enrichment_taxa_faprotax,
          enrichment_taxa_phyla,
          enrichment_taxa_funguild_guild,
          enrichment_taxa_funguild_trophicMode,
          enrichment_taxa_fusiondb_oxygen,
          enrichment_taxa_fusiondb_temperature
        ) %>%
        bind_rows() %>%
        group_by(db, kingdom) %>%
        nest() %>%
        rename(db_data = data)

      enrichment_taxa
    },

    enrichment_universe = {
      lineages %>%
        group_by(kingdom) %>%
        nest() %>%
        transmute(
          universe = data %>% map(~ .x$taxon)
        )
    },

    enrichment_taxa_sets = {
      enrichment_taxa_sets_lmm_abundance <-
        lmm_abundance %>%
        filter(term == "environment_group") %>%
        filter(q.value < 0.05) %>%
        group_by(kingdom) %>%
        nest() %>%
        transmute(
          set_group = "lmm diff abundant",
          set_name = "sign. impact environment_group",
          set_value = data %>% map(~ .x$taxon %>% unique())
        )
      
      enrichment_taxa_sets_modules <-
        wgcna_modules %>%
        unnest(modules) %>%
        mutate(module = unmerged_module) %>%
        group_by(kingdom, environment_group, module) %>%
        nest() %>%
        rename(set_name = module) %>%
        transmute(
          set_group = "consensus module",
          set_name = set_name %>% as.character(),
          set_value = data %>% map(~ .x$taxon %>% unique())
        )
      
      enrichment_taxa_sets_preserved_modules <-
        wgcna_preservations %>%
        unnest(preservation) %>%
        unnest(data) %>%
        mutate(is_preserved = as.character(p.value < 0.05)) %>%
        inner_join(enrichment_taxa_sets_modules %>% rename(module = set_name)) %>%
        group_by(kingdom, environment_group, is_preserved) %>%
        nest() %>%
        rename(set_name = is_preserved) %>%
        transmute(
          set_group = "preserved consensus module",
          set_name = set_name %>% as.character(),
          set_value = data %>% map(~ .x$set_value %>% reduce(union))
        )
      
      enrichment_taxa_sets_coabundaces <-
        graphs %>%
        filter(
          subset_name == "environment_group,kingdom" &
            cor_method == "sparcc" &
            norm_method == "raw"
        ) %>%
        separate(subset_value, into = c("environment_group", "kingdom"), sep = ",") %>%
        transmute(
          environment_group,
          kingdom,
          taxa = cor_graph %>% map(~ .x %>%
            activate(nodes) %>%
            pull(taxon))
        ) %>%
        unnest(taxa) %>%
        group_by(kingdom) %>%
        nest() %>%
        transmute(
          taxa = data %>% map(~ {
            .x %>%
              group_by(environment_group) %>%
              nest() %>%
              mutate(data = data %>% map(~ .x[["taxa"]])) %>%
              deframe()
          }),

          res = taxa %>% map(function(taxa) {
            intersections <-
              environment_groups %>%
              combn(2) %>%
              t() %>%
              as_tibble() %>%
              transmute(
                name = str_glue("{V1} ∩ {V2}"),
                value = V1 %>% map2(V2, ~ taxa[[.x]] %>% intersect(taxa[[.y]]))
              )

            exclusives <-
              environment_groups %>%
              tibble(name = .) %>%
              mutate(
                value = name %>% map(~ {
                  current_env <- .x
                  other_envs <- environment_groups %>% discard(~ .x == current_env)
                  other_taxa <- taxa[other_envs] %>% reduce(union)

                  taxa[current_env] %>%
                    setdiff(other_taxa) %>%
                    simplify()
                }),
                name = str_glue("only {name}")
              )

            all <- tibble(
              name = "all",
              value = taxa[environment_groups] %>% reduce(intersect) %>% list()
            )

            list(
              all,
              intersections,
              exclusives
            ) %>%
              bind_rows()
          })
        ) %>%
        select(kingdom, res) %>%
        unnest(res) %>%
        transmute(
          set_group = "environment coabundance",
          set_name = as.character(name),
          set_value = value
        )

      enrichment_taxa_sets <- list(
        enrichment_taxa_sets_lmm_abundance,
        enrichment_taxa_sets_modules,
        enrichment_taxa_sets_preserved_modules,
        enrichment_taxa_sets_coabundaces
      ) %>%
        reduce(full_join)

      enrichment_taxa_sets
    },

    enrichments = {
      enrichment_taxa_sets %>%
        left_join(enrichment_taxa) %>%
        left_join(enrichment_universe) %>%
        transmute(
          db,
          kingdom,
          environment_group,
          set_group,
          set_name,
          enrichment = list(set_value, db_data, universe) %>% pmap(possibly(~ {
            term2gene <- ..2 %>% select(term = entry, gene = taxon)
            test_taxa <- ..1 %>% unique()
            universe <- ..3

            enrich_msea(term2gene, test_taxa, length(universe))

            # statistical dependency between the two variables
            # term2gene$gene %>%
            #   union(test_taxa) %>%
            #   tibble(taxon = .) %>%
            #   mutate(
            #     in_ref = taxon %in% term2gene$gene,
            #     in_test = taxon %in% test_taxa
            #   ) %>%
            #   xtabs(~in_ref+in_test, data = .) %>%
            #   fisher.test(alternative = "greater") %>%
            #   tidy()

            # Over Representation Analysis
            # clusterProfiler::enricher(
            #   gene = test_taxa,
            #   TERM2GENE = term2gene,
            #   universe = universe,
            #   pvalueCutoff = 0.05,
            #   qvalueCutoff = 0.2,
            #   minGSSize = 1,
            #   maxGSSize = length(test_taxa),
            # ) %>%
            #   clusterProfiler::summarise()

            # taxon set enrichment analysis
            # clusterProfiler::GSEA(
            #   geneList = test_taxa %>% length() %>% runif() %>% set_names(test_taxa) %>% sort(decreasing = TRUE),
            #   TERM2GENE = term2gene
            # ) %>%
            #   clusterProfiler::summarise()
          }, NA))
        )
    },
    
    enrichment_preserved_plt = {
      enrichments %>%
        unnest(enrichment) %>%
        filter(q.value < 0.05 & set_group == "preserved consensus module") %>%
        group_by(kingdom, environment_group, set_name, db, term) %>%
        summarise(n_taxa = sum(n_shared)) %>%
        filter(n_taxa >= 5) %>%
        mutate(set_name = set_name %>% as.logical()) %>%
        
        # add missing zeros
        group_by(kingdom, term) %>%
        complete(kingdom, environment_group, set_name, db, term, fill = list(n_taxa = 0)) %>%
        
        ggplot(aes(term, n_taxa, color = environment_group, shape = set_name)) +
        geom_point(size = 3) +
        coord_flip() +
        scale_color_environment_group() +
        scale_y_continuous(expand = c(0.01, 0)) +
        facet_grid(kingdom ~ "", scales = "free_y", space = "free") +
        theme(panel.grid.major = element_line()) +
        labs(
          title = "Enrichment of preserved taxa modules",
          shape = "In habitat\npreserved module",
          color = "Environment",
          x = "Term",
          y = "# Taxa"
        )
    }
  )
}
