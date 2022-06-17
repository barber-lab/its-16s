get_wgcna_plan <- function() {
  drake_plan(
    wgcna_quantile = 0.95,

    wgcna_graphs = {
      graphs %>%
        filter(subset_name == "habitat,kingdom" & cor_method == "sparcc") %>%
        separate(subset_value, into = c("habitat", "kingdom"), sep = ",") %>%
        mutate(environment_group = habitat %>% map_chr(habitat_to_environment_group)) %>%
        discard(~ .x %>%
          unique() %>%
          length() == 1) %>%
        select(-where(is.list), -where(is.numeric), cor_graph)
    },

    wgcna_habitats = {
      wgcna_graphs %>%
        select(-cor_graph) %>%
        group_by(kingdom, environment_group) %>%
        nest() %>%
        transmute(habitats = data %>% map(~ .x[[1]]))
    },

    wgcna_taxa = {
      wgcna_graphs %>%
        group_by(kingdom, environment_group) %>%
        nest() %>%
        transmute(
          taxa = data %>% map(~ {
            n_subsets <- length(.x$cor_graph)

            .x$cor_graph %>%
              map(~ .x %>%
                activate(nodes) %>%
                pull(taxon)) %>%
              # taxa present in the majority of habitats
              simplify() %>%
              table() %>%
              keep(~ .x > n_subsets / 2) %>%
              names()
          })
        )
    },

    wgcna_matricies = {
      wgcna_graphs %>%
        inner_join(wgcna_taxa) %>%
        mutate(
          matrix = cor_graph %>% map2(taxa, possibly(~ {
            current_edges <-
              .x %>%
              activate(nodes) %>%
              mutate(name = taxon) %>%
              filter(name %in% .y) %>%
              activate(edges) %>%
              filter(from_taxon %in% .y & to_taxon %in% .y) %>%
              # Topological Overlap Matrices needed values in range [0, 1]
              mutate(scaled_estimate = 0.5 * (estimate + 1)) %>%
              as_adjacency_matrix(attr = "scaled_estimate", names = TRUE) %>%

              # complete missing edges
              # Function complete is not available for obejcts of class tbl_graph
              as.matrix() %>%
              as_tibble(rownames = "from_taxon") %>%
              pivot_longer(-from_taxon, names_to = "to_taxon", values_to = "estimate")

            missing_edges <-
              expand_grid(
                from_taxon = .y %>% setdiff(current_edges$from_taxon),
                to_taxon = .y %>% setdiff(current_edges$to_taxon),
                estimate = 0
              )

            current_edges %>%
              bind_rows(missing_edges) %>%
              pivot_wider(names_from = to_taxon, values_from = estimate, values_fill = list(estimate = 0)) %>%
              column_to_rownames("from_taxon") %>%
              as.matrix()
          }, NA))
        ) %>%
        filter(!is.na(matrix)) %>%
        select(-where(is.list), -where(is.numeric), matrix)
    },

    wgcna_tom_similarities = {
      wgcna_matricies %>%
        mutate(
          similarity = matrix %>% map(possibly(~ {
            res <- WGCNA::TOMsimilarity(.x, TOMType = "unsigned")
            colnames(res) <- colnames(.x)
            rownames(res) <- rownames(.x)
            res
          }, NA))
        ) %>%
        select(-where(is.list), -where(is.numeric), similarity)
    },

    wgcna_scale_quants = {
      wgcna_tom_similarities %>%
        group_by(kingdom, environment_group) %>%
        arrange(kingdom, environment_group) %>%
        filter(!is.na(similarity)) %>%
        transmute(
          habitat,
          is_ref = row_number() == 1,
          scale_quant = similarity %>% map_dbl(~ .x %>%
            as.dist() %>%
            quantile(probs = wgcna_quantile, type = 8) %>%
            as.numeric())
        )
    },

    wgcna_ref_scale_quants = {
      wgcna_scale_quants %>%
        filter(is_ref) %>%
        transmute(ref_scale_quant = scale_quant)
    },

    wgcna_consensus = {
      wgcna_tom_similarities %>%
        inner_join(wgcna_ref_scale_quants) %>%
        inner_join(wgcna_scale_quants) %>%
        # scale to make them comparable for consensus network
        mutate(
          similarity = list(similarity, ref_scale_quant, scale_quant, is_ref) %>% pmap(~ {
            if (is_ref) {
              return(..1)
            } # do not scale reference

            scale_powers <- log(ref_scale_quant) / log(scale_quant)
            ..1^scale_powers
          })
        ) %>%
        group_by(kingdom, environment_group) %>%
        nest() %>%
        transmute(
          consensus = data %>% map(possibly(~ {
            consensus_tom <-
              .x$similarity %>%
              WGCNA::pquantile.fromList(prob = 0.5)

            consensus_tree <-
              consensus_tom %>%
              {
                (1 - .)
              } %>% # similarity to distance
              as.dist() %>%
              hclust(method = "average")

            list(
              consensus_tom = consensus_tom,
              consensus_tree = consensus_tree
            )
          }, NA))
        ) %>%
        arrange(kingdom, environment_group) %>%
        select(-where(is.list), -where(is.numeric), consensus) %>%
        unnest_wider(consensus)
    },

    wgcna_modules = {
      wgcna_consensus %>%
        transmute(
          modules = list(consensus_tree, consensus_tom, kingdom, environment_group) %>%
            pmap(function(consTree, consensusTOM, cur_kingdom, cur_environment_group) {
              multiExpr <-
                wgcna_matricies %>%
                filter(
                  kingdom == cur_kingdom &
                    environment_group == cur_environment_group
                ) %>%
                transmute(habitat, data = matrix %>% map(~ list(data = .x))) %>%
                deframe()

              unmergedLabels <- cutreeDynamic(
                dendro = consTree,
                distM = 1 - consensusTOM,
                deepSplit = 2,
                cutHeight = 0.995,
                minClusterSize = 10,
                pamRespectsDendro = FALSE
              )

              unmergedColors <- labels2colors(unmergedLabels) %>% muted()

              # merge similar modules
              merge <- mergeCloseModules(multiExpr, unmergedLabels, cutHeight = 0.25, verbose = 3)
              moduleLabels <- merge$colors

              tibble(
                taxon = consTree$labels,
                unmerged_module = unmergedLabels,
                merged_module = moduleLabels
              ) %>%
                left_join(
                  lineages %>% transmute(taxon, phylum)
                )
            })
        )
    },

    wgcna_consensus_plots = {
      wgcna_consensus %>%
        inner_join(wgcna_modules) %>%
        transmute(
          consensus_plt = consensus_tree %>% map2(modules, ~ {
            dendro_plt <-
              .x %>%
              ggdendro::ggdendrogram() +
              scale_x_discrete(expand = c(0, 0)) +
              scale_y_continuous(expand = c(0, 0)) +
              theme_void()

            annot_plt <-
              .y %>%
              mutate_all(as.character) %>%
              # define colors
              mutate_at(colnames(.) %>% setdiff(c("taxon", "phylum")), ~ .x %>%
                labels2colors() %>%
                muted(60, 90)) %>%
              mutate(phylum = phyla_colors[phylum]) %>%
              pivot_longer(-taxon) %>%
              mutate(
                taxon = taxon %>% factor(levels = .x$labels[.x$order]),
                name = name %>% factor(levels = c("phylum", "merged_module", "unmerged_module"))
              ) %>%
              replace_na(list(value = phyla_colors["other"])) %>%
              ggplot(aes(taxon, name, fill = value)) +
              geom_tile() +
              scale_x_discrete(expand = c(0, 0)) +
              scale_y_discrete(expand = c(0, 0)) +
              theme_pub() +
              theme(axis.text.x = element_blank()) +
              labs(x = "", y = "")

            wrap_plots(dendro_plt, annot_plt, ncol = 1)
          })
        )
    },

    wgcna_preservations_res = {
      wgcna_tom_similarities %>%
        group_by(kingdom, environment_group) %>%
        nest() %>%
        rename(similarity = data) %>%
        inner_join(wgcna_modules, copy = TRUE) %>%
        transmute(
          preservation = similarity %>% map2(modules, ~ {
            multiExpr <-
              .x %>%
              transmute(habitat, data = similarity %>% map(~ list(data = .x))) %>%
              deframe()

            multiColor <-
              multiExpr %>%
              names() %>%
              tibble(habitat = .) %>%
              mutate(modules = {
                .y %>%
                  pull(merged_module) %>%
                  list()
              }) %>%
              deframe()

            enableWGCNAThreads(nThreads = 10)

            # BUGFIX: create variable weightsRef in global env (defaults.R)
            # with value NULL to allow function modulePreservation to run
            modulePreservation(
              multiData = multiExpr,
              multiColor = multiColor,
              dataIsExpr = FALSE, # data is similarity and not expression
              networkType = "unsigned",
              referenceNetworks = 1,
              nPermutations = 200,
              randomSeed = 1,
              quickCor = 0,
              verbose = 3,
              parallelCalculation = TRUE
            )
          })
        )
    },

    wgcna_preservations = {
      wgcna_preservations_res %>%
        left_join(wgcna_habitats) %>%
        transmute(
          preservation = preservation %>% map2(habitats, possibly(function(preservation, habitats){
            habitats %>%
              length() %>%
              seq() %>%
              discard(~ .x == 1) %>%
              tibble(test_habitat_nr = .) %>%
              mutate(
                from_habitat = habitats[1],
                to_habitat = habitats[-1]
              ) %>%
              mutate(
                data = test_habitat_nr %>% map(~ {c
                  tibble(
                    module = preservation$preservation$Z[[1]][[.x]] %>% rownames(),
                    size = preservation$preservation$Z[[1]][[.x]]$moduleSize,
                    z = preservation$preservation$Z[[1]][[.x]]$Zsummary.pres,
                    # see https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1001057#s4
                    p.value = 10^(preservation$preservation$log.p[[1]][[.x]]$log.p.cor.kME)
                  )
                })
              ) %>%
              select(-test_habitat_nr)
          }, NA))
        ) %>%
        filter(!is.na(preservation))
    }
  )
}
