#!/usr/bin/env R

antismash_files <- function(path) {
  
  list.files(path = path, full.names = T, pattern = ".json$", recursive = T) %>% 
    as_tibble() %>% 
    mutate(genome_id = str_split(string = value, 
                                 pattern = "/", 
                                 simplify = T)[,3]) 
  
}

is_there_results <- function(file) {
  
  flat <- fromJSON(file, flatten = T) %>% 
    .$records %>% 
    flatten(recursive = T) %>%
    # Will remove empty columns
    # If there are hits, areas must not be empty
    select(where(~!all(lengths(.) == 0)))
  
  # If areas was empty, it was removed
  if("areas"  %in% colnames(flat)) {
    
    flat %>% 
      select(areas)
    
  } else {
    NA
  }
}

parse_antismash <- function(files) {
  
  files %>% 
    mutate(bgcs = map(
      .x = value, 
      .f = ~ is_there_results(.x))
    )
}

group_bgcs <- function(df) {
  
  df %>% 
    unnest(bgcs) %>% 
    unnest(areas) %>% 
    select(genome_id, products) %>% 
    unnest(cols = c(products)) %>% 
    group_by(genome_id, products) %>% 
    count() %>% 
    arrange(desc(n))
  
}

generate_annotation <- function(bgcs_by_genome, genomes2find_BGCs, generalist_genus, specialist_genus) {
  
  type2complete <- tibble(genus = c("Erythrobacter","Lachnospiraceae"),
                          type = c("Specialist","Generalist"))
  
  g <- generalist_genus %>% 
    mutate(type = "Generalist") %>% 
    rename(genus = generalist)
  
  s <- specialist_genus %>% 
    mutate(type = "Specialist") %>% 
    rename(genus = specialists)
  
  
  metadata <- genomes2find_BGCs %>% 
    left_join(bind_rows(g,s, type2complete))
  
  metadata %>% 
    left_join(bgcs_by_genome, by = c("accession" = "genome_id")) %>% 
    mutate(n = replace_na(n, 0))
  
  
}

scale_fill_prevalence_type <- function(drop = FALSE, ...) {
  scale_fill_manual(
    values =  purrr::simplify(prevalence_type),
    na.value = "#000000",
    drop = drop,
    ...
  )
}

scale_fill_prevalence_group <- function(drop = FALSE, ...) {
  scale_fill_manual(
    values =  purrr::simplify(prevalence_group_colors),
    na.value = "#000000",
    drop = drop,
    ...
  )
}

scale_color_prevalence_type <- function(drop = FALSE, ...) {
  scale_color_manual(
    values =  purrr::simplify(prevalence_type),
    na.value = "#000000",
    drop = drop,
    ...
  )
}

symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                    symbols = c("****", "***", "**", "*", "ns"))

plot_bgcs_gen_size_cds <- function(generalist_specialists_bgcs, genomes2find_BGCs, cds_by_genome) {
  
  gen <- genomes2find_BGCs %>% 
    distinct() %>% 
    select(accession, genus, type, kingdom) %>% 
    left_join(genome_length, by = c("accession" = "genome_id")) %>% 
    select(-c(value, cmd)) %>% 
    unnest(g_size) %>% 
    mutate(g_size = as.numeric(g_size)) %>% 
    mutate(g_size = g_size/10^6)
  
  a <- generalist_specialists_bgcs %>% 
    rename(Groups = type) %>% 
    group_by(accession, Groups, kingdom) %>% 
    summarise(value = sum(n)) %>%
    ggplot(aes(Groups, value, color = Groups)) +
    geom_boxplot(width=0.5) +
    stat_compare_means(
      symnum.args = symnum.args,
      method = "wilcox", 
      comparisons = list(c("Specialist", "Generalist"))) +
    scale_y_log10() +
    facet_wrap(c("kingdom"), scales="free_y") +
    scale_color_prevalence_type() +
    annotation_logticks(sides = "l") +
    ylab("Number of BCGs") +
    guides(color = FALSE) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          strip.text = element_text(size = rel(1.1)),
          axis.text = element_text(size = rel(1.1)),
          legend.position = "none",
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
  
  b <- gen %>% 
    rename(Groups = type) %>% 
    ggplot(aes(Groups, g_size)) +
    geom_boxplot(aes(color = Groups), width=0.5) +
    facet_wrap(c("kingdom"), scales="free_y") +
    stat_compare_means(
      symnum.args = symnum.args,
      method = "wilcox", 
      comparisons = list(c("Generalist", "Specialist"))) +
    ylab("Genome size (MB)") +
    guides(color = FALSE) +
    scale_color_prevalence_type() +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          strip.text = element_text(size = rel(1.1)),
          axis.text = element_text(size = rel(1.1)),
          legend.position = "none",
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
  
  cds <- genomes2find_BGCs %>% 
    distinct() %>% 
    select(accession, genus, type, kingdom) %>% 
    left_join(cds_by_genome, by = c("accession" = "genome_id")) %>% 
    select(-value)
  
  c <- cds %>% 
    rename(Group = type) %>% 
    ggplot(aes(Group, cds_n)) +
    geom_boxplot(aes(color = Group), width=0.5) +
    facet_wrap(c("kingdom"), scales="free_y") +
    stat_compare_means(
      symnum.args = symnum.args,
      method = "wilcox",
      comparisons = list(c("Generalist", "Specialist"))) +
    ylab("CDS number") +
    guides(color = "none") +
    scale_color_prevalence_type() +
    scale_y_log10() +
    theme_bw() +
    scale_y_continuous(labels = scales::label_comma()) +
    theme(axis.title.x = element_blank(),
          strip.text = element_text(size = rel(1.1)),
          axis.text = element_text(size = rel(1.1)),
          legend.position = "none",
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
    
  
  p <- wrap_plots(a + labs(tag = "A"),
                  b + labs(tag = "B"),
                  c + labs(tag = "C"))
  
  
  ggsave(plot = p,
         filename = "plots/02-generalist-specialist-bgcs/figure1.png",
         device = "png", 
         width = 12, 
         height = 5, 
         dpi = 300)
  
  ggsave(plot = p,
         filename = "plots/02-generalist-specialist-bgcs/figure1.pdf",
         device = "pdf", 
         width = 12, 
         height = 5, 
         dpi = 300)
  
  ggsave(plot = b,
         filename = "plots/02-generalist-specialist-bgcs/gene_spec_gen_size.pdf",
         device = "pdf", 
         width = 12, 
         height = 5, 
         dpi = 300)
  
  return(list(a,b,c))
  
}

ratio_bgc_cds <- function(generalist_specialists_bgcs, cds_by_genome) {
  
  a <- generalist_specialists_bgcs %>% 
    select(accession, genus, type, kingdom, n) %>% 
    rename(bgcs_n = n) %>% 
    left_join(
      cds_by_genome %>% 
        rename(accession = genome_id) %>% 
        select(-value)
    ) %>% 
    mutate(ratio_bgc_cds = bgcs_n/cds_n*100) %>% 
    rename(Groups = type) %>%
    ggplot(aes(Groups, ratio_bgc_cds)) +
    geom_boxplot(aes(color = Groups), width=0.5) +
    facet_wrap(c("kingdom")) +
    stat_compare_means(
      symnum.args = symnum.args,
      method = "wilcox",
      comparisons = list(c("Generalist", "Specialist"))) +
    ylab("Ratio BCGs/CDS (%)") +
    guides(color = "none") +
    scale_color_prevalence_type() +
    scale_y_log10() +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 10),
          legend.position = "none",
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
  
  ggsave(plot = a,
         filename = "plots/02-generalist-specialist-bgcs/figure2.png",
         device = "png", 
         width = 7, 
         height = 5, 
         dpi = 300)
  
}

plot_cds <- function(genomes2find_BGCs, cds_by_genome, genome_length) {
  
  cds <- genomes2find_BGCs %>% 
    distinct() %>% 
    select(accession, genus, type, kingdom) %>% 
    left_join(cds_by_genome, by = c("accession" = "genome_id")) %>% 
    select(-value)
  
  a <- cds %>% 
    rename(Group = type) %>% 
    ggplot(aes(Group, cds_n)) +
    geom_boxplot(aes(color = Group), width=0.5) +
    facet_wrap(c("kingdom")) +
    stat_compare_means(
      symnum.args = symnum.args,
      method = "wilcox",
      comparisons = list(c("Generalist", "Specialist"))) +
    ylab("CDS") +
    guides(color = "none") +
    scale_color_prevalence_type() +
    scale_y_log10() +
    theme_bw() +
    scale_y_continuous(labels = scales::label_comma()) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 10),
          legend.position = "none",
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
  
  gen <- genomes2find_BGCs %>% 
    distinct() %>% 
    select(accession, genus, type, kingdom) %>% 
    left_join(genome_length, by = c("accession" = "genome_id")) %>% 
    select(-c(value, cmd)) %>% 
    unnest(g_size) %>% 
    mutate(g_size = as.numeric(g_size)) %>% 
    mutate(g_size = g_size/10^6)
  
  b <- gen %>% 
    rename(Groups = type) %>% 
    ggplot(aes(Groups, g_size)) +
    geom_boxplot(aes(color = Groups), width=0.5) +
    facet_wrap(c("kingdom"), scales="free_y") +
    stat_compare_means(
      symnum.args = symnum.args,
      method = "wilcox", 
      comparisons = list(c("Generalist", "Specialist"))) +
    ylab("Genome size (MB)") +
    guides(color = FALSE) +
    scale_color_prevalence_type() +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 10),
          legend.position = "none",
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
  
  
  p <- wrap_plots(a + labs(tag = "A"),
                  b + labs(tag = "B"))
  
  ggsave(filename = "plots/02-generalist-specialist-bgcs/figure3.png", 
         plot = p, 
         device = "png", 
         width = 9, 
         height = 5,
         dpi = 300)
}


plot_bcgs_gen_spec <- function(generalist_specialists_bgcs, bcgs_groups) {

  n <- 6
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  df <- generalist_specialists_bgcs %>% 
    rename(Groups = type) %>%
    group_by(Groups, kingdom, products) %>% 
    summarise(value = sum(n)) %>%
    filter(value > 0) %>% 
    mutate(perc = value / sum(value) * 100) %>% 
    left_join(bcgs_groups) %>%
    group_by(Groups, kingdom, cluster, perc) %>% 
    summarise(sum_perc = sum(perc))
  
  p <- df %>% 
    ggplot(aes(Groups, perc, fill = cluster)) +
    geom_col(aes(fill = fct_reorder(cluster, perc))) +
    facet_wrap(c("kingdom")) +
    scale_fill_manual(values = col_vector) +
    guides(fill=guide_legend(ncol=1)) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          strip.text = element_text(size = rel(1.1)),
          axis.text = element_text(size = rel(1.1)), 
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) +
    labs(y = "Biosynthetic cluster genes (%)",
         fill = "Products")
  
  ggsave(filename = "plots/02-generalist-specialist-bgcs/figure4.png", 
         plot = p, 
         device = "png", 
         width = 7, 
         height = 5,
         dpi = 300)
  
  ggsave(filename = "plots/02-generalist-specialist-bgcs/figure4.pdf", 
         plot = p, 
         device = "pdf", 
         width = 7, 
         height = 5,
         dpi = 300)
  
  return(p)

}

bgcs2tsv <- function(generalist_specialists_bgcs) {
  
  generalist_specialists_bgcs %>% 
    distinct(products) %>% 
    drop_na(products) %>% 
    arrange(products) %>% 
    write_tsv("tbl/01-bcgs-list/antismash-bcgs-products.tsv")
  
}

bcgs_statistics <- function(generalist_specialists_bgcs, bcgs_groups) {
  
  df <- generalist_specialists_bgcs %>% 
    rename(Groups = type) %>%
    drop_na(products) %>% 
    left_join(bcgs_groups) %>% 
    select(accession, Groups, kingdom, cluster, n) %>% 
    group_by(accession, Groups, kingdom, cluster) %>% 
    summarise(count = sum(n)) %>% 
    mutate(comp = str_glue("{cluster}-{Groups}"))
  
  
  p <- df %>% 
    ggplot(aes(comp, count, fill = cluster)) +
    geom_boxplot() +
    stat_compare_means(
      method = "wilcox",
      comparisons = list(
        c("RiPP-Specialist", "RiPP-Generalist"),
        c("NRPS-Specialist", "NRPS-Generalist"),
        c("Other-Specialist", "Other-Generalist"),
        c("PKS-Specialist", "PKS-Generalist"),
        c("Siderophore-Specialist", "Siderophore-Generalist"),
        c("Terpene-Specialist", "Terpene-Generalist")
      )) +
    facet_wrap(c("kingdom")) +
    coord_flip() +
    scale_y_log10()
  
  ggsave(filename = "plots/02-generalist-specialist-bgcs/bcgs_statistics.png", 
         plot = p, 
         device = "png", 
         width = 11, 
         height = 7)
}

import_genome_length <- function(path) {
  
  list.files(path, full.names = T) %>% 
    as_tibble() %>% 
    mutate(gen_len = map(value, ~ read_lines(.x))) %>% 
    unnest(gen_len) %>% 
    mutate(accession = str_extract(value, "GC[A|F]_[0-9]+.[0-9]"))
}

plot_main_fig2_paper <- function(generalist_specialists_bgcs, genomes2find_BGCs, cds_by_genome, genome_length, valid_amr) {
  
  prevalence_group_colors <- c(
    "generalist" = "#3a86ff",
    "specialist" = "#ffbe0b"
  )
  
  
  scale_color_prevalence_group <- function(drop = FALSE, ...) {
    scale_color_manual(
      values =  purrr::simplify(prevalence_group_colors),
      na.value = "#000000",
      drop = drop,
      ...
    )
  }
  
  symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                      symbols = c("****", "***", "**", "*", "ns"))
  
  a <- 
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
          symnum.args = symnum.args,
          comparisons = list(c("generalist", "specialist"))
        ) +
        annotation_logticks(sides = "l") +
        scale_color_prevalence_group() +
        theme_bw() +
        theme(
          strip.text = element_text(size = rel(1.1)),
          axis.text = element_text(size = rel(1.1)), 
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_blank(),
          axis.ticks.x=element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) +
        guides(color = FALSE) +
        labs(
          x = "",
          y = "Abundance (TSS)"
        )
    })
  
  
  
  cds <- genomes2find_BGCs %>% 
    distinct() %>% 
    select(accession, genus, type, kingdom) %>% 
    left_join(cds_by_genome, by = c("accession" = "genome_id")) %>% 
    select(-value) %>% 
    left_join(genome_length %>% 
                rename(accession = genome_id) %>% 
                unnest(g_size) %>% 
                mutate(g_size = as.numeric(g_size)) %>% 
                mutate(g_size = g_size/10^6) %>% 
                distinct(accession, g_size)) %>% 
    mutate(cds_norm = cds_n/g_size)
  
  b <- cds %>% 
    rename(Group = type) %>% 
    ggplot(aes(Group, cds_n)) +
    geom_boxplot(aes(color = Group), width=0.5) +
    facet_wrap(c("kingdom"), scales="free_y") +
    stat_compare_means(
      symnum.args = symnum.args,
      method = "wilcox",
      comparisons = list(c("Generalist", "Specialist"))) +
    ylab("Number of CDS") +
    guides(color = "none") +
    scale_color_prevalence_type() +
    scale_y_log10() +
    theme_bw() +
    annotation_logticks(sides = "l") +
    theme(axis.title.x = element_blank(),
          strip.text = element_text(size = rel(1.1)),
          axis.text = element_text(size = rel(1.1)),
          legend.position = "none",
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))

    
  c <- generalist_specialists_bgcs %>% 
    left_join(genome_length %>% 
                rename(accession = genome_id) %>% 
                unnest(g_size) %>% 
                mutate(g_size = as.numeric(g_size)) %>% 
                mutate(g_size = g_size/10^6) %>% 
                distinct(accession, g_size)) %>% 
    rename(Groups = type) %>% 
    mutate(n = n/g_size) %>% 
    group_by(accession, Groups, kingdom) %>% 
    summarise(value = sum(n)) %>%
    ggplot(aes(Groups, value, color = Groups)) +
    geom_boxplot(width=0.5) +
    scale_y_log10() +
    stat_compare_means(
      symnum.args = symnum.args,
      method = "wilcox", 
      comparisons = list(c("Specialist", "Generalist"))) +
    facet_wrap(c("kingdom"), scales="free_y") +
    scale_color_prevalence_type() +
    annotation_logticks(sides = "l") +
    ylab("Number of BCGs/Genome Length") +
    guides(color = FALSE) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          strip.text = element_text(size = rel(1.1)),
          axis.text = element_text(size = rel(1.1)),
          legend.position = "none",
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
  
  amr <- valid_amr %>% 
    rename(accession = genome_id) %>% 
    left_join(genomes2find_BGCs %>% 
                distinct(accession, genus, type, kingdom))
  
  amr2 <- amr %>% 
    rename(Groups = type,
           element_type = `Element type`,
           element_sub = `Element subtype`,
           class = Class) %>% 
    select(accession, Groups, kingdom, element_type, element_sub, class) %>% 
    mutate(count = 1) %>% 
    group_by(accession, Groups, kingdom, element_type, element_sub, class) %>% 
    summarise(n = sum(count)) %>% 
    left_join(genome_length %>% 
                rename(accession = genome_id) %>% 
                unnest(g_size) %>% 
                mutate(g_size = as.numeric(g_size)) %>% 
                mutate(g_size = g_size/10^6) %>% 
                distinct(accession, g_size)) %>% 
    mutate(norm_amr = n/g_size)
  
  symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                      symbols = c("****", "***", "**", "*", "ns"))
  
  d <- amr2 %>% 
    filter(element_type == "AMR") %>% 
    ggplot(aes(Groups, norm_amr, color = Groups)) +
    geom_point(size = 0.1) +
    geom_boxplot(width=0.5) +
    stat_compare_means(
      method = "wilcox", 
      symnum.args = symnum.args,
      comparisons = list(c("Generalist", "Specialist")),
      method.args = list(alternative = "less")) +
    facet_wrap(c("element_type")) +
    scale_y_log10() +
    scale_y_continuous(trans='log10') +
    scale_color_prevalence_type() +
    annotation_logticks(sides = "l") +
    ylab("Number of AMR genes/Genome length") +
    guides(color = FALSE) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          strip.text = element_text(size = rel(1.1)),
          axis.text = element_text(size = rel(1.1)),
          legend.position = "none",
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) 
  
  generalists_common_coabundance_graph <- readRDS("cache/15-microbiome-rserver-hki-rdata/generalists_common_coabundance_graph.rds")
  
  kingdoms_colors <- c(
    "Bacteria" = "#3f78a9",
    "Fungi" = "#A93F55",
    "all" = "black"
  )
  
  
  scale_color_kingdom <- function(drop = FALSE, na.value = "#d9d9d9", ...) {
    scale_color_manual(
      values =  purrr::simplify(kingdoms_colors),
      na.value = na.value,
      drop = drop,
      ...
    )
  }
  
  environment_group_colors <- c(
    "soil" = "#8b786d",
    "aquatic" = "#78a1bb",
    "host" = "#e0aa5a",
    "all" = "black"
  )
  
  environment_group_sets_colors <- c(
    "host-aquatic-soil" = "#000000",
    "aquatic-soil" = "#155A82",
    "host-aquatic" = "#83155F",
    "host-soil" = "#826516"
  ) %>%
    append(environment_group_colors)
  
e <- 
    generalists_common_coabundance_graph %>%
    activate(edges) %>%
    mutate(estimate_direction = estimate_direction %>% recode("positive" = "Positive", "negative" = "Negative")) %>%
    ggraph(layout = layout_in_circle(generalists_common_coabundance_graph, order = order(V(generalists_common_coabundance_graph)$taxon_group))) +
    geom_edge_link(aes(color = environment_group, alpha = involves_prevalence_group_taxon)) +
    geom_node_point(aes(color = kingdom, size = is_generalist)) +
    scale_edge_color_manual(values = environment_group_sets_colors, labels = str_to_sentence) +
    scale_color_kingdom(na.value = "black") +
    scale_size_manual(values = c(`TRUE` = 1.5, `FALSE` = 0.2), labels = c("No", "Yes")) +
    scale_edge_alpha_manual(values = c(`TRUE` = 0.5, `FALSE` = 0.05), guide = "none") +
    coord_fixed() +
    facet_edges(~estimate_direction) +
    guides(
      color = guide_legend(override.aes = list(size = 6)),
      size = guide_legend(override.aes = list(size = c(2,6)),
      edge_color = guide_legend(override.aes = list(edge_width = 13)))
           ) +
    theme_void() +
    labs(
      color = "Kingdom",
      edge_color = "Common in environments",
      size = "Generalist") +
    theme(
      legend.text = element_text(size = rel(1.1)),
      legend.title = element_text(size = rel(1.1)),
      strip.text.x = element_text(size = 12, margin = margin(b = 5))
    )
  

layout <- "
AABB
CCD#
EEEF
"
  
  f <-
    wrap_plots(
      wrap_plots(a$Bacteria + labs(tag = "A"), a$Fungi),
      b + labs(tag = "B"),
      c + labs(tag = "C"),
      d + labs(tag = "D"),
      e + labs(tag = "E"),
    guides = "collect") +
    guide_area() +
    plot_layout(design = layout)
  
  ggsave(filename = "plots/04-manuscript-ap/fig2-v1.11.pdf",
         plot = f,
         device = "pdf",
         width = 30,
         height = 35,
         units = "cm",
         dpi = 300)
  
  ggsave(filename = "plots/04-manuscript-ap/fig2-v1.11.png",
         plot = f,
         width = 30,
         height = 35,
         units = "cm",
         dpi = 300)
  
  return(f)
 
}

cds_to_plot <- function(genomes2find_BGCs, cds_by_genome, genome_length) {
  
  genomes2find_BGCs %>% 
    distinct() %>% 
    select(accession, genus, type, kingdom) %>% 
    left_join(cds_by_genome, by = c("accession" = "genome_id")) %>% 
    select(-value) %>% 
    left_join(genome_length %>% 
                rename(accession = genome_id) %>% 
                unnest(g_size) %>% 
                mutate(g_size = as.numeric(g_size)) %>% 
                mutate(g_size = g_size/10^6) %>% 
                distinct(accession, g_size)) %>% 
    mutate(cds_norm = cds_n/g_size) %>% 
    rename(Group = type)
}

permutation_test <- function(df, y = "cds_n", x = "Group", n_iter = 1000) {
  
  p_values <- vector("numeric", n_iter)
  
  # Loop through the iterations
  for (i in 1:n_iter) {
    
    # Shuffle the labels randomly
    shuffled_labels <- sample(df[[x]])
    
    # Combine the shuffled labels with the variable data
    shuffled_data <- data.frame(values = df[[y]], group = shuffled_labels)

    #Perform the Wilcoxon rank-sum test and store the p-value
    test_result <- wilcox.test(values ~ group, data = shuffled_data)
    p_values[i] <- test_result$p.value
    
  }
  
  observed_p_value <- wilcox.test(df[[y]] ~ df[[x]])$p.value

  # Calculate the proportion of permutation p-values smaller than the observed p-value
  proportion_smaller <- mean(p_values <= observed_p_value)

  # If the proportion of permuted p-values smaller than or equal to the observed 
  # p-value is 0, it suggests that none of the permuted datasets produced a result
  # as extreme as the observed result. This implies that the observed difference 
  # between the groups defined by x (as measured by the test statistic, typically 
  # the Wilcoxon rank-sum test statistic in this case) is so extreme that it is 
  # unlikely to have occurred by random chance alone.
  
  # Therefore, a proportion of 0 indicates strong evidence against the null 
  # hypothesis of no difference between the groups, implying that there are 
  # indeed differences between the groups.
  
  return(proportion_smaller)
}

bgcs_df_paper <- function(generalist_specialists_bgcs, genome_length) {
  
  generalist_specialists_bgcs %>% 
    left_join(genome_length %>% 
                rename(accession = genome_id) %>% 
                unnest(g_size) %>% 
                mutate(g_size = as.numeric(g_size)) %>% 
                mutate(g_size = g_size/10^6) %>% 
                distinct(accession, g_size)) %>% 
    rename(Groups = type) %>% 
    mutate(n = n/g_size) %>% 
    group_by(accession, Groups, kingdom) %>% 
    summarise(value = sum(n))
}

amr_df_paper <- function(valid_amr, genomes2find_BGCs, genome_length) {
  
valid_amr %>% 
    rename(accession = genome_id) %>% 
    left_join(genomes2find_BGCs %>% 
                distinct(accession, genus, type, kingdom)) %>% 
    rename(Groups = type,
           element_type = `Element type`,
           element_sub = `Element subtype`,
           class = Class) %>% 
    select(accession, Groups, kingdom, element_type, element_sub, class) %>% 
    mutate(count = 1) %>% 
    group_by(accession, Groups, kingdom, element_type, element_sub, class) %>% 
    summarise(n = sum(count)) %>% 
    left_join(genome_length %>% 
                rename(accession = genome_id) %>% 
                unnest(g_size) %>% 
                mutate(g_size = as.numeric(g_size)) %>% 
                mutate(g_size = g_size/10^6) %>% 
                distinct(accession, g_size)) %>% 
    mutate(norm_amr = n/g_size) %>% 
    filter(element_type == "AMR")
}


compare_bcgs_products_bac_gen <- function(generalist_specialists_bgcs, genome_length, bcgs_groups) {
  
  bcgs_df <- 
    generalist_specialists_bgcs %>% 
    left_join(genome_length %>% 
                rename(accession = genome_id) %>% 
                unnest(g_size) %>% 
                mutate(g_size = as.numeric(g_size)) %>% 
                mutate(g_size = g_size/10^6) %>% 
                distinct(accession, g_size)) %>% 
    rename(Groups = type) %>% 
    mutate(n = n/g_size) %>% 
    mutate(products = replace_na(products, "No-BCGs")) %>% 
    group_by(accession, Groups, kingdom, products) %>% 
    summarise(value = sum(n)) %>% 
    filter(kingdom == "Bacteria" & value > 0) %>% 
    left_join(bcgs_groups)
  
# p_products <- 
  bcgs_df %>% 
  ggplot(aes(Groups, value)) +
  geom_boxplot(aes(color=cluster)) +
  stat_compare_means(
    method = "wilcox",
    symnum.args = symnum.args,
    label = "..p.adj..",
    comparisons = list(c("Specialist", "Generalist")),
    method.args = list(alternative = "greater"),
  ) +
  #facet_wrap(~Groups, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_y_log10()

# If the p-value is less than a chosen significance level (commonly 0.05), it 
# suggests that there is statistically significant evidence to reject the null 
# hypothesis that the means are equal, in favor of the alternative hypothesis 
# that the mean of the "Generalist" group is greater than the mean of the 
# "Specialist" group.

ggsave(filename = "plots/04-manuscript-ap/bcgs_products_bacteria.pdf", 
       plot = p_products, width = 9, height = 7, dpi = 300)

ggsave(filename = "plots/04-manuscript-ap/bcgs_products_bacteria.png", 
       plot = p_products, width = 9, height = 7, dpi = 300)

}
