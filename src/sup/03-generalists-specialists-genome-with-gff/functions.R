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
  # %>% 
  # Remove taxa with no antismash hits
  #drop_na(products)
  
  
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
    #unite(col = Categories, c("kingdom", "type"), sep = " ") %>% 
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
  #In order to avoid horizontal white line on geom_col, summarize it first!
  
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
      #symnum.args = symnum.args,
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