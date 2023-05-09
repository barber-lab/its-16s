#!/usr/bin/env R

list_amr_files <- function(path) {
  
  list.files(path = path, full.names = T, pattern = ".tsv$", recursive = T) %>% 
    as_tibble() %>% 
    mutate(genome_id = str_extract(value, "GCA_\\d+\\.\\d+")) 
  
}

parse_amr <- function(path) {
  
  df <- read_tsv(path)
  
  if (nrow(df != 0)) {
    return(df)
  } else {
    return(NA)
  }
  
}

import_amr_results <- function(amr_results_files) {
  
  amr_results_files %>%
    mutate(amr = map(value, ~ parse_amr(.x)))
  
}

tidy_amr <- function(get_amr) {
  
  get_amr %>% 
    filter(!map_lgl(amr, is_logical)) %>% 
    unnest(amr) %>% 
    mutate(genome_id = str_extract(value, "G.._\\d+\\.\\d+")) %>% 
    mutate(`Element type` = recode(`Element type`, "STRESS" = "Stress"))
}

scale_color_prevalence_type <- function(drop = FALSE, ...) {
  scale_color_manual(
    values =  purrr::simplify(prevalence_type),
    na.value = "#000000",
    drop = drop,
    ...
  )
}

amr_plots <- function(valid_amr, genomes2find_BGCs) {
  
  df <- valid_amr %>% 
    rename(accession = genome_id) %>% 
    left_join(genomes2find_BGCs %>% 
                distinct(accession, genus, type, kingdom))
  
  df2 <- df %>% 
    rename(Groups = type,
           element_type = `Element type`,
           element_sub = `Element subtype`,
           class = Class) %>% 
    select(accession, Groups, kingdom, element_type, element_sub, class) %>% 
    mutate(count = 1) %>% 
    group_by(accession, Groups, kingdom, element_type, element_sub, class) %>% 
    summarise(n = sum(count))
  
  # Type 
  ############
  
  p <- df2 %>% 
    filter(!element_type == "VIRULENCE") %>% 
    ggplot(aes(Groups, n, color = Groups)) +
    geom_boxplot(width=0.5) +
    geom_point() +
    stat_compare_means(
      method = "wilcox", 
      comparisons = list(c("Specialist", "Generalist"))) +
    scale_y_log10() +
    facet_wrap(c("element_type"), scales="free_y") +
    scale_color_prevalence_type() +
    annotation_logticks(sides = "l") +
    ylab("Number of AMR genes") +
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

  # Subtype  
  ############
  
  q <- df2 %>% 
    ggplot(aes(Groups, n, color = Groups)) +
    geom_boxplot(width=0.5) +
    geom_point() +
    stat_compare_means(
      method = "wilcox", 
      comparisons = list(c("Specialist", "Generalist"))) +
    scale_y_log10() +
    facet_wrap(c("element_sub"), scales="free_y") +
    scale_color_prevalence_type() +
    annotation_logticks(sides = "l") +
    ylab("Number of AMR genes") +
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
  
  
  ggsave("plots/07-amr/comp-amr-elem-type.png",
         plot = p, 
         device = "png", 
         width = 9, 
         height = 7, 
         dpi = 150)
  
  ggsave("plots/07-amr/comp-amr-elem-type.pdf",
         plot = p, 
         device = "pdf", 
         width = 9, 
         height = 7, 
         dpi = 150)
  
  ggsave("plots/07-amr/comp-amr-elem-subtype.png",
         plot = q, 
         device = "png", 
         width = 9, 
         height = 7, 
         dpi = 150)
  
  ggsave("plots/07-amr/comp-amr-elem-subtype.pdf",
         plot = q, 
         device = "pdf", 
         width = 9, 
         height = 7, 
         dpi = 150)
  
# Class  
############
  df3 <- df %>% 
    rename(Groups = type,
           element_type = `Element type`,
           element_sub = `Element subtype`,
           class = Class) %>% 
    select(accession, Groups, kingdom, class) %>% 
    mutate(count = 1) %>% 
    group_by(accession, Groups, kingdom, class) %>% 
    summarise(n = sum(count)) %>% 
    drop_na(class)
  
  r <- df3 %>% 
    ggplot(aes(Groups, n, color = Groups)) +
    geom_boxplot(width=0.5) +
    geom_point() +
    stat_compare_means(
      method = "wilcox", 
      comparisons = list(c("Specialist", "Generalist"))) +
    scale_y_log10() +
    facet_wrap(c("class"), scales="free_y") +
    scale_color_prevalence_type() +
    annotation_logticks(sides = "l") +
    ylab("Number of AMR genes") +
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
  
  ggsave("plots/07-amr/comp-amr-class.png",
         plot = r, 
         device = "png", 
         width = 20, 
         height = 20, 
         dpi = 150)
  
  ggsave("plots/07-amr/comp-amr-class.png.pdf",
         plot = r, 
         device = "pdf", 
         width = 20, 
         height = 20, 
         dpi = 150)
  
list(p,q,r)
}

export_amr <- function(valid_amr) {
  
  dir.create("tbl/07-amr", recursive = T)
  
  valid_amr %>% 
    select(genome_id, `Element type`) %>% 
    filter(!`Element type` == "VIRULENCE") %>% 
    mutate(n = 1) %>% 
    group_by(genome_id, `Element type`) %>% 
    summarise(sum_n = sum(n)) %>% 
    pivot_wider(names_from = `Element type`, values_from = sum_n, values_fill = 0) %>% 
    write_tsv("tbl/07-amr/amr_res.tsv")
}
