library(tidyverse)
library(phyloseq)
library(furrr)
library(decontam)
library(ggpubr)
library(patchwork)
library(rstatix)

# Taxonomic ranks 
groups <- list(
  "species" = c("kingdom", "phylum", "class", "order", "family", "genus", "species")
)

control_pattern <- c("ExtNeg|PCRNeg|bactexControl|waterCtrl|poolextractCtrl|poolPCRwaterCtrl|NEG|bactextractCtrl")
fix_proj <- c("PRJNA522449" = "PRJNA241408_PRJNA522449",
              "PRJNA241408" = "PRJNA241408_PRJNA522449")

unite_ids <- aqua_analysed %>% 
  filter(amplicon == "unite") %>% 
  pull(samples)

control_metadata <- list.files("raw/06-aqua-projects-with-controls", full.names = T) %>% 
  as_tibble() %>% 
  mutate(meta = map(value, ~ read_csv(.x)  %>%
                      select(any_of(c("Run", "LibraryLayout", "BioProject", "Library Name", "Sample Name")))
  )) %>% 
  unnest(meta) %>% 
  mutate(BioProject = recode(BioProject, !!!fix_proj)) %>% 
  rename(sample_name = "Sample Name", libnName = "Library Name") %>% 
  mutate(ctrl_type = map(libnName, ~ str_extract(.x, control_pattern))) %>%
  unnest(ctrl_type) %>% 
  mutate(ctrl_type = replace_na(ctrl_type, "sample")) %>% 
  filter(Run  %in% unite_ids) %>% 
  rename(samples = Run)

control_samples <- control_metadata %>% 
  filter(!ctrl_type == "sample") %>% 
  select(samples, BioProject, ctrl_type)

real_samples <- control_metadata %>% 
  filter(ctrl_type == "sample") %>% 
  select(samples, BioProject, ctrl_type)

unite_tax <- read_tsv("ref/04-silva132_97-unite8.2_97/unite_tax_original.txt", col_names = F) %>% 
  mutate(X2 = str_remove_all(X2, "[a-z]__")) %>% 
  separate(col = X2, into = groups$species, sep = ";") %>% 
  rename(otu = X1)

its_res <- list.files(path = "cache/08-decontamination/19-abundances", 
                                 full.names = T, pattern = "unite") %>% 
  as_tibble() %>% 
  mutate(data = map(value, ~ read_tsv(.x, skip = 1))) %>% 
  unnest(data) %>% 
  select(-value)

# Add taxonomy to abundance
its_tax <- its_res %>% 
  rename(otu = "#OTU ID") %>% 
  pivot_longer(cols = -otu, names_to = "samples", values_to = "abundance") %>% 
  filter(abundance > 0) %>% 
  left_join(unite_tax)


its_control <- control_samples %>% 
  left_join(its_tax) %>% 
  drop_na(otu)

its_real <- real_samples %>% 
  left_join(its_tax) %>% 
  drop_na(otu)

#################### RUN DECONTAMINATION BY BIOPROJECT AND CONTROL TYPE #########################

tax <- its_control %>% 
  group_by(BioProject, ctrl_type) %>% 
  nest(.key = "tax") %>% 
  mutate(tax = map2(tax, BioProject, ~ .x %>% 
                      bind_rows(its_real)
                    )) 

tax$tax[[1]] %>% 
  filter(str_detect(genus, "Cortinarius"))

otu <- its_control %>% 
  group_by(BioProject, ctrl_type) %>% 
  nest(.key = "otu") %>% 
  mutate(otu = map2(otu, BioProject, ~ .x %>% 
                      bind_rows(its_real) %>% 
                      filter(BioProject == .y | is.na(BioProject))  %>%
                      select(samples, otu, abundance) %>% 
                      pivot_wider(id_cols = otu, 
                                  names_from = samples, 
                                  values_from = abundance, 
                                  values_fill = 0) %>% 
                      column_to_rownames("otu") %>%
                      as.matrix() %>%
                      otu_table(taxa_are_rows = TRUE)
                    ))

meta <- its_control %>% 
  group_by(BioProject, ctrl_type) %>% 
  select(samples, BioProject, ctrl_type) %>% 
  nest(.key = "meta") %>% 
  mutate(meta = map2(meta, BioProject, ~ .x %>% 
                       bind_rows(its_real) %>%
                       mutate(BioProject = replace_na(BioProject, .y)) %>% 
                       mutate(ctrl_type = replace_na(ctrl_type, "control")) %>% 
                       filter(BioProject == .y | is.na(BioProject)) %>%
                       select(BioProject, samples, ctrl_type) %>% 
                       distinct() %>% 
                       as.data.frame() %>%
                       column_to_rownames("samples") %>%
                       sample_data()
  )) 

detect_contaminants <- function(ps, thrsld) {
  
  sample_data(ps)$is.neg <- sample_data(ps)$ctrl_type == "control"
  contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=thrsld)
} 

plan(multisession, workers = 5)

df_phylo <- meta %>% 
  left_join(tax) %>% 
  left_join(otu) %>% 
  mutate(phylo = future_pmap(.l = list(otu, tax, meta), 
                             .f = ~ phyloseq(..1, ..2, ..3))) %>% 
  expand_grid(threshold = c(0.1, 0.5)) %>% 
  mutate(cont = future_map2(phylo, threshold, ~ detect_contaminants(.x, .y)))

generalist_specialist <- read_tsv("raw/generalists/generalist.txt") %>%
  rename(genus = generalist) %>%
  mutate(type = "Generalist") %>%
  bind_rows(
    read_tsv("raw/specialists/specialists.txt") %>%
      rename(genus = specialists) %>%
      mutate(type = "Specialist"))

plts_dir <- "cache/08-decontamination/20-plots"

dir.create(path = plts_dir, recursive = T)

summarise_by_contaminant <- function(otu, cont_res) {
  otu %>% 
    left_join(
      cont_res %>% 
        rownames_to_column("otu")
    ) 
}

genera_special_contamination <- function(df) {
  
  # Calculate tss for all
  not_decon <- df %>% 
    group_by(samples, genus, contaminant) %>%
    summarise(abundance = sum(abundance)) %>% 
    group_by(samples) %>% 
    mutate(tss = abundance/sum(abundance)) %>% 
    mutate(group = "Not decontaminated") %>% 
    left_join(generalist_specialist) %>%
    filter(type == "Generalist") %>% 
    drop_na(type)
  
  # Calculate tss for non contaminants
  decon <- df %>% 
    group_by(samples, genus, contaminant) %>%
    summarise(abundance = sum(abundance)) %>% 
    group_by(samples) %>% 
    mutate(tss = abundance/sum(abundance)) %>% 
    mutate(group = "Decontaminated") %>% 
    left_join(generalist_specialist) %>%
    filter(type == "Generalist") %>% 
    drop_na(type) %>% 
    filter(!contaminant) # Remove after normalization
  
  # summarise by contamination status
  bind_rows(decon, not_decon)  %>% 
    arrange(samples, genus) %>% 
    group_by(samples, genus, type, group) %>% 
    summarise(tss = sum(tss)) %>% 
    arrange(samples, genus)
}

plot_generalist_genus <- function(df, proj, control, threshold) {
  
  comp_16s <- df %>%
    ungroup() %>%
    distinct(genus, group) %>%
    mutate(label = str_glue("{genus} {group}")) %>%
    pivot_wider(id_cols = genus, names_from = group, values_from = label) %>%
    rename(not_decon = "Not decontaminated") %>%
    mutate(comp = map2(Decontaminated, not_decon, ~ c(.x, .y))) %>%
    pull(comp)
  
  p <- df %>% 
    mutate(paired = str_glue("{samples}_{genus}")) %>% 
    mutate(genus = str_glue("{genus} {group}")) %>%
    ggplot(aes(genus, tss)) +
    geom_point(aes(genus, tss), alpha = 0.5) +
    geom_line(aes(group=paired), alpha = 0.5, colour = "darkgrey") +
    geom_boxplot(aes(fill=group)) +
    stat_compare_means(aes(label=p.signif),
                       comparisons = comp_16s, 
                       step.increase = 0,
                       method = "wilcox.test") +
    scale_y_log10() +
    coord_flip() +
    labs(title = proj, subtitle = paste("Control: ",control, " Threshold: ", threshold))
  
  ggsave(filename = paste0(plts_dir, "/", proj, "_", control, "_", threshold, "_", "generalists_genus_its.png"),
         plot = p,
         width = 9,
         height = 7,
         dpi=150)
}

otu <- its_control %>%
  group_by(BioProject, ctrl_type) %>%
  nest(.key = "otu") %>%
  mutate(otu = map2(otu, BioProject, ~ .x %>%
                      bind_rows(its_real) %>%
                      filter(BioProject == .y | is.na(BioProject))
                    
  ))

phylo2 <- df_phylo %>% 
  select(BioProject, ctrl_type, cont, threshold) %>% 
  left_join(otu) %>% 
  mutate(otu_by_cont = map2(otu, cont, ~ summarise_by_contaminant(.x, .y))) %>% 
  mutate(gen_spe = map(otu_by_cont, ~ genera_special_contamination(.x))) %>% 
  mutate(plt = pmap(list(gen_spe, BioProject, ctrl_type, threshold), ~ plot_generalist_genus(..1, ..2, ..3, ..4)))

gen_spec_by_control_grouped <-
  phylo2 %>% 
  select(BioProject, ctrl_type, gen_spe, threshold) %>% 
  unnest(gen_spe) %>% 
  distinct() %>% 
  mutate(paired = str_glue("{samples}_{genus}")) %>% 
  ggplot(aes(group, tss)) +
  geom_point(alpha = 0.5)+ 
  geom_line(aes(group=paired), alpha = 0.5, colour = "darkgrey") + 
  geom_boxplot(aes(fill=group)) +
  stat_compare_means(aes(label=p.signif),comparisons = list(c("Not decontaminated", "Decontaminated"))) +
  scale_y_log10() +
  facet_wrap(c("ctrl_type", "BioProject", "threshold"))

ggsave(filename = paste0(plts_dir, "/", "gen_spec_by_control_its.png"),
       plot = gen_spec_by_control_grouped,
       width = 9,
       height = 7,
       dpi=150)

# Figure for the paper
###############################################################################

groups_fill <- c(
  "Contaminated"   = "#d8b365",
  "Decontaminated" = "#5ab4ac"
)

contaminants_by_proj <- function(df, proj, boxp_wid, x_margin, dot_size, habitat) {
  
  comp_16s <- df %>%
    ungroup() %>%
    distinct(ctrl_type, group) %>%
    mutate(label = str_glue("{ctrl_type} {group}")) %>%
    pivot_wider(id_cols = ctrl_type, names_from = group, values_from = label) %>%
    mutate(comp = map2(Decontaminated, Contaminated, ~ c(.x, .y))) %>%
    pull(comp)
  
  xlabel <- df %>%
    ungroup() %>%
    distinct(ctrl_type, group) %>%
    mutate(label = str_glue("{ctrl_type} {group}")) %>%
    mutate(label = map2(label, ctrl_type, ~ case_when(
      str_detect(.x, "Decontaminated") ~ "",
      str_detect(.x, "Contaminated") ~ .y))) %>% 
    unnest(label) %>% 
    pull(label)
  
  df %>% 
    mutate(paired = str_glue("{samples}_{genus}_{ctrl_type}")) %>% 
    mutate(ctrl_type = str_glue("{ctrl_type} {group}")) %>%
    ggplot(aes(ctrl_type, tss)) +
    geom_line(aes(group=paired), alpha = 0.5, colour = "darkgrey") +
    geom_point(alpha = 0.5, size = dot_size) +
    geom_boxplot(aes(fill=group), 
                 position = position_dodge2(preserve = "single", width = boxp_wid), 
                 outlier.size = dot_size) +
    stat_compare_means(aes(label=p.signif),comparisons = comp_16s, step.increase = 0) +
    labs(x = "", y = "TSS", fill = "Groups", title = proj, subtitle = habitat) +
    scale_y_continuous(expand = expansion(mult = c(0, 1000))) +
    scale_y_log10() +
    theme_bw() +
    scale_x_discrete(labels = xlabel) +
    theme(
      strip.text = element_text(size = rel(1.1)),
      axis.text = element_text(size = rel(1.1)), 
      axis.text.x = element_text(angle = 90, hjust = 0, vjust = x_margin),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_blank(),
      axis.ticks.x = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA))
}

param1 <- c(
  "BioProject" = "PRJNA282687",
  "boxp_wid" = 0.5,
  "x_margin" = -0.7,
  "dot_size" = 0.5,
  "habitat" = "Aquatic: Large lakes"
)

param2 <- c(
  "BioProject" = "PRJEB30970",
  "boxp_wid" = 0.5,
  "x_margin" = 6,
  "dot_size" = 0.5,
  "habitat" = "Aquatic: Other freshwater"
)

params <- paste("param", seq(1:2), sep = "") %>% 
  str_split(" ", simplify = TRUE) %>% 
  map_dfr(~eval(parse(text = .))) %>% 
  mutate_at(.vars = c("boxp_wid", "x_margin", "dot_size"), .funs = as.numeric)

proj_its_general <- phylo2 %>% 
  select(BioProject, ctrl_type, gen_spe, threshold) %>% 
  filter(threshold == 0.1) %>% 
  unnest(gen_spe) %>% 
  mutate(group = recode(group, "Not decontaminated" = "Contaminated")) %>% 
  group_by(BioProject) %>% 
  nest() %>% 
  left_join(params) 

proj_its_general_70 <- df_its %>% 
  mutate(ctrl_type = "DNA concentration") %>%
  select(BioProject, ctrl_type, df_decont, thrshld) %>% 
  filter(thrshld == 0.1) %>% 
  unnest(df_decont) %>% 
  mutate(group = recode(group, "Not decontaminated" = "Contaminated", "decontaminated" = "Decontaminated")) %>% 
  group_by(BioProject) %>% 
  nest() %>% 
  left_join(params) 


plt_decontaminated <- function(data, proj, habitat, boxp_wid, dot_size) {
  
  data %>%
    mutate(paired = str_glue("{samples}_{genus}")) %>% 
    mutate(group = recode(group, "Contaminated" = "Raw data")) %>% 
    mutate(group = factor(group, levels = c("Raw data", "Decontaminated"))) %>% 
    ggplot(aes(group, tss)) +
    geom_line(aes(group=paired), alpha = 0.5, colour = "darkgrey", size = 0.1) + 
    geom_point(alpha=0.5, size = dot_size)+ 
    geom_boxplot(aes(fill=group), 
                 position = position_dodge2(preserve = "single", width = boxp_wid), 
                 width = boxp_wid,
                 outlier.size = dot_size) +
    stat_compare_means(aes(label=p.signif), comparisons = list(c("Raw data", "Decontaminated"))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.5))) +
    labs(x = "", y = "TSS", fill = "Groups", title = proj, subtitle = habitat) +
    scale_y_log10() +
    facet_wrap(c("ctrl_type")) +
    theme_bw() +
    theme(
      strip.text = element_text(size = rel(1.1)),
      axis.text = element_text(size = rel(1.1)), 
      axis.text.x = element_blank(),
      legend.title = element_text(size = rel(1.3)), 
      legend.text = element_text(size = rel(1.3)), 
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(face = "bold"),
      panel.border = element_rect(colour = "black", fill = NA))
  
}

plt_its <- proj_its_general %>% 
  bind_rows(proj_its_general_70) %>% 
  mutate(plt_proj = pmap(list(data, BioProject, habitat, boxp_wid, dot_size), ~ plt_decontaminated(.x, .y, ..3, ..4, ..5)))
