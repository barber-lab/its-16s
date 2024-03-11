library(tidyverse)
library(phyloseq)
library(decontam)
library(ggpubr)

valid <- read_csv("raw/samples_used1580.csv")

groups <- list(
  "species" = c("kingdom", "phylum", "class", "order", "family", "genus", "species")
)

generalist_specialist <- read_tsv("raw/generalists/generalist.txt") %>%
  rename(genus = generalist) %>%
  mutate(type = "Generalist") %>%
  bind_rows(
    read_tsv("raw/specialists/specialists.txt") %>%
      rename(genus = specialists) %>%
      mutate(type = "Specialist"))

proj70_ids <- valid %>% 
  filter(bioproject_id == "PRJEB30970" & amplicon == "bac_16s") %>% 
  distinct(run_id) %>% 
  pull(run_id)

silva_tax <- read_tsv("ref/04-silva132_97-unite8.2_97/silva_tax.txt", col_names = F) %>% 
  mutate(X2 = str_remove_all(X2, "[A-Z]_[0-9]__")) %>% 
  separate(col = X2, into = groups$species, sep = ";") %>% 
  rename(otu = X1)

res_16s <- list.files(path = "cache/08-decontamination/19-abundances", 
                      full.names = T, pattern = "silva") %>% 
  as_tibble() %>% 
  mutate(data = map(value, ~ read_tsv(.x, skip = 1))) %>% 
  unnest(data) %>% 
  select(-value)

# Add taxonomy to abundance
tax_16s <- res_16s %>% 
  rename(otu = "#OTU ID") %>% 
  pivot_longer(cols = -otu, names_to = "samples", values_to = "abundance") %>% 
  filter(abundance > 0) %>% 
  left_join(silva_tax) %>% 
  filter(samples  %in% proj70_ids)

tax70_16s <- tax_16s %>% 
  select(-c("abundance", "samples")) %>%
  distinct() %>%
  column_to_rownames("otu") %>%
  as.matrix() %>%
  tax_table()

otu70 <- 
  tax_16s %>% 
  select(samples, otu, abundance) %>% 
  pivot_wider(id_cols = otu, 
              names_from = samples, 
              values_from = abundance, 
              values_fill = 0) %>% 
  column_to_rownames("otu") %>%
  as.matrix() %>%
  otu_table(taxa_are_rows = TRUE)

dna70 <- read_csv("raw/07-dna-concentratations-PRJEB30970/dna.csv") 

alias_pattern <- dna70 %>% 
  distinct(samples) %>% 
  pull(samples) %>% 
  paste(collapse = "|")

meta70 <- valid %>% 
  select(run_id, sample_alias, bioproject_id, amplicon) %>% 
  filter(bioproject_id == "PRJEB30970" & amplicon == "bac_16s") %>% 
  left_join(dna70 %>% rename(sample_alias = samples)) %>% 
  as.data.frame() %>%
  column_to_rownames("run_id") %>%
  sample_data()

df_by_contaminants <- function(contam70, otu_tax) {
  
  not_decontaminated <- otu_tax %>% 
    left_join(
      contam70 %>% 
        rownames_to_column("otu") %>% 
        as_tibble() %>% 
        select(-contains(c("freq", "prev")), -"p")
    ) %>% 
    mutate(contaminant = replace_na(contaminant, FALSE)) %>%
    group_by(samples, genus, contaminant) %>% 
    summarise(abundance = sum(abundance)) %>% 
    group_by(samples) %>% 
    mutate(tss = abundance/sum(abundance)) %>% 
    left_join(generalist_specialist) %>% 
    drop_na(type) %>% 
    mutate(group = "Not decontaminated")
  
  decontaminated <- otu_tax %>% 
    left_join(
      contam70 %>% 
        rownames_to_column("otu") %>% 
        as_tibble() %>% 
        select(-contains(c("freq", "prev")), -"p")
    ) %>% 
    mutate(contaminant = replace_na(contaminant, FALSE)) %>%
    group_by(samples, genus, contaminant) %>% 
    summarise(abundance = sum(abundance)) %>% 
    group_by(samples) %>% 
    mutate(tss = abundance/sum(abundance)) %>% 
    left_join(generalist_specialist) %>% 
    drop_na(type) %>%
    filter(contaminant == "FALSE") %>% 
    mutate(group = "decontaminated")
  
  bind_rows(not_decontaminated, decontaminated) %>% 
    group_by(samples, genus, group) %>% 
    summarise(tss = sum(tss))
}

plts_dir <- "cache/08-decontamination/20-plots"

plot_decont70 <- function(df, thrshld, amplicon, proj) {
  
  comp <- df %>%
    ungroup() %>%
    distinct(genus, group) %>%
    mutate(label = str_glue("{genus} {group}")) %>%
    pivot_wider(id_cols = genus, names_from = group, values_from = label) %>%
    rename(not_decon = "Not decontaminated") %>%
    mutate(comp = map2(decontaminated, not_decon, ~ c(.x, .y))) %>%
    pull(comp)
  
  p <-
    df %>% 
    mutate(paired = str_glue("{samples}_{genus}")) %>% 
    mutate(genus = str_glue("{genus} {group}")) %>%
    ggplot(aes(genus, tss)) +
    geom_point(aes(genus, tss), alpha = 0.5) +
    geom_line(aes(group=paired), alpha = 0.5, colour = "darkgrey") +
    geom_boxplot(aes(fill=group)) +
    stat_compare_means(aes(label=p.signif),comparisons = comp, step.increase = 0) +
    scale_y_log10() +
    coord_flip() +
    labs(title = proj, subtitle = paste("Threshold: ", thrshld))
  
  ggsave(filename = paste0(plts_dir, "/", "generalists_", proj, "_", thrshld, "_", amplicon, ".png"),
         plot = p,
         width = 10,
         height = 7,
         dpi=150)
  

}

df_16s_70 <- 
  tibble(thrshld = c(0.1, 0.5)) %>% 
  expand_grid(BioProject = "PRJEB30970") %>% 
  expand_grid(otu = list(otu70)) %>% 
  expand_grid(tax = list(tax70_16s)) %>% 
  expand_grid(meta = list(meta70)) %>% 
  expand_grid(otu_tax = list(tax_16s)) %>% 
  mutate(phylo = pmap(list(tax, otu, meta), ~ phyloseq(..1, ..2, ..3))) %>% 
  mutate(cont = map2(.x = phylo, 
                     .y = thrshld,
                     .f = ~ isContaminant(seqtab = .x, 
                                          method="frequency", 
                                          conc = "dna", 
                                          threshold = .y))) %>% 
  mutate(df_decont = map2(cont, otu_tax, ~ df_by_contaminants(.x, .y))) %>% 
  mutate(plt = pmap(list(df_decont,thrshld, "16s", BioProject),  ~ plot_decont70(.x, .y, ..3, ..4)))

################################################################################
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

its70_ids <- valid %>% 
  filter(bioproject_id == "PRJEB30970" & amplicon == "fun_its") %>% 
  distinct(run_id) %>% 
  pull(run_id)

its70 <- its_res %>% 
  rename(otu = "#OTU ID") %>% 
  pivot_longer(cols = -otu, names_to = "samples", values_to = "abundance") %>% 
  filter(abundance > 0) %>% 
  left_join(its_tax) %>% 
  filter(samples  %in% its70_ids)

its_tax_70 <- its70 %>% 
  select(-c("abundance", "samples")) %>%
  distinct() %>%
  column_to_rownames("otu") %>%
  as.matrix() %>%
  tax_table()

its_otu70 <- 
  its70 %>% 
  select(samples, otu, abundance) %>% 
  pivot_wider(id_cols = otu, 
              names_from = samples, 
              values_from = abundance, 
              values_fill = 0) %>% 
  column_to_rownames("otu") %>%
  as.matrix() %>%
  otu_table(taxa_are_rows = TRUE)

alias_pattern <- dna70 %>% 
  distinct(samples) %>% 
  pull(samples) %>% 
  paste(collapse = "|")

its_meta70 <- valid %>% 
  select(run_id, sample_alias, bioproject_id, amplicon) %>% 
  filter(bioproject_id == "PRJEB30970" & amplicon == "fun_its") %>% 
  left_join(dna70 %>% rename(sample_alias = samples)) %>% 
  as.data.frame() %>%
  column_to_rownames("run_id") %>%
  sample_data()

df_its <- 
  tibble(thrshld = c(0.1, 0.5)) %>% 
  expand_grid(BioProject = "PRJEB30970") %>% 
  expand_grid(otu = list(its_otu70)) %>% 
  expand_grid(tax = list(its_tax_70)) %>% 
  expand_grid(meta = list(its_meta70)) %>% 
  expand_grid(otu_tax = list(its_tax)) %>% 
  mutate(phylo = pmap(list(tax, otu, meta), ~ phyloseq(..1, ..2, ..3))) %>% 
  mutate(cont = map2(.x = phylo, 
                     .y = thrshld,
                     .f = ~ isContaminant(seqtab = .x, 
                       method="frequency", 
                       conc = "dna", 
                       threshold = .y))) %>% 
  mutate(df_decont = map2(cont, otu_tax, ~ df_by_contaminants(.x, .y))) %>% 
  mutate(plt = pmap(list(df_decont,thrshld, "its", BioProject),  ~ plot_decont70(.x, .y, ..3, ..4)))
