#!/usr/bin/env R

list_files <- function(path = path) {
  
  list.files(path = path, 
             full.names = T) %>% 
    as_tibble() 
  # %>% 
  #   slice(1:2)
}

mapped_reads_from_bed <- function(tbl) {
  
  tbl %>% 
    mutate(sample_id = str_split(value, "/", simplify=T)[,3]) %>%
    mutate(sample_id = str_remove(sample_id, ".bed")) %>%
    mutate(import_bed = map(
      .x = value,
      .f = ~ read_tsv(file = .x, col_names = F) %>%
        rename(reference = X1,
               read = X4) %>%
        select(reference, read))) %>%
    select(-value) %>% 
    unnest(import_bed)
}

import_reads_length <- function(tbl) {
  
  tbl %>% 
    mutate(sample_id = str_split(value, "/", simplify=T)[,3]) %>% 
    mutate(sample_id = str_remove(sample_id, ".txt")) %>% 
    mutate(import_length = map(
      .x = value, 
      .f = ~ read_tsv(file = .x, col_names = F) %>% 
        separate(col = X1, 
                 into = c("read", "other1", "other2"), 
                 sep = " ") %>% 
        select(read, X2) %>% 
        rename(read_length = X2))) %>% 
    select(-value) %>% 
    unnest(import_length)
  
}

aln_covered_bases <- function(tbl) {
  
  tbl %>% 
    mutate(sample_id = str_split(value, "/", simplify=T)[,3]) %>% 
    mutate(sample_id = str_remove(sample_id, ".txt")) %>% 
    mutate(import_cov = map(
      .x = value, 
      .f = ~ read_tsv(file = .x, col_names = T))) %>% 
    select(-value) %>% 
    unnest(import_cov)
}

merge_metrics <- function(import_bed, get_reads_length, import_aln_coverage, runs){
  
  merge_metrics <- import_bed %>% 
    left_join(get_reads_length) %>% 
    left_join(
      import_aln_coverage %>% 
        rename(reference = "#ID") %>% 
        select(sample_id, reference, Covered_bases)
    ) %>% 
    mutate(reads_cov = Covered_bases/read_length*100) %>% 
    mutate(reference = str_remove(reference, "_[:digit:]+")) %>% 
    mutate(reference = str_remove(reference, "_amp"))
  
  merge_metrics %>%
    left_join(runs %>% 
                select(run_id, bioproject_id), by = c("sample_id" = "run_id"))
  
}

plot_reads_coverage <- function(proj16s) {
  
  rm_regions <- proj16s %>%
    rowid_to_column("index") %>% 
    group_by(bioproject_id, reference) %>% 
    count() %>% 
    filter(n < 100) %>% 
    select(-n) %>% 
    mutate(rm = TRUE)
  
  proj16s_fltrd <- proj16s %>%
    group_by(bioproject_id, reference) %>% 
    nest() %>%
    left_join(rm_regions) %>%
    filter(is.na(rm)) %>%
    select(-rm) %>% 
    unnest()
  
  primers <- proj16s_fltrd %>% 
    ggplot(aes(x=reads_cov)) +
    geom_histogram(mapping = aes(fill=reference), 
                   alpha = 0.8) +
    facet_wrap(c("bioproject_id"), scales = "free_y") +
    xlab("Reads coverage (%)") +
    ylab("Count") +
    xlim(c(0,200))
  
  ggsave(filename = "plots/01-detect-primers/bac_16s_reads_coverage_with_NAs.png",
         plot = primers,
         device = "png", 
         width = 12, 
         height = 9, 
         dpi = 150)  
}

plot_reference_coverage <- function(import_aln_coverage, runs) {
  
  proj16s <- import_aln_coverage %>%
    rename(reference = "#ID") %>%
    left_join(runs %>% select(run_id, bioproject_id), by = c("sample_id" = "run_id")) %>%
    relocate(sample_id, bioproject_id) %>%
    mutate(reference = str_remove(reference, "_[:digit:]+")) %>%
    mutate(reference = str_remove(reference, "_amp"))
  
  primers <- proj16s %>%
    ggplot(aes(x=Covered_percent)) +
    geom_histogram(aes(fill=reference)) +
    facet_wrap(c("bioproject_id"), scales = "free_y") +
    xlab("Coverage (%)") +
    ylab("Count")
  
  ggsave(filename = "plots/01-detect-primers/bac_16s_reference_coverage_NAs.png",
         plot = primers,
         device = "png",
         width = 12,
         height = 9,
         dpi = 150)

}

extracted_its <- function(path, runs) {
  
  its <- read_delim(file = path, 
                    col_names = F)
  
  # Parse fasta headers from grep
  seqs <- its %>%
    select(X1) %>% 
    mutate(sample_id = str_split(X1, "\\.", simplify=T)[,1]) %>% 
    mutate(ITS = str_extract(X1, "ITS[1-2]"))
  
  # Merge run_id and bioprojects
  sampleByITS <- seqs %>% 
    drop_na(ITS) %>% 
    group_by(sample_id, ITS) %>%
    count() %>% 
    #pivot_wider(names_from = its, values_from = n, values_fill = 0) %>% 
    left_join(runs %>% select(run_id, bioproject_id), by = c("sample_id" = "run_id"))
  
  # Primers used
  primers <- sampleByITS %>% 
    ggplot(aes(x=n)) +
    geom_histogram(aes(fill=ITS)) +
    facet_wrap(c("bioproject_id"), scales = "free") +
    xlab("Sequences number") +
    ylab("Count")
  
  # Plot
  ggsave(filename = "plots/01-detect-primers/fungi_its_amplified_NAs.png",
         device = "png", 
         width = 12, 
         height = 9, 
         dpi = 150)
}