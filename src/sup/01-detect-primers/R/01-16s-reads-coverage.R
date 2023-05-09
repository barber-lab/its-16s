library(tidyverse)
#/home/ailtonpcf/draco/proj/04-global-microbiome/
# Retrieve reference sequence and read aligned to the reference
bed <- list.files(path = "cache/10-bam2bed", 
                  full.names = T) %>% 
  as_tibble() %>% 
  # slice(1) %>% 
  mutate(sample_id = str_split(value, "/", simplify=T)[,3]) %>% 
  mutate(sample_id = str_remove(sample_id, ".bed")) %>% 
  mutate(import_bed = map(
    .x = value, 
    .f = ~ read_tsv(file = .x, col_names = F) %>% 
      rename(reference = X1,
             read = X4) %>% 
      select(reference, read))) %>% 
  select(-value)

# Retrieve reads length
read_length <- list.files(path = "cache/09-reads-length", 
                          full.names = T) %>% 
  as_tibble() %>% 
  # slice(1) %>% 
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
  select(-value)

# Retrieve the coverage stats from the alignment
# Most important information: Covered_bases (Alignment length between reference and read)
coverage <- list.files(path = "cache/08-coverage", 
                  full.names = T) %>% 
  as_tibble() %>% 
  # slice(1) %>% 
  mutate(sample_id = str_split(value, "/", simplify=T)[,3]) %>% 
  mutate(sample_id = str_remove(sample_id, ".txt")) %>% 
  mutate(import_cov = map(
    .x = value, 
    .f = ~ read_tsv(file = .x, col_names = T))) %>% 
  select(-value)

# coverage %>% 
#   unnest(import_cov)

merge_metrics <- bed %>% 
  unnest(import_bed) %>%
  left_join(
    read_length %>% 
      unnest(import_length)
  ) %>% 
  left_join(
    coverage %>% 
      unnest(import_cov) %>% 
      rename(reference = "#ID") %>% 
      select(sample_id, reference, Covered_bases)
  ) %>% 
  mutate(reads_cov = Covered_bases/read_length*100) %>% 
  mutate(reference = str_remove(reference, "_[:digit:]+")) %>% 
  mutate(reference = str_remove(reference, "_amp")) 

runs <- read_csv(file = "raw/runs.csv")

proj16s <- merge_metrics %>%
  left_join(runs %>% 
              select(run_id, bioproject_id), by = c("sample_id" = "run_id"))

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

ggsave(filename = "plots/01-detect-primers/bac_16s_reads_coverage2.png",
       plot = primers,
       device = "png", 
       width = 12, 
       height = 9, 
       dpi = 150)  

