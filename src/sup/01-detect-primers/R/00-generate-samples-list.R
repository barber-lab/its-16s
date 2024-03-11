library(tidyverse)

runs <- read_csv("raw/runs.csv") %>% 
  as_tibble()
bioprojects <- read_csv("raw/bioprojects.csv") %>% 
  as_tibble()
biosamples <- read_tsv("raw/biosamples.tsv") %>% 
  as_tibble()

proj_noRef <- c(
  "PRJNA282687",
  "PRJNA324410",
  "PRJNA271113",
  "PRJNA359237",
  "PRJNA473079",
  "PRJNA263505",
  "PRJNA406830",
  "PRJNA432446",
  "PRJNA496065"
)

runs %>%
  filter(bioproject_id  %in% proj_noRef) %>% 
  filter(str_detect(amplicon, "bac_16s")) %>% 
  filter(str_detect(library_layout, "SINGLE")) %>%
  select(run_id) %>% 
  write_tsv(file = "draco/proj/04-global-microbiome/raw/proj_no_literature/bacteria_single.tsv",col_names = T)

runs %>%
  filter(bioproject_id  %in% proj_noRef) %>% 
  filter(str_detect(amplicon,"bac_16s")) %>% 
  filter(str_detect(library_layout, "PAIRED")) %>%
  select(run_id) %>% 
  write_tsv(file = "draco/proj/04-global-microbiome/raw/proj_no_literature/bacteria_paired.tsv",col_names = T)

runs %>%
  filter(bioproject_id  %in% proj_noRef) %>% 
  filter(str_detect(amplicon, "fun_its")) %>% 
  filter(str_detect(library_layout, "SINGLE")) %>%
  select(run_id) %>% 
  write_tsv(file = "draco/proj/04-global-microbiome/raw/proj_no_literature/fungi_single.tsv",col_names = T)

runs %>%
  filter(bioproject_id  %in% proj_noRef) %>% 
  filter(str_detect(amplicon, "fun_its")) %>% 
  filter(str_detect(library_layout, "PAIRED")) %>%
  select(run_id) %>% 
  write_tsv(file = "draco/proj/04-global-microbiome/raw/proj_no_literature/fungi_paired.tsv",col_names = T)
