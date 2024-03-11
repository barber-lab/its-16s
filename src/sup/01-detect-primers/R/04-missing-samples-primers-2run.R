library(tidyverse)

## After the first round of analysis to discover which set of primers were used
## for bacteria and fungi, some projects with papers already published didn't 
## tell us neither the primer neither the region. This code is to run missing
## bioprojects, but keeping previous results.

# Import projects metadata done manually
source("src/01-detect-primers/R/03-16s-its-primers-bioprojects.R")

# Bioprojects analysed that have no reference published
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

# 1580 samples described in the paper as passing qc
samples <- read_csv("raw/runs.csv")

# Get bioprojects without primers description
bac16s <- primers %>% 
  filter(is.na(primer_bac)) %>% 
  pull(bioproject_id)

# Get bioprojects without primers description
its <- primers %>% 
  filter(is.na(primer_euk)) %>% 
  pull(bioproject_id)

bac_pd <- samples %>% 
  filter(bioproject_id  %in% c(proj_noRef, bac16s)) %>%
  filter(str_detect(amplicon, "bac")) %>% 
  filter(str_detect(library_layout, "PAIRED")) %>% 
  distinct(run_id) %>%
  write_tsv("raw/proj_no_literature_incomplete/bacteria_paired.tsv")

bac_sg <- samples %>% 
  filter(bioproject_id  %in% c(proj_noRef, bac16s)) %>%
  filter(str_detect(amplicon, "bac")) %>% 
  filter(str_detect(library_layout, "SINGLE")) %>% 
  distinct(run_id) %>%
  write_tsv("raw/proj_no_literature_incomplete/bacteria_single.tsv")

fun_pd <- samples %>% 
  filter(bioproject_id  %in% c(proj_noRef, its)) %>%
  filter(str_detect(amplicon, "fun")) %>% 
  filter(str_detect(library_layout, "PAIRED")) %>% 
  distinct(run_id) %>%
  write_tsv("raw/proj_no_literature_incomplete/fungi_paired.tsv")

fun_sg <- samples %>% 
  filter(bioproject_id  %in% c(proj_noRef, its)) %>%
  filter(str_detect(amplicon, "fun")) %>% 
  filter(str_detect(library_layout, "SINGLE")) %>% 
  distinct(run_id) %>%
  write_tsv("raw/proj_no_literature_incomplete/fungi_single.tsv")
