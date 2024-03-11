#!/usr/bin/env R

# Defaults file to bootstrap a standard global environment

# Set the correct targets project
Sys.setenv(TAR_PROJECT = "generalist_specialist_genome_metadata")

# parallelization of targets package
options(
    tidyverse.quiet = TRUE,
    clustermq.scheduler = "multiprocess"
)

library(targets)
library(tarchetypes)
library(tidyverse)
library(here)
library(jsonlite)
library(ggpubr)
library(ggsignif)
library(patchwork)

setwd(here::here())

# do not recursively call itself
# discards additional source files
# discard the targets file itself

list.files(path = "src/03-generalists-specialists-genome", 
           pattern = ".R$", 
           full.names = TRUE, 
           recursive = F) %>%
  discard(~ .x == "src/03-generalists-specialists-genome/defaults.R") %>%
  discard(~ .x == "src/03-generalists-specialists-genome/_targets.R") %>%
  discard(~ .x == "src/03-generalists-specialists-genome/analyze.R") %>% 
  walk(source)

options(
  repos = "https://packagemanager.rstudio.com/cran/__linux__/focal/2022-01-28",
  Ncpus = 9
)

prevalence_group_colors <- c(
  "Bacteria Generalist" = "#3a86ff",
  "Bacteria Specialist" = "#ffbe0b",
  "Fungi Generalist" = "#3a86ff",
  "Fungi Specialist" = "#ffbe0b"
)

prevalence_type <- c(
  "Generalist" = "#3a86ff",
  "Specialist" = "#ffbe0b"
)

kingdoms_colors <- c(
  "Bacteria" = "#3f78a9",
  "Fungi" = "#A93F55",
  "all" = "black"
)



message("Global env bootstraped.")


