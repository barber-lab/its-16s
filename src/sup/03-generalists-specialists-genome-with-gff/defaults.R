#!/usr/bin/env R
# Defaults file to bootstrap a standard global environment

# Set the correct targets project
work_dir <- "03-generalists-specialists-genome-with-gff"
Sys.setenv(TAR_PROJECT = "generalist_specialist_genome_gff")

# parallelization of targets package
# options(
#     tidyverse.quiet = TRUE,
#     clustermq.scheduler = "multiprocess"
# )

library(targets)
library(tarchetypes)
library(tidyverse)
library(here)
library(jsonlite)
library(ggpubr)
library(ggsignif)
library(patchwork)
library(future)
library(future.callr)
library(RColorBrewer)

setwd(here::here())

# do not recursively call itself
# discards additional source files
# discard the targets file itself

list.files(path = paste("src/", work_dir, sep = ""), 
           pattern = ".R$", 
           full.names = TRUE, 
           recursive = F) %>%
  discard(~ .x == paste("src", work_dir, "defaults.R", sep = "/")) %>%
  discard(~ .x == paste("src", work_dir, "_targets.R", sep = "/")) %>%
  discard(~ .x == paste("src", work_dir, "analyze_local.R", sep = "/")) %>%
  discard(~ .x == paste("src", work_dir, "analyze_server.R", sep = "/")) %>%
  discard(~ .x == paste("src", work_dir, "install.R", sep = "/")) %>% 
  walk(source)

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

message("Global env bootstraped")


