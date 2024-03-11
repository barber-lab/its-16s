#!/usr/bin/env R

# Defaults file to bootstrap a standard global environment

# Set the correct targets project
Sys.setenv(TAR_PROJECT = "aln-metrics")

# parallelization of targets package
options(
    tidyverse.quiet = TRUE,
    clustermq.scheduler = "multiprocess"
)

library(targets)
library(clustermq)
library(tarchetypes)
library(tidyverse)
library(here)

setwd(here::here())

# do not recursively call itself
# discards additional source files
# discard the targets file itself
list.files("src/0r-microbiome", pattern = ".R$", full.names = TRUE, recursive = F) %>%
  discard(~ .x == "src/aln-metrics/defaults.R") %>%
  discard(~ .x == "src/aln-metrics/_targets.R") %>%
  discard(~ .x == "src/aln-metrics/analyze.R") %>%
  discard(~ .x == "src/aln-metrics/utils.R") %>%
  discard(~ .x == "src/aln-metrics/install.R") %>%
  walk(source)

message("Global env bootstraped.")


