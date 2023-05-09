#!/usr/bin/env R

# Defaults file to bootstrap a standard global environment

# Set the correct targets project
work_dir <- "04-quantify-primers-effect"
Sys.setenv(TAR_PROJECT = "primers_effect")

# parallelization of targets package
options(
    tidyverse.quiet = TRUE,
    clustermq.scheduler = "multiprocess"
)

library(targets)
library(tarchetypes)
library(here)
library(reshape)
library(patchwork)
#library(MatrixCorrelation)
library(broom)
library(tidyverse)
library(ggcorrplot2)

# if (!requireNamespace("devtools")) install.packages("devtools")
# devtools::install_github("caijun/ggcorrplot2")

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
  discard(~ .x == paste("src", work_dir, "analyze.R", sep = "/")) %>% 
  walk(source)

options(
  repos = "https://packagemanager.rstudio.com/cran/__linux__/focal/2022-01-28",
  Ncpus = 9
)

message("Global env bootstraped.")
