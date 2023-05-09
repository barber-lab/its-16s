library(tidyverse)

# Primers from paper
################################################################################
proj1 = c("bioproject_id" = "PRJEB31590", 
          "primer_bac"    = "V3-V4", 
          "primer_euk"    = "ITS1")

proj2 = c("bioproject_id" = "PRJEB30970", 
          "primer_bac"    = "V3-V4", 
          "primer_euk"    = "ITS1")

proj3 = c("bioproject_id" = "PRJNA287840", 
          "primer_bac"    = "V3-V4", 
          "primer_euk"    = "18S-V1-V3")

proj4 = c("bioproject_id" = "PRJEB11419", 
          "primer_bac"    = "V4", 
          "primer_euk"    = "ITS1")

proj5 = c("bioproject_id" = "PRJEB32265", 
          "primer_bac"    = "V3-V4", 
          "primer_euk"    = "ITS1")

# primer_bac not in the paper > NA filled with own analysis
proj6 = c("bioproject_id" = "PRJNA478949", 
          "primer_bac"    = "V4", 
          "primer_euk"    = "ITS1")
# In the methods they amplify the region ITS1 but in the paper they said ITS2
# No region no primer used for bacteria

proj7 = c("bioproject_id" = "PRJNA241408_PRJNA522449", 
          "primer_bac"    = "V4", 
          "primer_euk"    = "ITS1")

proj8 = c("bioproject_id" = "PRJEB23282", 
          "primer_bac"    = "V4-V5", 
          "primer_euk"    = "ITS2")

proj9 = c("bioproject_id" = "PRJNA606949", 
          "primer_bac"    = "V3-V4", 
          "primer_euk"    = "ITS1")
# https://doi.org/10.1016/j.soilbio.2020.107953
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3395698/

proj10 = c("bioproject_id" = "PRJNA450848", 
           "primer_bac"    = "V4", 
           "primer_euk"    = "ITS1")

# primer_bac not in the paper > NA filled with own analysis
proj11 = c("bioproject_id" = "PRJNA418896", 
           "primer_bac"    = "V3-V4", 
           "primer_euk"    = "ITS2")

proj12 = c("bioproject_id" = "PRJNA550037", 
           "primer_bac"    = "V3-V4", 
           "primer_euk"    = "ITS2")

# primer_bac not in the paper > NA filled with own analysis
proj13 = c("bioproject_id" = "PRJNA561410_PRJNA561568", 
           "primer_bac"    = "V3-V4", 
           "primer_euk"    = "ITS2")

proj14 = c("bioproject_id" = "PRJNA415280_PRJNA415285", 
           "primer_bac"    = "V4", 
           "primer_euk"    = "ITS1")

proj15 = c("bioproject_id" = "PRJNA517449", 
           "primer_bac"    = "V4", 
           "primer_euk"    = "ITS2")

proj16 = c("bioproject_id" = "PRJNA526458", 
           "primer_bac"    = "V3-V4", 
           "primer_euk"    = "ITS1")

proj17 = c("bioproject_id" = "PRJNA492720", 
           "primer_bac"    = "V4", 
           "primer_euk"    = "ITS2")

proj18 = c("bioproject_id" = "PRJNA528359", 
           "primer_bac"    = "V4", 
           "primer_euk"    = "ITS1")

# Primers from my analysis
################################################################################

proj19 = c("bioproject_id" = "PRJNA282687", 
           "primer_bac"    = "V1-V2", 
           "primer_euk"    = "ITS1")

proj20 = c("bioproject_id" = "PRJNA324410", 
           "primer_bac"    = "V4", 
           "primer_euk"    = "ITS1")

proj21 = c("bioproject_id" = "PRJNA496065", 
           "primer_bac"    = "V1-V2", 
           "primer_euk"    = "ITS1")

proj22 = c("bioproject_id" = "PRJNA271113", 
           "primer_bac"    = "V1-V2", 
           "primer_euk"    = "ITS1")

proj23 = c("bioproject_id" = "PRJNA359237", 
           "primer_bac"    = "V1-V2", 
           "primer_euk"    = "ITS1-ITS2")

proj24 = c("bioproject_id" = "PRJNA473079", 
           "primer_bac"    = "V4", 
           "primer_euk"    = "ITS2")

proj25 = c("bioproject_id" = "PRJNA263505", 
           "primer_bac"    = "V3-V4", 
           "primer_euk"    = "ITS1")

proj26 = c("bioproject_id" = "PRJNA406830", 
           "primer_bac"    = "V4-V5", 
           "primer_euk"    = "ITS2")

proj27 = c("bioproject_id" = "PRJNA432446", 
           "primer_bac"    = "V3-V4", 
           "primer_euk"    = "ITS1")

primers <- paste("proj", seq(1:27), sep = "") %>% 
  str_split(" ", simplify = TRUE) %>% 
  map_dfr(~eval(parse(text = .)))

primers %>% 
  mutate(set1bac = if_else(str_detect(primer_bac, "V5"), "NA", "V1-V4")) %>% 
  mutate(set2bac = if_else(str_detect(primer_bac, "V4|V5"), "V3-V5", "NA")) %>%
  mutate(set1fun = if_else(str_detect(primer_euk, "ITS1"), "ITS1", "NA")) %>% 
  mutate(set2fun = if_else(str_detect(primer_euk, "ITS2"), "ITS2", "NA")) %>% 
  write_csv(file = "tbl/primers_comparison_16s_its.csv")
