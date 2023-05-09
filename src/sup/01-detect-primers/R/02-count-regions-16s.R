library(tidyverse)

raw_coverage <- list.files(path = "/home/ailtonpcf/draco/proj/04-global-microbiome/cache/08-coverage", 
                           full.names = T) %>%
  as_tibble() %>% 
  mutate(coverage = map(value, ~ read_tsv(.x)))

cov <- raw_coverage %>%
  mutate(sample_id = str_split(value, "/", simplify=T)[,9]) %>%
  mutate(sample_id = str_remove(sample_id, ".txt")) %>% 
  select(-value) %>% 
  unnest(coverage)

runs <- read_csv(file = "/home/ailtonpcf/draco/proj/04-global-microbiome/raw/runs.csv")

proj16s <- cov %>%
  rename(rRNA = "#ID") %>% 
  left_join(runs %>% select(run_id, bioproject_id), by = c("sample_id" = "run_id")) %>% 
  relocate(sample_id, bioproject_id) 
# %>% 
#   mutate(rRNA = str_remove(rRNA, "_[:digit:]+")) %>% 
#   mutate(rRNA = str_remove(rRNA, "_amp")) 

proj16s %>% 
  group_by(rRNA) %>% 
  count() %>% 
  print(n=Inf)

proj16s %>% 
  filter(str_detect(bioproject_id, "3505")) %>% 
  filter(str_detect(sample_id, "8912")) %>% 
  arrange(desc(Covered_percent))

proj16s %>% 
  filter(str_detect(bioproject_id, "1113")) %>% 
  arrange(desc(Covered_percent))


