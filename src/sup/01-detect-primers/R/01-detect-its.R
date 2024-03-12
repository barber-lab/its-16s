library(tidyverse)

# ITS sequences after ITSx extraction
its <- read_delim(file = "cache/00-tmp/extracted_its_seqs_NAs_primers.txt", 
                col_names = F)

# Parse fasta headers from grep
seqs <- its %>%
  select(X1) %>% 
  mutate(sample_id = str_split(X1, "\\.", simplify=T)[,1]) %>% 
  mutate(ITS = str_extract(X1, "ITS[1-2]"))

# Import bioprojects and run_ids
runs <- read_csv(file = "raw/runs.csv")

# Merge run_id and bioprojects
sampleByITS <- seqs %>% 
  drop_na(ITS) %>% 
  group_by(sample_id, ITS) %>%
  count() %>% 
  #pivot_wider(names_from = its, values_from = n, values_fill = 0) %>% 
  left_join(runs %>% select(run_id, bioproject_id), by = c("sample_id" = "run_id"))

# Total counts by project and primer
sampleByITS %>% 
  ungroup() %>% 
  select(-sample_id) %>% 
  group_by(bioproject_id, ITS) %>% 
  summarise(total = sum(n)) %>% 
  filter(!total < 100) %>% 
  arrange(total) %>% 
  ggplot(aes(fct_reorder(bioproject_id, total), total)) +
  geom_col(aes(fill = ITS), position = position_dodge(width=0.6), width = 0.5) +
  scale_y_continuous(trans = "log10") +
  ylab("Log10(Reads number)") +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank())
  
ggsave(filename = "plots/01-detect-primers/fungi_its_sum_bioproject.png",
       device = "png", 
       width = 12, 
       height = 9, 
       dpi = 150)
  

# Primers used
primers <- sampleByITS %>% 
  ggplot(aes(x=n)) +
  geom_histogram(aes(fill=ITS)) +
  facet_wrap(c("bioproject_id")) +
  xlab("Sequences number") +
  ylab("Count")

ggsave(filename = "plots/01-detect-primers/fungi_its_amplified.png",
       device = "png", 
       width = 12, 
       height = 9, 
       dpi = 150)

  









