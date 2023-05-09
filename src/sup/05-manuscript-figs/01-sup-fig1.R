library(tidyverse)
library(RColorBrewer)
library(patchwork)

# 16S

raw_coverage <- list.files(path = "cache/08-coverage",
                           full.names = T) %>%
  as_tibble() %>%
  #slice(1:3) %>% 
  mutate(coverage = map(value, ~ read_tsv(.x)))

cov <- raw_coverage %>%
  mutate(sample_id = str_split(value, "/", simplify=T)[,3]) %>%
  mutate(sample_id = str_remove(sample_id, ".txt")) %>%
  select(-value) %>%
  unnest(coverage)

runs <- read_csv(file = "/home/ailtonpcf/draco/proj/04-global-microbiome/raw/runs.csv")

# Color palete
n <- 10
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# Merge dfs
proj16s <- cov %>%
  rename(rRNA = "#ID") %>%
  left_join(runs %>% select(run_id, bioproject_id), by = c("sample_id" = "run_id")) %>%
  relocate(sample_id, bioproject_id) %>%
  mutate(rRNA = str_remove(rRNA, "_[:digit:]+")) %>%
  mutate(rRNA = str_remove(rRNA, "_amp")) %>% 
  mutate(bioproject_id = str_replace(bioproject_id, "_", "\n"))

#Plot1
p1 <- proj16s %>%
  ggplot(aes(x=Covered_percent)) +
  geom_boxplot(aes(color=rRNA)) +
  facet_wrap(c("bioproject_id"), scales = "free_y") +
  scale_color_manual(values = col_vector) +
  guides(fill=guide_legend(ncol=1)) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  labs(x = "Coverage (%)",
       y = "Count",
       color = "rRNA", 
       tag = "A")

p2 <- proj16s %>%
  ggplot(aes(x=Covered_percent)) +
  geom_density(aes(color=rRNA), alpha = 0.75) +
  facet_wrap(c("bioproject_id"), scales = "free_y") +
  scale_color_manual(values = col_vector) +
  guides(fill=guide_legend(ncol=1)) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  labs(x = "Coverage (%)",
       y = "Count",
       color = "rRNA", 
       tag = "A")

p3 <- proj16s %>%
  ggplot(aes(x=Covered_percent)) +
  geom_histogram(aes(fill=rRNA), alpha = 0.75) +
  facet_wrap(c("bioproject_id"), scales = "free_y") +
  scale_fill_manual(values = col_vector) +
  guides(fill=guide_legend(ncol=1)) +
  theme_bw() +
  theme(axis.text = element_text(size = rel(1.1)),
        strip.text = element_text(size = rel(1.1)),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  labs(x = "Coverage (%)",
       y = "Number of reads",
       fill = "rRNA", 
       tag = "A")

# ITS

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

# Color palete
n <- 2
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


# Total counts by project and primer
p4 <- sampleByITS %>% 
  ungroup() %>% 
  select(-sample_id) %>% 
  group_by(bioproject_id, ITS) %>% 
  summarise(total = sum(n)) %>% 
  filter(!total < 100) %>% 
  arrange(total) %>% 
  ggplot(aes(fct_reorder(bioproject_id, total), total)) +
  geom_col(aes(fill = ITS), position = "fill", width = 0.5) +
  coord_flip() +
  scale_fill_manual(values = col_vector) +
  guides(fill=guide_legend(ncol=1)) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        strip.text = element_text(size = rel(1.1)),
        axis.text = element_text(size = rel(1.1)), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        axis.ticks.x=element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  labs(y = "Reads (%)",
       fill = "ITS", 
       tag = "B")


wp3 <- p3 + p4 + plot_layout(widths = c(3,1), guides = "collect")

ggsave(filename = "/home/ailtonpcf/draco/proj/04-global-microbiome/plots/04-manuscript-ap/sup-fig1.png",
       plot = wp3,
       device = "png",
       width = 13,
       height = 9,
       dpi = 300)

ggsave(filename = "/home/ailtonpcf/draco/proj/04-global-microbiome/plots/04-manuscript-ap/sup-fig1.pdf",
       plot = wp3,
       device = "pdf",
       width = 13,
       height = 9,
       dpi = 300)




