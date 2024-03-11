library(tidyverse)
library(RColorBrewer)

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
  mutate(rRNA = str_remove(rRNA, "_amp"))

#Plot1
p1 <- proj16s %>%
  ggplot(aes(x=Covered_percent)) +
  geom_boxplot(aes(color=rRNA)) +
  facet_wrap(c("bioproject_id"), scales = "free_y") +
  scale_color_manual(values = col_vector) +
  guides(fill=guide_legend(ncol=1)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 10),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  labs(x = "Coverage (%)",
       y = "Count",
       color = "rRNA")

ggsave(filename = "/home/ailtonpcf/draco/proj/04-global-microbiome/plots/04-manuscript-ap/01-sup-16s_boxplot.png",
       plot = p1,
       device = "png",
       width = 12,
       height = 9,
       dpi = 300)

ggsave(filename = "/home/ailtonpcf/draco/proj/04-global-microbiome/plots/04-manuscript-ap/01-sup-16s_boxplot.pdf",
       plot = p1,
       device = "pdf",
       width = 12,
       height = 9,
       dpi = 300)


p2 <- proj16s %>%
  ggplot(aes(x=Covered_percent)) +
  geom_density(aes(color=rRNA), alpha = 0.75) +
  facet_wrap(c("bioproject_id"), scales = "free_y") +
  scale_color_manual(values = col_vector) +
  guides(fill=guide_legend(ncol=1)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 10),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  labs(x = "Coverage (%)",
       y = "Count",
       color = "rRNA")

ggsave(filename = "/home/ailtonpcf/draco/proj/04-global-microbiome/plots/04-manuscript-ap/01-sup-16s_density.png",
       plot = p2,
       device = "png",
       width = 12,
       height = 9,
       dpi = 300)

ggsave(filename = "/home/ailtonpcf/draco/proj/04-global-microbiome/plots/04-manuscript-ap/01-sup-16s_density.pdf",
       plot = p2,
       device = "pdf",
       width = 12,
       height = 9,
       dpi = 300)

p3 <- proj16s %>%
  ggplot(aes(x=Covered_percent)) +
  geom_histogram(aes(fill=rRNA), alpha = 0.75) +
  facet_wrap(c("bioproject_id"), scales = "free_y") +
  scale_fill_manual(values = col_vector) +
  guides(fill=guide_legend(ncol=1)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 10),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  labs(x = "Coverage (%)",
       y = "Count",
       fill = "rRNA")

ggsave(filename = "/home/ailtonpcf/draco/proj/04-global-microbiome/plots/04-manuscript-ap/01-sup-16s_histogram.png",
       plot = p3,
       device = "png",
       width = 12,
       height = 9,
       dpi = 300)

ggsave(filename = "/home/ailtonpcf/draco/proj/04-global-microbiome/plots/04-manuscript-ap/01-sup-16s_histogram.pdf",
       plot = p3,
       device = "pdf",
       width = 12,
       height = 9,
       dpi = 300)
