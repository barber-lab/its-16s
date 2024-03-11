# Most ecological papers use Levins' niche breadth index to define generalist taxon
source("src/defaults.R")

library(MicroNiche)

loadd(sub_abundances)
loadd(selected_generalists_specialists)
loadd(samples)

taxa_data <-
  sub_abundances %>%
  filter(norm_method == "tss" & subset_value == "all") %>%
  unnest(data) %>%
  # pseudo counts to normalize for seqdepth
  mutate(abundance = round(abundance * 10e3)) %>%
  select(sample_id, taxon, abundance) %>%
  pivot_wider(names_from = sample_id, values_from = abundance, values_fill = list(abundance = 0)) %>%
  as.data.frame()

levins_indicies <-
  tibble(sample_grouping = c("bioproject_id", "environment_group", "habitat")) %>%
  mutate(
    levins_index = sample_grouping %>% map(possibly(~ {
      sample_data <-
        taxa_data %>%
        colnames() %>%
        tibble(sample_id = .) %>%
        left_join(samples) %>%
        pluck(.x)

      n_sample_groups <- sample_data %>%
        unique() %>%
        length()

      levins.Bn(taxa_data, n_sample_groups, sample_data) %>%
        as_tibble(rownames = "taxon")
    }, NA))
  )

#saveRDS(object = levins_indicies, file = "rdata/levins_indicies.rds")
levins_indicies <- readRDS("rdata/levins_indicies.rds")

prevalence_group_colors_upper <- c(
  "Generalist" = "#3a86ff",
  "Specialist" = "#ffbe0b"
)

scale_color_prevalence_group_upper <- function(drop = FALSE, ...) {
  scale_color_manual(
    values =  purrr::simplify(prevalence_group_colors_upper),
    na.value = "#000000",
    drop = drop,
    ...
  )
}

p1 <- levins_indicies %>%
  unnest(levins_index) %>%
  left_join(selected_generalists_specialists) %>%
  filter(!is.na(prevalence_group)) %>%
  mutate(
    sample_grouping = sample_grouping %>%
      recode("bioproject_id" = "Project", "habitat" = "Habitat", "environment_group" = "Environment") %>%
      factor(levels = c("Project", "Habitat", "Environment"))
  ) %>%
  mutate(prevalence_group = recode(prevalence_group, "generalist" = "Generalist", "specialist" = "Specialist")) %>% 
  ggplot(aes(prevalence_group, Bn, color = prevalence_group)) +
  geom_boxplot() +
  scale_color_prevalence_group_upper() +
  stat_compare_means(
    method = "wilcox",
    comparisons = list(c("Generalist", "Specialist"))
  ) +
  facet_wrap(~sample_grouping) +
  guides(color = "none") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        strip.text = element_text(size = rel(1.2)),
        legend.text = element_text(size=rel(1.2)),
        axis.text = element_text(size=rel(1.2)),
        axis.title = element_text(size=rel(1.2)),
        panel.grid.minor.x = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  labs(
    x = "",
    y = TeX("Levins' niche breadth index B_n"),
    tag = "A"
  )

# SNB
################################################################################
library(ggbeeswarm)

generalists <- read_tsv(here::here("raw/snb/generalist_edit.txt"), skip = 1,
                        col_names = "taxa")  %>% 
  mutate(gen_spec = "generalist")

specialists <- read_tsv(here::here("raw/snb/specialists_edit.txt"), skip = 1,
                        col_names = "taxa") %>% 
  mutate(gen_spec = "specialist")

gen_spec <- bind_rows(generalists, specialists) %>% 
  mutate(taxa = str_replace_all(taxa, " ", "_"))

snb <- read_tsv(here::here("raw/snb/bas_snb.txt")) %>% 
  filter(rank == "genus") %>% 
  mutate(taxa = str_extract(`taxonomic lineage`, "genus\\.(.+)"),
         taxa = str_remove(taxa, "genus\\.")
  ) %>% 
  relocate(taxa, .after = rank)

plot_df <- gen_spec %>% 
  left_join(snb)

p2 <- snb %>% 
  ggplot(aes(x = "1", y = SNB)) +
  geom_violin() +
  geom_quasirandom(data = plot_df, aes(x = "1", y = SNB, color = gen_spec), 
                   width = 0.07, size = 3) +
  scale_color_manual(values = c("#3a86ff", "#ffbe0b"), name = "", 
                     labels = c("Generalist", "Specialist")) +
  labs(x = "", tag = "B") +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        legend.position=c(.7, .9),
        legend.key = element_blank(),
        legend.background=element_blank(),
        strip.text = element_text(size = rel(1.2)),
        legend.text = element_text(size=rel(1.2)),
        axis.text = element_text(size=rel(1.2)),
        axis.title = element_text(size=rel(1.2)),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

p <- p1 + p2 + plot_layout(widths = c(2,1))

ggsave(filename = "results/sup-fig3.png", 
       plot = p, 
       device = "png",
       width = 11, 
       height = 7, 
       dpi = 300)

ggsave(filename = "results/sup-fig3.pdf", 
       plot = p, 
       device = "pdf",
       width = 11, 
       height = 7, 
       dpi = 300)
