source("src/04-quantify-primers-effect/defaults.R")
library(ggpubr)

# Import results generated before
tar_load(fungal_primers_heat_rds)
tar_load(bacterial_primers_heat_rds)

venn <- read_rds("cache/15-microbiome-rserver-hki-rdata/venn_no_aqua_fixed.rds")
prevalent_taxa <- readRDS(file = "cache/15-microbiome-rserver-hki-rdata/02supfig-prevalent_taxa.rds")
sub_abundances <- readRDS(file = "cache/15-microbiome-rserver-hki-rdata/02supfig-sub_abundances.rds")
primers_by_sample <- readRDS(file = "cache/15-microbiome-rserver-hki-rdata/02supfig-primers_by_sample.rds")

environment_group_colors <- c(
  "soil" = "#8b786d",
  "aquatic" = "#78a1bb",
  "host" = "#e0aa5a",
  "all" = "black"
)


scale_color_environment_group <- function(drop = FALSE, ...) {
  scale_color_manual(
    values = environment_group_colors,
    na.value = "#d9d9d9",
    drop = drop,
    ...
  )
}

symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                    symbols = c("****", "***", "**", "*", "ns"))

bacv1_v4 <-
  prevalent_taxa %>%
  nest(-kingdom) %>%
  filter(kingdom == "Bacteria") %>% 
  mutate(
    data = data %>% map(~ {
      .x %>%
        mutate(
          group = case_when(
            soil & host & aquatic ~ "Common\nin all",
            soil & !host & !aquatic ~ "Unique to\nsoil",
            !soil & !host & aquatic ~ "Unique to\naquatic",
            !soil & host & !aquatic ~ "Unique to\nhost"
          ),
          environment_group = group %>% recode(
            "Common\nin all" = "all", "Unique to\nsoil" = "soil", "Unique to\naquatic" = "aquatic",
            "Unique to\nhost" = "host"
          )
        ) %>%
        left_join(
          sub_abundances$data[[2]] %>% select(sample_id, taxon, abundance)
        ) %>%
        filter(abundance > 0.01e-2) %>%
        filter(!is.na(group)) %>%
        left_join(sub_abundances$data[[2]]) %>%
        left_join(primers_by_sample %>% select(sample_id, contains("bac"))) %>%
        filter(!is.na(set1bac)) %>%
        ggplot(aes(group, abundance, color = environment_group)) +
        geom_boxplot(width = 0.5) +
        stat_compare_means(
          method = "wilcox",
          symnum.args = symnum.args,
          comparisons = list(
            c("Common\nin all", "Unique to\naquatic"),
            c("Common\nin all", "Unique to\nhost"),
            c("Common\nin all", "Unique to\nsoil")
          )
        ) +
        scale_y_log10(expand = c(0, 0.5)) +
        scale_color_environment_group() +
        annotation_logticks(sides = "l") +
        guides(color = FALSE) +
        theme_bw() +
        theme(
          strip.text = element_text(size = rel(1.1)),
          axis.text = element_text(size = rel(1.1)), 
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_blank(),
          axis.ticks.x=element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) +
        labs(x = "", y = "Abundance (TSS)") +
        facet_wrap(c("set1bac"))
        
    })
  ) %>%
  deframe() %>% 
  pluck(1)

bacv3_v5 <-
  prevalent_taxa %>%
  nest(-kingdom) %>%
  filter(kingdom == "Bacteria") %>% 
  mutate(
    data = data %>% map(~ {
      .x %>%
        mutate(
          group = case_when(
            soil & host & aquatic ~ "Common\nin all",
            soil & !host & !aquatic ~ "Unique to\nsoil",
            !soil & !host & aquatic ~ "Unique to\naquatic",
            !soil & host & !aquatic ~ "Unique to\nhost"
          ),
          environment_group = group %>% recode(
            "Common\nin all" = "all", "Unique to\nsoil" = "soil", "Unique to\naquatic" = "aquatic",
            "Unique to\nhost" = "host"
          )
        ) %>%
        left_join(
          sub_abundances$data[[2]] %>% select(sample_id, taxon, abundance)
        ) %>%
        filter(abundance > 0.01e-2) %>%
        filter(!is.na(group)) %>%
        left_join(sub_abundances$data[[2]]) %>%
        left_join(primers_by_sample %>% select(sample_id, contains("bac"))) %>%
        filter(!is.na(set2bac)) %>% 
        ggplot(aes(group, abundance, color = environment_group)) +
        geom_boxplot(width = 0.5) +
        stat_compare_means(
          method = "wilcox",
          symnum.args = symnum.args,
          comparisons = list(
            c("Common\nin all", "Unique to\naquatic"),
            c("Common\nin all", "Unique to\nhost"),
            c("Common\nin all", "Unique to\nsoil")
          )
        ) +
        scale_y_log10(expand = c(0, 0.5)) +
        scale_color_environment_group() +
        annotation_logticks(sides = "l") +
        guides(color = FALSE) +
        theme_bw() +
        theme(
          strip.text = element_text(size = rel(1.1)),
          axis.text = element_text(size = rel(1.1)), 
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_blank(),
          axis.ticks.x=element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) +
        labs(x = "", y = "Abundance (TSS)") +
        facet_wrap(c("set2bac"))
      
    })
  ) %>%
  deframe() %>% 
  pluck(1)

funits1 <-
  prevalent_taxa %>%
  nest(-kingdom) %>%
  filter(kingdom == "Fungi") %>% 
  mutate(
    data = data %>% map(~ {
      .x %>%
        mutate(
          group = case_when(
            soil & host & aquatic ~ "Common\nin all",
            soil & !host & !aquatic ~ "Unique to\nsoil",
            !soil & !host & aquatic ~ "Unique to\naquatic",
            !soil & host & !aquatic ~ "Unique to\nhost"
          ),
          environment_group = group %>% recode(
            "Common\nin all" = "all", "Unique to\nsoil" = "soil", "Unique to\naquatic" = "aquatic",
            "Unique to\nhost" = "host"
          )
        ) %>%
        left_join(
          sub_abundances$data[[2]] %>% 
            drop_na(environment_group) %>% 
            select(sample_id, taxon, abundance)
        ) %>%
        drop_na(environment_group) %>%
        filter(abundance > 0.01e-2) %>%
        filter(!is.na(group)) %>%
        left_join(sub_abundances$data[[2]]) %>%
        left_join(primers_by_sample %>% select(sample_id, contains("fun"))) %>%
        filter(!is.na(set1fun)) %>% 
        ggplot(aes(group, abundance, color = environment_group)) +
        geom_boxplot(width = 0.5) +
        stat_compare_means(
          method = "wilcox",
          symnum.args = symnum.args,
          comparisons = list(
            c("Common\nin all", "Unique to\naquatic"),
            c("Common\nin all", "Unique to\nhost"),
            c("Common\nin all", "Unique to\nsoil")
          )
        ) +
        scale_y_log10(expand = c(0, 0.5)) +
        scale_color_environment_group() +
        annotation_logticks(sides = "l") +
        guides(color = FALSE) +
        theme_bw() +
        theme(
          strip.text = element_text(size = rel(1.1)),
          axis.text = element_text(size = rel(1.1)), 
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_blank(),
          axis.ticks.x=element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) +
        labs(x = "", y = "Abundance (TSS)") +
        facet_wrap(c("set1fun"))
      
    })
  ) %>%
  deframe() %>% 
  pluck(1)

funits2 <-
  prevalent_taxa %>%
  nest(-kingdom) %>%
  filter(kingdom == "Fungi") %>% 
  mutate(
    data = data %>% map(~ {
      .x %>%
        mutate(
          group = case_when(
            soil & host & aquatic ~ "Common\nin all",
            soil & !host & !aquatic ~ "Unique to\nsoil",
            !soil & !host & aquatic ~ "Unique to\naquatic",
            !soil & host & !aquatic ~ "Unique to\nhost"
          ),
          environment_group = group %>% recode(
            "Common\nin all" = "all", "Unique to\nsoil" = "soil", "Unique to\naquatic" = "aquatic",
            "Unique to\nhost" = "host"
          )
        ) %>%
        left_join(
          sub_abundances$data[[2]] %>% 
            drop_na(environment_group) %>% 
            select(sample_id, taxon, abundance)
        ) %>%
        filter(abundance > 0.01e-2) %>%
        filter(!is.na(group)) %>%
        left_join(sub_abundances$data[[2]]) %>%
        left_join(primers_by_sample %>% select(sample_id, contains("fun"), env_check)) %>%
        filter(!is.na(set2fun)) %>%
        #filter(!(is.na(env_check) & environment_group == "aquatic")) %>%
        ggplot(aes(group, abundance, color = environment_group)) +
        geom_boxplot(width = 0.5) +
        stat_compare_means(
          method = "wilcox",
          symnum.args = symnum.args,
          comparisons = list(
            c("Common\nin all", "Unique to\nhost"),
            c("Common\nin all", "Unique to\nsoil")
          )
        ) +
        scale_y_log10(expand = c(0, 0.5)) +
        scale_color_environment_group() +
        annotation_logticks(sides = "l") +
        guides(color = FALSE) +
        theme_bw() +
        theme(
          strip.text = element_text(size = rel(1.1)),
          axis.text = element_text(size = rel(1.1)), 
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_blank(),
          axis.ticks.x=element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) +
        labs(x = "", y = "Abundance (TSS)") +
        facet_wrap(c("set2fun"))
      
    })
  ) %>%
  deframe() %>% 
  pluck(1)

p <- wrap_plots(
  venn[[1]],
  venn[[2]],
  venn[[3]],
  venn[[4]],
  (bacterial_primers_heat_rds + labs(tag = "B")),
  (fungal_primers_heat_rds + labs(tag = "C")),
  bacv1_v4 + labs(tag = "D"),
  bacv3_v5,
  funits1,
  funits2,
  ncol = 2,
  widths = c(1,1)
)

ggsave(filename = "plots/04-manuscript-ap/sup-fig2.pdf",
       plot = p,
       device = "pdf",
       width = 33,
       height = 40,
       units = "cm",
       dpi = 300)

ggsave(filename = "plots/04-manuscript-ap/sup-fig2.png",
       plot = p,
       device = "png",
       width = 33,
       height = 40,
       units = "cm",
       dpi = 300)
