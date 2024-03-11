library(tidyverse)
source("src/08-decontamination/lib/07-decontaminated-its.R")
source("src/08-decontamination/lib/08-decontaminated-16s.R")

# Load scripts 16s and its in the same environment first, then create the figure

plt <- proj_16s_general %>% 
  bind_rows(proj_16s_general_70) %>% 
  mutate(amplicon = "16S") %>% 
  mutate(BioProject = str_glue("{BioProject}\n{habitat}")) %>% 
  bind_rows(
    proj_its_general %>% 
      bind_rows(proj_its_general_70) %>% 
      mutate(amplicon = "ITS") %>% 
      mutate(BioProject = str_glue("{BioProject}\n{habitat}"))
  ) %>% 
  mutate(plt_proj = pmap(list(data, amplicon, BioProject, boxp_wid, dot_size), ~ plt_decontaminated(.x, .y, ..3, ..4, ..5)))

layout <- '
AAAA
BCCG
DEFF
'

paper_plt_merged <- 
  wrap_plots(
    guides = "collect",
    plt$plt_proj[[1]] + facet_wrap(c("ctrl_type"), nrow = 1) + labs(tag = "A"),
    plt$plt_proj[[2]] + labs(tag = "B"),
    plt$plt_proj[[3]] + labs(tag = "C"),
    plt$plt_proj[[4]] + labs(tag = "D"),
    plt$plt_proj[[6]] + labs(tag = "E"),
    plt$plt_proj[[5]] + labs(tag = "F"),
    guide_area(),
    design = layout
  )

ggsave(filename = paste0(plts_dir, "/", "01-paper-generalists-16s-its-merged.pdf"),
       plot = paper_plt_merged,
       width = 12,
       height = 14,
       dpi=150)

ggsave(filename = paste0(plts_dir, "/", "01-paper-generalists-16s-its-merged.png"),
       plot = paper_plt_merged,
       width = 12,
       height = 14,
       dpi=150)
