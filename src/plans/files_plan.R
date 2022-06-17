#
# Write targets to files on network storage to work with them from otehr computers
#

get_files_plan <- function() {
  drake_plan(
    tmp_object_dir = "/sbidata/dloos/prj/AX3-its-16s-stat/tmp/objects/",
    results_dir = "/sbidata/dloos/prj/AX3-its-16s-stat/results/",
    sub_abundances_file = write_rds(sub_abundances, paste0(tmp_object_dir, "sub_abundances.rds")),
    lineages_file = write_rds(lineages, paste0(tmp_object_dir, "lineages.rds")),
    generalists_and_specialists_file = write_rds(generalists_and_specialists, paste0(tmp_object_dir, "generalists_and_specialists.rds")),
    selected_generalists_specialists_file = write_rds(selected_generalists_specialists, paste0(tmp_object_dir, "selected_generalists_specialists.rds")),
    bioproject_prevalences_file = write_rds(bioproject_prevalences, paste0(tmp_object_dir, "bioproject_prevalences.rds")),
    bioproject_environment_groups_file = write_rds(bioproject_environment_groups, paste0(tmp_object_dir, "bioproject_environment_groups.rds")),
    fig1_overview_file = ggsave(paste0(results_dir, "fig1.png"), fig1_overview, width = 6, height = 9),
    stepwise_constrained_environment_group_plot_file = ggsave(paste0(results_dir, "fig4_biome.png"), stepwise_constrained_environment_group_plot, width = 6, height = 6),
    stepwise_constrained_plot_file = ggsave(paste0(results_dir, "fig4_pooled.png"), stepwise_constrained_plot, width = 6, height = 6),
    stepwise_constrained_ordinations_contribution_table_file = write_xlsx(stepwise_constrained_ordinations_contribution_table, paste0(results_dir, "stepwise_constrained_ordinations_contribution.xlsx")),
    selected_generalists_specialists_table_file = write_xlsx(selected_generalists_specialists_table, paste0(results_dir, "generalists_specialists.xlsx")),
    levins_plot_file = ggsave(paste0(results_dir, "levins.png"), levins_plot, width = 6, height = 3, dpi = 600, bg = "white"),
    fig2c_file = ggsave(paste0(results_dir, "fig2c.png"), generalists_and_specialists_abundances_plt$Bacteria, width = 3, height = 2, dpi = 600, bg = "white"),
    fig2d_file = ggsave(paste0(results_dir, "fig2d.png"), generalists_and_specialists_abundances_plt$Fungi, width = 3, height = 2, dpi = 600, bg = "white"),
    fig3_alphadiv_file = ggsave(paste0(results_dir, "fig3.png"), fig3_alphadiv, width = 8, height = 8, dpi = 600, bg = "white"),
    fig6_heatmap_file = target({
      png(paste0(results_dir, "fig6.png"), width = 10, height = 10, units = "in", res = 300)
      draw(fig6_heatmap, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = TRUE)
      dev.off()
    }),
    correlation_common_taxa_across_environment_group_topology_plt_file = target(
      {
        ggsave(paste0(results_dir, "correlation_common_taxa_across_environment_group_topology_plt.png"),
          correlation_common_taxa_across_environment_group_topology_plt,
          width = 6, height = 4
        )
      },
      hpc = FALSE
    ),
    generalists_common_coabundance_graph_plt_file = target(
      {
        ggsave(paste0(results_dir, "generalists_common_coabundance_graph_plt.png"),
          generalists_common_coabundance_graph_plt,
          width = 6, height = 4
        )
      },
      hpc = FALSE
    ),

    all_files = target({
      dependencies <- c(
        sub_abundances_file,
        lineages_file,
        generalists_and_specialists_file,
        selected_generalists_specialists_file,
        bioproject_prevalences_file,
        bioproject_environment_groups_file,
        fig1_overview_file,
        fig2c_file,
        fig2d_file,
        fig3_alphadiv_file,
        levins_plot_file,
        fig6_heatmap_file,
        stepwise_constrained_environment_group_plot_file,
        stepwise_constrained_plot_file,
        stepwise_constrained_ordinations_contribution_table_file,
        selected_generalists_specialists_table_file,
        correlation_common_taxa_across_environment_group_topology_plt_file,
        generalists_common_coabundance_graph_plt_file
      )

      # dummy output to not actually store the outputs
      "done"
    })
  )
}
