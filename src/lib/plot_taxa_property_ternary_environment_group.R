#' Plot a taxa property in a ternary plot
#'
#' Axes of ternary plot: Enviornment group
#' Point color: Kingdom
#'
#' @param properties tibble with column data. e.g. foundations or keystones
#' @param property column name of `properties` holding values of taxa to plot
#' @param taxa_facetting column name of `data` of `properties` used to facett plots``
#' @param pooled_lineages tibble with column taxon and kingdom
plot_taxa_property_ternary_environment_group <- function(properties, property, taxa_facetting = "dominance", pooled_lineages) {
  properties %>%
    tidyr::unnest(data) %>%
    dplyr::filter(taxa_grouping == taxa_facetting & samples_grouping == "environment_group") %>%
    pivot_wider(names_from = samples_group, values_from = property) %>%
    left_join(pooled_lineages, by = "taxon") %>%
    ggplot(aes(x = soil, y = aquatic, z = host)) +
    geom_point(aes_string(color = "kingdom")) +
    facet_wrap(~taxa_group) +
    scale_color_kingdom() +
    ggtern::coord_tern() +
    ggtern::theme_custom(
      col.T = environment_group_colors$aquatic,
      col.R = environment_group_colors$host,
      col.L = environment_group_colors$soil,
      col.grid.minor = "grey"
    ) +
    ggplot2::theme(
      panel.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 18),
      plot.title = element_text(size = 20)
    ) +
    labs(
      title = taxa_facetting,
      color = "Kingdom"
    )
}
