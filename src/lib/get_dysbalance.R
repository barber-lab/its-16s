get_dysbalance <- function(distance, max_balanced_quantile = 0.5, max_clusters = 5) {
  # prevent more clusters than samples
  max_clusters <- min(
    max_clusters,
    distance %>% as.matrix() %>% dim() %>% pluck(1)
  )

  clustering_search <-
    distance %>%
    fpc::pamk(diss = TRUE, krange = 1:max_clusters)

  clustering <- clustering_search$pamobject
  clustering_tbl <-
    clustering$clustering %>%
    as_tibble(rownames = "sample_id") %>%
    rename(cluster = value)

  dists_to_medoid_tbl <-
    clustering_tbl %>%
    inner_join(
      clustering$medoids %>% tibble(cluster = length(.) %>% seq(), medoid = .),
      by = "cluster"
    ) %>%
    inner_join(
      distance %>%
        as.matrix() %>%
        as_tibble(rownames = "item1") %>%
        tidyr::gather(item2, distance, -item1) %>%
        rename(sample_id = item1, medoid = item2),
      by = c("sample_id", "medoid")
    ) %>%
    dplyr::transmute(
      sample_id,
      distance_to_medoid = ifelse(sample_id == medoid, 0, distance),
      cluster = cluster %>% as.factor()
    )

  # calculate dysbalances
  dysbalance_tbl <-
    dists_to_medoid_tbl %>%
    group_by(cluster) %>%
    mutate(
      rank = rank(distance_to_medoid),
      norm_rank = (rank - 1) / n(),
      dysbalance = norm_rank %>% {
        ifelse(. <= max_balanced_quantile, "balanced", "dysbalanced")
      }
    ) %>%
    ungroup()

  dysbalance_thresholds_tbl <-
    dysbalance_tbl %>%
    group_by(cluster, dysbalance) %>%
    filter(dysbalance != "medium") %>%
    dplyr::summarise(
      min_dist_to_medoid = min(distance_to_medoid),
      max_dist_to_medoid = max(distance_to_medoid)
    )

  plt <-
    dists_to_medoid_tbl %>%
    inner_join(dysbalance_tbl, by = c("sample_id", "distance_to_medoid", "cluster")) %>%
    mutate_at("cluster", as.numeric) %>%
    ggplot(aes(x = cluster)) +
    geom_rect(
      data = dysbalance_thresholds_tbl %>% ungroup() %>% mutate_at("cluster", as.numeric),
      mapping = aes(
        xmin = cluster - 0.5,
        xmax = cluster + 0.5,
        ymin = min_dist_to_medoid,
        ymax = max_dist_to_medoid,
        fill = dysbalance
      )
    ) +
    geom_violin(aes(y = distance_to_medoid, group = cluster)) +
    geom_boxplot(aes(y = distance_to_medoid, group = cluster), width = 0.5) +
    scale_x_continuous(breaks = seq(clustering_search$nc)) +
    scale_fill_dysbalance() +
    labs(
      x = "Cluster",
      y = "distance to medoid (bray)"
    )

  list(
    cluster_search = clustering_search,
    dysbalance_thresholds_tbl = dysbalance_thresholds_tbl,
    dysbalance_tbl = dysbalance_tbl,
    dists_to_medoid_tbl = dists_to_medoid_tbl,
    plt = plt
  )
}
