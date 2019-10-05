#' Plot CCF clusters data (tile).
#' 
#' @description This function creates a \code{ggplot}-style
#' heatmap of the input CCF cluster of each clone in the data.
#' The heatmap is annotated for the drivers status of each
#' clone (with/ without driver). The CCF values are used to
#' colour the heatmap (`NA` values are in light gray).
#'
#' @param x A \code{ctree} tree.
#' @param patient A patient id.
#' @param ... Extra parameters, not used.
#'
#' @return A \code{ggplot} plot.
#' @export
#'
#' @examples
#' data(ctree_input)
#' 
#' x = ctrees(
#' ctree_input$CCF_clusters, 
#' ctree_input$drivers,
#' ctree_input$samples,
#' ctree_input$patient,
#' ctree_input$sspace.cutoff,
#' ctree_input$n.sampling,
#' ctree_input$store.max
#' )
#' 
#' plot_CCF_clusters(x[[1]])
plot_CCF_clusters = function(x, ...)
{
  p = x$patient
  sm = x$samples
  cl = x$CCF

  # Values CCF
  cl_tab = cl %>%
    select(!!sm, cluster) %>%
    reshape2::melt(id = 'cluster') %>%
    rename(
      region = variable,
      CCF = value) %>%
    as_tibble() %>%
    mutate(CCF = ifelse(CCF == 0, NA, CCF))
  
  # Annotations
  cl_tab_anno = cl %>%
     select(cluster, nMuts, is.driver, is.clonal)
  
  # Cluster ordering by sum of CCF
  cluster_ordering = cl_tab %>%
    group_by(cluster) %>%
    summarise(tot = sum(CCF, na.rm = TRUE)) %>%
    arrange(tot) %>%
    pull(cluster)
  
  # Combined
  cl_tab = cl_tab %>%
    left_join(cl_tab_anno %>% select(cluster, is.driver), by = 'cluster') %>%
    mutate(is.driver = ifelse(is.na(CCF), FALSE, is.driver))
  
  # Factor to sort
  cl_tab$cluster = factor(cl_tab$cluster, levels = cluster_ordering)
  
  ggplot(
    cl_tab,
    aes(
      x = region,
      y = cluster,
      z = CCF,
      fill = CCF,
      color = is.driver)
  ) +
    geom_tile(aes(width = .8, height = .8), size = 1) +
    geom_text(aes(label = CCF), color = 'orange') +
    scale_fill_distiller(palette = 'Blues', na.value = 'gainsboro', direction = 1) +
    scale_color_manual(
      values = c(`TRUE` = 'red', `FALSE` = NA)
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = 'bottom',
      legend.key.size = unit(3, 'mm')
    ) +
    labs(
      title = paste("Clones for", p),
      y = 'Cluster',
      x = 'Region',
      subtitle = paste("Clonal cluster = ", cl$cluster[cl$is.clonal])
    ) +
    guides(
      fill = guide_colourbar(barwidth = 6),
      color = guide_legend("With driver")
      )
}

