#' Plot the information transfer of a tree (drivers ordering).
#' 
#' @description 
#' 
#' The information transfer of the tree is the set of orderings
#' associated to the internal annotated driver events.
#' 
#' This function plots these orderings following a topological
#' sort of the node of the corresponding clone tree.
#'
#' @param x A \code{ctree} tree.
#' @param node_palette A function that can return, for an input number,
#' a number of colours. 
#' @param tree_layout Layout of this model, as of \code{ggraph}.
#' @param ... Other parameters, not used in this case.
#'
#' @return A \code{ggplot} object for the plot.
#' 
#' @import tidygraph
#' @import crayon
#' @import ggraph
#' 
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
#' plot_information_transfer(x[[1]])
#' 
#' # Change layout -- use igraph's "kk" layout
#' plot_information_transfer(x[[1]], tree_layout = 'kk')
plot_information_transfer = function(x,
                                     node_palette = colorRampPalette(RColorBrewer::brewer.pal(n = 9, "Set1")),
                                     tree_layout = 'tree',
                                     ...
)
{
  tree = x
  
  # Get the tidygraphs that we use here
  tb_tree = as_tbl_graph(tree$transfer$drivers) %>%
    activate(nodes) %>%
    rename(driver = name)
  tb_tree_all_tree = tree$tb_adj_mat
  
  # Color the nodes by cluster id, as in the plot of the tree
  # use a topological sort to pick the colors in the same order
  clones_orderings = igraph::topo_sort(
    igraph::graph_from_adjacency_matrix(DataFrameToMatrix(tree$transfer$clones)),
    mode = 'out'
  )$name
  
  nDrivers = length(clones_orderings) - 1 # avoid GL
  
  drivers_colors = c('white', node_palette(nDrivers))
  names(drivers_colors) = clones_orderings
  
  tb_tree = tb_tree %>%
    activate(nodes) %>%
    left_join(tree$drivers %>% rename(driver = variantID), by = 'driver') %>%
    mutate(cluster = ifelse(driver == "GL", 'GL', cluster))
  
  # Plot call
  layout <- create_layout(tb_tree, layout = tree_layout)
  
  ggraph(layout) +
    geom_edge_diagonal(
      arrow = arrow(length = unit(2.5, 'mm')),
      end_cap = circle(2.5, 'mm'),
      start_cap  = circle(2.5, 'mm'),
      color = 'steelblue'
    ) +
    geom_node_point(
      aes(colour = cluster, fill = cluster),
      alpha = .5,
      size = 6
    ) +
    geom_node_text(aes(label = driver), color = 'black', size = 3) +
    coord_cartesian(clip = 'off') +
    theme_void(base_size = 8) +
    theme(
      legend.position = 'bottom',
      legend.key.size = unit(3, "mm")
    ) +
    scale_fill_manual(values = drivers_colors) +
    scale_color_manual(
      values = drivers_colors,
      guide = guide_legend(override.aes = list(shape = 21, size = 3), alpha = 1)
    ) +
    labs(
      title = paste(tree$patient),
      subtitle = paste0('Information transfer')
    ) +
    guides(
      colour = guide_legend('Clone'),
      fill = FALSE
    )
}
