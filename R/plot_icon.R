#' Plot the information transfe among clones in icon format (tiny).
#' 
#' @description 
#' 
#' The information transfer of a tree is the set of orderings
#' associated to the internal annotated driver events.
#' 
#' This function plots the clones with drivers (so their orderings) following a topological
#' sort of the node of the corresponding clone tree, and in tiny
#' icon format. Some graohics changes with respect ot a standard plot.
#'
#' @param x A \code{ctree} tree.
#' @param node_palette A function that can return, for an input number,
#' a number of colours. 
#' @param ... Other parameters, not used in this case.
#'
#' @return A \code{ggplot} object for the plot.
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
#' plot_icon(x[[1]])
plot_icon = function(x,
                     node_palette = colorRampPalette(RColorBrewer::brewer.pal(n = 9, "Set1")),
                     ...
)
{
    # Get the tidygraph
    tree = x
    tb_tree = tree$tb_adj_mat
    
    # Color the nodes by cluster id
    tb_node_colors = tb_tree %>% filter(is.driver) %>% pull(cluster)
    
    tb_node_colors = node_palette(length(tb_node_colors))
    tb_node_colors = c(tb_node_colors, `GL` = 'white')
    names(tb_node_colors) = c(tb_tree %>% filter(is.driver) %>% pull(cluster), 'GL')
    
    # Graph from transfer
    tb_icon = as_tbl_graph(tree$transfer$clones) %>%
      rename(cluster = name) %>%
      activate(edges) %>%
      mutate(cluster = .N()$cluster[from])
    
    # Plot call
    layout <- create_layout(tb_icon, layout = 'tree')
    
    ggraph(layout) +
      geom_edge_link(aes(colour = cluster),
                     width = 1) +
      geom_node_point(aes(colour = cluster),
                      size = 2.5) +
      coord_cartesian(clip = 'off') +
      theme_void(base_size = 8) +
      theme(legend.position = 'none') +
      scale_color_manual(values = tb_node_colors) +
      scale_edge_color_manual(values = tb_node_colors) +
      guides(color = FALSE,
             size = FALSE)

}
