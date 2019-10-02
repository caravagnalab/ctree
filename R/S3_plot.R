#' S3 method that plots a REVOLVER tree.
#'
#' @param x A REVOLVER tree (object of class \code{"rev_phylo"}).
#' @param cex Cex for the plot.
#' @param node_palette A function that applied to a number will return a set of colors.
#' By default this is a \code{colorRampPalette} applied to 9 colours of the \code{RColorBrewer}
#' palette \code{Set1}. Colors are generated following a topological sort of the information
#' transfer, which is obtained from \code{igraph}.
#' @param tree_layout A layout that can be used by \code{tidygraph}, which wraps \code{igraph}'s
#' layouts. By default this is a `tree` layout.
#' @param information_transfer If `TRUE`, the plot will show only the information
#' transfer of the tree. The colouring of the nodes of the trees will match the colouring of
#' the drivers.
#' @param icon If `TRUE` the icon tree version of a tree is plot. This type of view does not show
#' a clone unless it has a driver annotated. If this option is true together with \code{information_transfer}
#' an erorr is generated.
#' @param ... Extra parameters
#'
#' @return The plot. If `add_information_transfer = TRUE` the object is a combined figure from
#' package \code{ggpubr}, otherwise it is a single \code{ggplot} object.
#'
#' @export plot.rev_phylo
#'
#' @import crayon
#' @import igraph
#' @import tidygraph
#' @import ggraph
#' @import ggpubr
#' @import RColorBrewer
#'
#' @examples
#' data(CRC.cohort)
#' plot(CRC.cohort$phylogenies[['adenoma_3']][[1]])
plot.ctree = function(x,
                      node_palette = colorRampPalette(RColorBrewer::brewer.pal(n = 9, "Set1")),
                      tree_layout = 'tree',
                      ...)
{
   # Get the tidygraph
    tree = x
    tb_tree = tree$tb_adj_mat
    
    cex = 1
    
    # TODO Color edges as of information transfer
    #  - get path
    #  - modify edges etc.
    # tree$transfer
    
    # Color the nodes by cluster id, using a topological sort
    # to pick the colors in the order of appeareance in the tree
    clones_orderings = igraph::topo_sort(igraph::graph_from_adjacency_matrix(DataFrameToMatrix(tree$transfer$clones)),
                                         mode = 'out')$name
    
    nDrivers = length(clones_orderings) - 1 # avoid GL
    
    drivers_colors = c('white', node_palette(nDrivers))
    names(drivers_colors) = clones_orderings
    
    # Add non-driver nodes, with the same colour
    non_drivers = tb_tree %>%
      activate(nodes) %>%
      filter(!is.driver) %>%
      pull(cluster) # GL is not selected because is NA for is.driver
    
    non_drivers_colors = rep("gainsboro", length(non_drivers))
    names(non_drivers_colors) = non_drivers
    
    tb_node_colors = c(drivers_colors, non_drivers_colors)
    
    # Plot call
    layout <- create_layout(tb_tree, layout = tree_layout)
    
    mainplot = ggraph(layout) +
      geom_edge_link(
        arrow = arrow(length = unit(2 * cex, 'mm')),
        end_cap = circle(5 * cex, 'mm'),
        start_cap  = circle(5 * cex, 'mm')
      ) +
      geom_label_repel(
        aes(
          label = driver,
          x = x,
          y = y,
          colour = cluster
        ),
        na.rm = TRUE,
        nudge_x = .3,
        nudge_y = .3,
        size = 2.5 * cex
      ) +
      geom_node_point(aes(colour = cluster,
                          size = nMuts),
                      na.rm = TRUE) +
      geom_node_text(aes(label = cluster),
                     colour = 'black',
                     vjust = 0.4) +
      coord_cartesian(clip = 'off') +
      # theme_graph(base_size = 8 * cex, base_family = '') +
      theme_void(base_size = 8 * cex) +
      theme(legend.position = 'bottom',
            legend.key.size = unit(3 * cex, "mm")) +
      scale_color_manual(values = tb_node_colors) +
      scale_size(range = c(3, 10) * cex) +
      guides(color = FALSE,
             size = guide_legend("Clone size", nrow = 1)) +
      labs(title = paste(tree$patient),
           subtitle = paste0('Scores ',
                             format(tree$score, scientific = T),
                             '.'))
    
    return(mainplot)
}