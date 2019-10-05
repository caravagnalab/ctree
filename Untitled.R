require(tidyverse)
require(tidygraph)
require(ggraph)
require(crayon)
require(igraph)
require(ggrepel)


library(R.utils)
sourceDirectory('./R_2/')
sourceDirectory('./R')

data('ctree_input')

x = ctrees(
  ctree_input$CCF_clusters,
  ctree_input$drivers,
  ctree_input$samples,
  ctree_input$patient,
  ctree_input$sspace.cutoff,
  ctree_input$n.sampling,
  ctree_input$store.max
)

plot(x[[2]])

plot_information_transfer(x[[1]])
plot_icon(x[[1]])
plot_CCF_clusters(x[[1]])
plot_clone_size(x[[1]])
