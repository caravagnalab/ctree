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

CCF_clusters = ctree_input$CCF_clusters
drivers = ctree_input$drivers
samples = ctree_input$samples
patient = 'erwerfwer'

sspace.cutoff = 10000
n.sampling = 5000
store.max = 100

ctree =
    CCF_clusters,
    drivers,
    samples,
    patient,
    M,
    score,
    annotation = paste0("Clone tree for patient ", patient))


x = ctree(CCF_clusters,
      drivers,
      samples,
      patient,
      M = TREES[[1]],
      score = 2)

print(x)
plot(x)

x = ctrees(
  ctree_input$CCF_clusters,
  ctree_input$drivers,
  ctree_input$samples,
  ctree_input$patient,
  ctree_input$sspace.cutoff,
  ctree_input$n.sampling,
  ctree_input$store.max
)
