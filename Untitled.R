require(tidyverse)
require(ggraph)
require(crayon)

library(R.utils)
sourceDirectory('./R_2/')
sourceDirectory('./R')

load('ctree_input.RData')

CCF_clusters = ctree_input$CCF_clusters
drivers = ctree_input$drivers
samples = ctree_input$samples
patient = 'erwerfwer'

ctree =
    CCF_clusters,
    drivers,
    samples,
    patient,
    M,
    score,
    annotation = paste0("Clone tree for patient ", patient))

