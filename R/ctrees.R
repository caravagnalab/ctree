ctrees = function(CCF_clusters,
                              drivers,
                              samples,
                              patient,
                              sspace.cutoff = 10000,
                              n.sampling = 5000,
                              store.max = 100)
{
  
  # TODO - check input formats
  
  pio::pioHdr(
    paste("ctree ~ sample clone trees for", patient),
    toPrint = c(
      `MonteCarlo search if there are more than` = sspace.cutoff,
      `Number of Montecarlo samples, if not exhaustive` = n.sampling,
      `Maximumum number of trees to store, if multiple are available` = store.max
    ),
    prefix = '\t'
  )
  
  # Sample structure for all trees
  structures = trees_sampler(
    CCF_clusters,
    drivers,
    samples,
    patient,
    sspace.cutoff,
    n.sampling,
    store.max
  )
  
  TREES = structures[[1]]
  SCORES = structures[[2]]
  
  # Trees assembly 
  cat("Total number of trees with non-zero scores is", length(TREES))
  
  if(length(TREES) > store.max)
  {
    cat(red(' -- too many, storing top', store.max, '\n'))
    
    TREES = TREES[1:store.max]
    SCORES = SCORES[1:store.max]
    
  } else cat(green(' -- storing all\n'))
  
  LSCORES = as.data.frame(SCORES)
  LSCORES = split(LSCORES, f = LSCORES[,1])
  
  LSCORES = lapply(LSCORES, function(w) w[sample(1:nrow(w), replace = FALSE), , drop = FALSE])
  
  # Shuffle indexes of equal-scoring fits to avoid a sistematic bias in the way I coded this
  permuted.indexes = as.integer(rev(unlist(lapply(LSCORES, rownames))))
  names(permuted.indexes) = NULL
  
  TREES = TREES[permuted.indexes]
  SCORES = SCORES[permuted.indexes]
  
  #####################################################################
  
  easypar::run(
    FUN = function(i)
    {
      ctree(
        CCF_clusters,
        drivers,
        samples,
        patient,
        M = TREES[[i]],
        score = SCORES[i],
        annotation = paste('ctree rank ', i, '/', length(TREES), ' for ', patient, sep ='')
      )
    },
    PARAMS = lapply(seq_along(TREES), list),
    parallel = FALSE
  )
  
  
}
