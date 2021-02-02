trees_sampler = function(CCF_clusters,
                         drivers,
                         samples,
                         patient,
                         sspace.cutoff = 10000,
                         n.sampling = 5000,
                         store.max = 100)
{
  TREES = SCORES = NULL
  
  # We will need to use the CCF clusters for this patient
  clusters = CCF_clusters
  nclusters = nrow(clusters)
  
  clonal.cluster = clusters %>% filter(is.clonal) %>% pull(cluster)
  
  df_clusters = data.frame(clusters[, samples],
                           row.names = clusters$cluster)
  
  print(clusters)
  cat('\n')
  
  if (nclusters == 1)
  {
    cli::cli_alert_warning('Model with 1 node, trivial trees returned')
    
    M = matrix(0, ncol = 1, nrow = 1)
    colnames(M) = rownames(M) = rownames(clusters)
    
    TREES = append(TREES, list(M))
    SCORES = c(SCORES, 1)
  }
  else
  {
    # ################## Generate all trees that are compatible with the observed CCFs, we do this
    # ################## by analyzing one sample at a time.
    alternative = NULL
    alternative$models = ClonEvol_surrogate(clusters, samples, clonal.cluster, min.CCF = 0.01)
    clonevol.obj = alternative
    
    # remove the trees which have no edges (returned for samples with only 1 cluster for instance)
    numSol = sapply(clonevol.obj$models, function(w) {
      sum(sapply(w, nrow) > 0)
    })
    
    # pio::pioStr(
    #   " Trees per region", 
    #   paste(numSol, collapse = ', '),
    #   prefix = crayon::green(clisymbols::symbol$tick),
    #   suffix = '\n')
    
    cli::cli_alert_success(
      "Trees per region {.field {paste(numSol, collapse = ', ')}}"
    )
    
    ################## Build all possible clonal trees
    # 1) hash them
    # 2) create a consensus as the union of all trees
    # 3) generate or sample a large number of possible trees, where a parent x --> y is assigned
    #    with probability proportional to how often the edge is detected
    
    CLONAL.TREES = hashTrees(clonevol.obj, samples)
    
    CONSENSUS = consensusModel(clonevol.obj, samples)
    CONSENSUS.TREE = CONSENSUS$S
    WEIGHTS.CONSENSUS.TREE = CONSENSUS$weights
    
    # print(CONSENSUS.TREE)
    # print(WEIGHTS.CONSENSUS.TREE)
    
    # # Sampling is carried out if there are more than 'sspace.cutoff' trees, in that case we
    # # sample 'n.sampling' possible trees. Otherwise all possible trees are generated.
    TREES = all.possible.trees(
      CONSENSUS.TREE,
      WEIGHTS.CONSENSUS.TREE,
      sspace.cutoff,
      n.sampling
    )
    
    # ################## Ranking trees. A tree is good according to the following factors:
    # # 1) the MI among the variables x and y, if they are connected by an edge x --> y [TODO: consider if we really need MI]
    # # 2) the Multinomial probability of edge x --> y in the trees determined by the CCF
    # # 3) for every edge  x --> y, the number of times that the CCF of x is greater than the CCF of y
    # # 3) for every node  x --> y1 ... yK, the number of times that the CCF of x is greater than the sum of the CCFs of y1 ... yK
    # binary.data = revolver:::binarize(x$dataset, samples)
    
    # 1) MI from binarized data -- different options, with a control sample which avoids 0log(0)
    # • a=0:maximum likelihood estimator (see entropy.empirical)
    # • a=1/2:Jeffreys’ prior; Krichevsky-Trovimov (1991) entropy estimator
    # • a=1:Laplace’s prior
    # • a=1/length(y):Schurmann-Grassberger (1996) entropy estimator
    # • a=sqrt(sum(y))/length(y):minimax prior
    binary.data = t(df_clusters)
    binary.data[binary.data > 0] = 1
    
    MI.table = computeMI.table(binary.data,
                               MI.Bayesian.prior = 0,
                               add.control = TRUE)
    # Old parameter now fixed to FALSE
    use.MI = FALSE
    if (!use.MI)
      MI.table[TRUE] = 1
    
    
    # Steps 1 and 2 are collapsed, multiply MI by the Multinomial probability
    MI.table = weightMI.byMultinomial(MI.table, WEIGHTS.CONSENSUS.TREE)
    
    # 3) Get penalty for direction given CCFs -- this is done for all possible edges in the data
    # CCF = clusters[, samples, drop = FALSE]
    # penalty.CCF.direction = edge.penalty.for.direction(TREES, CCF)
    penalty.CCF.direction = 1
    
    # 4) Compute the branching penalty  --  this is done for each tree that we are considering
    # pio::pioStr(
    #   " Pigeonhole Principle",
    #   prefix = crayon::green(clisymbols::symbol$tick),
    #   suffix = '\n')
    
    penalty.CCF.branching = node.penalty.for.branching(TREES, df_clusters)
    
    cli::cli_h3("\nRanking trees")
    # pio::pioStr(
    #   " Ranking trees",
    #   prefix = crayon::green(clisymbols::symbol$tick),
    #   suffix = '\n')
    
    RANKED = rankTrees(TREES, MI.table, penalty.CCF.branching)
    TREES = RANKED$TREES
    SCORES = RANKED$SCORES
    
    TREES = TREES[SCORES > 0]
    SCORES = SCORES[SCORES > 0]
    
    TREES = lapply(TREES, DataFrameToMatrix)
  }
  
  return(list(adj_mat = TREES, scores = SCORES))
}
