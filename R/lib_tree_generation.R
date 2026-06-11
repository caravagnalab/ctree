hashTrees = function(clonevol.obj, sample.groups)
{
  trees = map.trees = NULL
  for(s in sample.groups)
  {
    l = lapply(
      clonevol.obj$models[[s]],
      function(x, region)
      {
        # Nothing to hash here -- empty tree
        if(nrow(x) == 0) return(NULL)

        # hash ;)
        model = DataFrameToMatrix(x)
        model = sort(MatrixToEdges(model))

        return(list(model = paste(model, collapse = ":"), region = region))
      },
      region = s)

    trees = rbind(trees, Reduce(rbind, l))
  }

  rownames(trees) = NULL
  trees = data.frame(
    model = unlist(trees[, 'model']),
    model = unlist(trees[, 'region'])
  )
  colnames(trees) = c('model', 'region')

  trees = trees[trees$model != "", , drop = FALSE]

  trees$model = as.factor(trees$model)

  map.trees = split(trees, f = trees$model)
  map.trees = lapply(map.trees, function(x)paste(x[, 'region'], collapse = ','))

  return(map.trees[map.trees != ""])

}



consensusModel = function(clonevol.obj, sample.groups)
{
  S = lapply(clonevol.obj$models, function(x) Reduce(rbind, x))
  S = lapply(S, unique)
  S = Reduce(rbind, S)

  counts = DataFrameToEdges(S)
  counts = table(counts)

  counts = cbind(edgesToDataFrame(names(counts)), unlist(counts))

  counts$counts = NULL

  counts = split(counts, f = counts$to)
  counts = lapply(counts, function(x) {y = x$Freq/ sum(x$Freq); names(y) = x$from; return(y)} )

  S = DataFrameToMatrix(S)
  S = MatrixToDataFrame(S)

  return(list(S = S, weights = counts))
}


all.possible.trees = function(
  G,
  W,
  sspace.cutoff = 10000,
  n.sampling = 1000
  )
{
  M = DataFrameToMatrix(G)
  r = root(M)
  parents = sapply(colnames(M), parent_of, model = M)

  parents = parents[!unlist(lapply(parents, is.null))]

  # singletons -- template
  singletons = parents[unlist(lapply(parents, function(x) length(x) == 1))]
  nsingletons = length(singletons)

  sglt = expand.grid(singletons, stringsAsFactors = FALSE)
  if(ncol(sglt) == 0) {
    sglt = NULL
  }

  alternatives = parents[unlist(lapply(parents, function(x) length(x) > 1))]
  nalternatives = length(alternatives)

  altn = NULL
  if(nalternatives > 0)
  {
    combalternatives = prod(unlist(lapply(alternatives, length)))

    ex_search = ifelse(
      combalternatives < sspace.cutoff,
      "exahustive" %>% crayon::bold(),
      paste0('Monte Carlo for ', n.sampling, 'distinct trees') %>% crayon::bold()
    )

    cli::cli_alert_info(
      "Total {.field {combalternatives}} tree structures - search is {.count {ex_search}}")

    if(combalternatives > sspace.cutoff)
    {
      return(
        weighted.sampling(
          DataFrameToMatrix(G),
          W,
          n.sampling
        )
      )
    }
    altn = expand.grid(alternatives, stringsAsFactors = FALSE)
  }
  else cat(red('There are no alternatives!\n'))

  # all combinations
  if(is.null(altn) && is.null(sglt)) stop('Error -- no trees?')

  if(is.null(sglt) && !is.null(altn)) comb = altn
  if(!is.null(sglt) && is.null(altn)) comb = sglt
  if(!is.null(sglt) && !is.null(altn)) comb =  cbind(altn, sglt)

  progress_bar = getOption('ctree.progressBar', default = TRUE)
  if (progress_bar) cli::cli_progress_bar("Generating trees", total = nrow(comb))

  models = NULL
  for(i in 1:nrow(comb))
  {
    if (progress_bar)
      cli::cli_progress_update()

    tree = data.frame(from = unlist(comb[i, ]), to = colnames(comb), stringsAsFactors = FALSE)
    test.tree = DataFrameToMatrix(tree)

    if(length(root(test.tree)) > 1 ) next;

    if(!igraph::is_dag(igraph::graph_from_adjacency_matrix(test.tree))) next;

    models = append(models, list(tree))
  }

  return(models)
}


rankTrees = function(TREES, MI.table, structural.score)
{
  progress_bar = getOption('ctree.progressBar', default = TRUE)
  if (progress_bar) cli::cli_progress_bar("Ranking trees", total = length(TREES))

  MI.TREES = NULL
  for(i in 1:length(TREES))
  {
    if (progress_bar)
      cli::cli_progress_update()

    M = DataFrameToMatrix(TREES[[i]])
    M = M[colnames(MI.table), colnames(MI.table)]

    M.entries = MI.table[which(M > 0, arr.ind = TRUE)]

    val = NA
    if(all(is.null(structural.score))) val = prod(M.entries)
    else val = prod(M.entries) * structural.score[i]

    if(any(M.entries == 0))
    {
      n = sum(M.entries == 0)
      M.entries[M.entries == 0] = 1e-9

      warning("Used MI correction for", n, "entries equal 0; set equal to 1e-9.")
    }

    MI.TREES = c(MI.TREES, val)
  }

  o = order(MI.TREES, decreasing = TRUE)
  MI.TREES = MI.TREES[o]
  TREES = TREES[o]
  structural.score = structural.score[o]

  return(list(TREES = TREES, SCORES = MI.TREES, PENALTIES = structural.score))
}


computeMI.table = function(binary.data, MI.Bayesian.prior = 0, add.control = FALSE)
{
  if(add.control) binary.data = rbind(binary.data, wt = 0)

  # • a=0:maximum likelihood estimator (see entropy.empirical)
  # • a=1/2:Jeffreys’ prior; Krichevsky-Trovimov (1991) entropy estimator
  # • a=1:Laplace’s prior
  # • a=1/length(y):Schurmann-Grassberger (1996) entropy estimator
  # • a=sqrt(sum(y))/length(y):minimax prior

  MI.table = matrix(
    apply(
      expand.grid(colnames(binary.data), colnames(binary.data)),
      1,
      function(x) {
        # Counting process.. with a Bayesian prior
        i = x[1]
        j = x[2]
        jo11 = (binary.data[, i] %*% binary.data[, j])/nrow(binary.data)
        jo10 = (binary.data[, i] %*% (1-binary.data[, j]))/nrow(binary.data)
        jo01 = ((1-binary.data[, i]) %*% binary.data[, j])/nrow(binary.data)
        jo00 = 1 - (jo10 + jo01 + jo11)
        entropy::mi.Dirichlet(matrix(c(jo11, jo10, jo01, jo00), nrow = 2), a = MI.Bayesian.prior)
      }),
    byrow = TRUE, ncol = ncol(binary.data))
  colnames(MI.table) = rownames(MI.table) = colnames(binary.data)

  return(MI.table)
}


weightMI.byMultinomial = function(MI.table, W)
{
  Coeff = matrix(0, ncol = ncol(MI.table), nrow = nrow(MI.table))
  rownames(Coeff) = rownames(MI.table)
  colnames(Coeff) = colnames(MI.table)

  for(j in names(W))
    Coeff[ names(W[[j]]) , j] = unlist(W[[j]])

  return(matrixcalc::hadamard.prod(MI.table, Coeff))
}

weighted.sampling = function(G, W, n)
{
  S = S.hashcodes = NULL

  sampleT = function()
  {
    r = root(G)

    E = setdiff(colnames(G), r)
    E = sample(E, length(E))

    tree = NULL
    repeat {
      parent_draw = sapply(E, function(node){
        parents = W[[node]]
        sample(names(parents), prob = unlist(parents), size = 1)
      })

      tree = data.frame(from = unlist(parent_draw), to = names(parent_draw), stringsAsFactors = FALSE)
      R = c(r, reach(tree, r))

      if(all(colnames(G) %in% R)) break
    }

    return(tree)
  }

  c = 0
  repeat{
    Tree = sampleT()
    hash = paste(sort(DataFrameToEdges(Tree)), collapse = ':')

    if(!(hash %in% S.hashcodes))
    {
      S = append(S, list(Tree))
      S.hashcodes = c(S.hashcodes, hash)
      c = c + 1
    }

    if(c == n) break;
  }

  return(S)

}


# for every node  x --> y1 ... yK, the number of times that the CCF of x is greater than the sum of the CCFs of y1 ... yK
node.penalty.for.branching = function(TREES, ccf)
{
  progress_bar = getOption('ctree.progressBar', default = TRUE)
  if (progress_bar) cli::cli_progress_bar("Computing branching penalties", total = length(TREES))

  scores = NULL
  for(x in 1:length(TREES))
  {
    if (progress_bar) cli::cli_progress_update()

    t = DataFrameToMatrix(TREES[[x]])
    nodes = rownames(ccf)

    c = sapply(nodes, function(n) {
      cl = children(t, n)
      if(length(cl) == 0) return(1)

      1 - sum(as.numeric(ccf[n, ] < colSums(ccf[cl, , drop = FALSE])))/ncol(ccf)
    })

    scores = c(scores, prod(c))
  }

   return(scores)
}
