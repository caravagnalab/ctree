



#' Summary of a \code{"ctree"} object.
#' @description 
#' 
#' Reports some summary statistics for a \code{"rev_phylo"} object, 
#' which is a bit more than just using \code{print.}
#'
#' @param object A \code{ctree} tree.
#' @param ... Extra parameters
#'
#' @return None.
#' @export summary.ctree
#' @exportS3Method summary ctree
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
#' summary(x[[1]])
summary.ctree <- function(object, ...) {
  print.ctree(object, ...)
  
  pio::pioStr("CCF clusters:",
              "",
              prefix = '\n',
              suffix = '\n\n')
  print(object$CCF)
  
  pio::pioStr("Drivers:",
              "",
              prefix = '\n',
              suffix = '\n\n')  
  print(object$drivers)
  
  stat = stats.rev_phylo(object)
  
  pio::pioStr("Pigeonhole principle:",
              crayon::green(
                nrow(stat$CCF.pigeonhole) * ncol(stat$CCF.pigeonhole) - stat$violations['pp']
              ),
              crayon::red(stat$violations['pp']),
              prefix = '\n',
              suffix = '\n\n')  
  print.noquote(stat$CCF.pigeonhole)
  
  pio::pioStr("Goodness-of-fit:",
              stat$gofit,
              "",
              prefix = '\n',
              suffix = '\n\n')  
}




stats.rev_phylo = function(x)
{
  
  M = MatrixToDataFrame(x$adj_mat)
  M = M[M$from != 'GL', , drop = FALSE]
  
  if (nrow(M) == 0) {
    return(
      list(
        probs = NULL,
        violations = NULL,
        gofit = 1,
        Suppes = NULL,
        CCF.crossing.rule = NULL,
        CCF.pigeonhole = NULL
      )
    )
    
  }
  
  # CCF data for this patient
  CCF = x$CCF[, x$samples, drop = FALSE]
  
  ## Branching penalty is the Pigeonhole Principle
  Mmatrix = DataFrameToMatrix(M)
  
  branching.penalty = sapply(x$CCF %>% pull(cluster),
                             function(n) {
                               # print(n)
                               cl = children(Mmatrix, n)
                               if (length(cl) == 0)
                                 return(NULL)
                               CCF[n, , drop = FALSE] >= colSums(CCF[cl, , drop = FALSE])
                             })
  branching.penalty = branching.penalty[!unlist(lapply(branching.penalty, is.null))]
  branchp = Reduce(rbind, branching.penalty)
  rownames(branchp) = names(branching.penalty)
  
  # cr.viol = sum(as.integer(!direction.penalty))
  ph.viol = sum(as.integer(!branchp))
  
  # violations = c(tp = tp.viol, pr = pr.viol, cr = cr.viol, pp = ph.viol)
  violations = c(pp = ph.viol)
  
  # Some alternative Goodness-Of-Fit measures. In the end, we take the
  # exact definition of the penalty for phylogenetic trees.
  # gofit = 2 * (x$numNodes - 1) * x$numRegions + x$numNodes  * x$numRegions
  # gofit = 1 - sum(violations)/gofit
  # gofit = nrow(branchp) * x$numRegions
  # gofit = 1 - ph.viol/gofit
  capture.output({
    gofit = node.penalty.for.branching(list(M), CCF)
  })
  
  
  return(
    list(
      violations = violations,
      gofit = gofit,
      CCF.pigeonhole = branchp
    )
  )
  
}



#
# #' S3 method calling \code{\link{revolver_plt_tree}}.
# #'
# #' @param x An object of class \code{"rev_phylo"}.
# #' @param file Output file, or \code{NA}.
# #' @param palette RColorBrewer palette to colour clusters.
# #' @param cex Cex for the graph.
# #' @param alpha Transparency.
# #' @param verbose Output.
# #' @param ... Extra parameters
# #'
# #' @return Nothing
# #' @export plot.rev_phylo
# #' @import crayon
# #'
# #' @examples
# #' data(CRC.cohort)
# #' plot(CRC.cohort$phylogenies[['adenoma_3']][[1]])
# plot.rev_phylo = function(x,
#                           file = NA,
#                           palette = 'Set1',
#                           cex = 1,
#                           alpha = 0.7)
# {
#   revolver_plt_tree(
#     x,
#     file = file,
#     # edge.width = edge.width,
#     # edge.label = edge.label,
#     palette = palette,
#     cex = cex,
#     alpha = alpha
#   )
#   invisible(NULL)
# }

# # a-la-dictionary
# driver = function(x, c) {
#   x$dataset[which(x$dataset$cluster == c &
#                     x$dataset$is.driver), 'variantID']
# }
# driver.indexOf = function(x, d) {
#   x$dataset[which(x$dataset$variantID == d &
#                     x$dataset$is.driver), 'cluster']
# }

# # Compute the frontier of "var" in a model. The frontier is the set of
# # mutations in Gamma that are directly reachable via
# # the transitive closure of the --> relation. So, these are the events
# # selected by "var"'s evolutionary trajectories.
# information.transfer = function(x,
#                                 transitive.closure = FALSE,
#                                 indistinguishable = FALSE)
# {
#   aux = function(r)
#   {
#     r.d = driver(x, r)
#
#     if (length(r.d) > 0)
#       return(r) # stop recursion
#
#     c = children(model, r)
#
#     if (is.null(c))
#       return(NULL) # leaf
#
#     # recursion, reduction etc.
#     return(Reduce(union,
#                   sapply(c, aux)))
#   }
#
#   model = x$adj_mat
#   nodes.drivers = x$dataset %>%
#     filter(is.driver) %>%
#     pull(cluster)
#
#   # Root
#   df = expand.grid(
#     from = x$root,
#     to = aux(x$root),
#     stringsAsFactors = FALSE
#   )
#
#   # and all other stuff
#   for (n in nodes.drivers)
#   {
#     df.n = NULL
#
#     for (c in children(model, n))
#       df.n = rbind(df.n,
#                    expand.grid(
#                      from = n,
#                      to = aux(c),
#                      stringsAsFactors = FALSE
#                    ))
#
#     df = rbind(df, df.n)
#   }
#
#   # Then, for all clones, we expand the drivers that they have
#   expanded = apply(df,
#                    1,
#                    function(w) {
#                      e.from = driver(x, w['from'])
#                      if (length(e.from) == 0)
#                        e.from = w['from']
#
#                      e.to = driver(x, w['to'])
#
#                      expand.grid(from = e.from,
#                                  to = e.to,
#                                  stringsAsFactors = FALSE)
#                    })
#   expanded = Reduce(rbind, expanded)
#
#   # Then we see if we want to expand also the indistinguishible
#   if (indistinguishable) {
#     for (n in nodes.drivers) {
#       expanded = rbind(expanded,
#                        expand.grid(
#                          from = driver(x, n),
#                          to = driver(x, n),
#                          stringsAsFactors = FALSE
#                        ))
#     }
#   }
#
#   # If we need transitive closures, we compute them here
#   if (transitive.closure)
#   {
#     df = DataFrameToMatrix(df)
#     df = nem::transitive.closure(df, mat = T, loops = FALSE)
#
#     expanded = DataFrameToMatrix(expanded)
#     expanded = nem::transitive.closure(expanded, mat = T, loops = FALSE)
#
#     df = MatrixToDataFrame(df)
#     expanded = MatrixToDataFrame(expanded)
#   }
#
#   return(list(clones = df, drivers = expanded))
# }


