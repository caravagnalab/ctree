make_test_tree = function() {
  data("ctree_input", package = "ctree", envir = environment())

  withr::local_options(ctree.progressBar = FALSE)

  trees = suppressWarnings(suppressMessages(
    ctrees(
      ctree_input$CCF_clusters,
      ctree_input$drivers,
      ctree_input$samples,
      ctree_input$patient,
      ctree_input$sspace.cutoff,
      ctree_input$n.sampling,
      ctree_input$store.max
    )
  ))

  trees[[1]]
}
