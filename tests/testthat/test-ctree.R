test_that("ctrees() generates a ranked list of ctree objects from ctree_input", {
  data("ctree_input")

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

  expect_type(trees, "list")
  expect_true(length(trees) > 0)
  expect_true(all(vapply(trees, inherits, logical(1), what = "ctree")))

  # Trees are ranked by decreasing score
  scores = vapply(trees, function(x) x$score, numeric(1))
  expect_equal(scores, sort(scores, decreasing = TRUE))
})

test_that("ctree() can be constructed manually from an adjacency matrix", {
  data("ctree_input")

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

  x = trees[[1]]

  # Remove the GL node added automatically by ctree()
  M = x$adj_mat
  M = M[rownames(M) != "GL", colnames(M) != "GL", drop = FALSE]

  y = ctree(
    ctree_input$CCF_clusters,
    ctree_input$drivers,
    ctree_input$samples,
    ctree_input$patient,
    M,
    score = 123456,
    annotation = "Manual clone tree"
  )

  expect_s3_class(y, "ctree")
  expect_equal(y$score, 123456)
  expect_equal(y$annotation, "Manual clone tree")
  expect_true(is_tree(y$adj_mat, allow.empty = nrow(y$CCF) == 1))
})

test_that("ctree() rejects input without drivers", {
  data("ctree_input")

  M = matrix(0, ncol = 1, nrow = 1)
  colnames(M) = rownames(M) = "1"

  expect_error(
    ctree(
      ctree_input$CCF_clusters,
      ctree_input$drivers[0, ],
      ctree_input$samples,
      ctree_input$patient,
      M,
      score = 1
    ),
    "No drivers"
  )
})

test_that("ctree() rejects a non-matrix adjacency input", {
  data("ctree_input")

  expect_error(
    ctree(
      ctree_input$CCF_clusters,
      ctree_input$drivers,
      ctree_input$samples,
      ctree_input$patient,
      M = data.frame(from = "1", to = "2"),
      score = 1
    ),
    "adjacency matrix"
  )
})

test_that("print.ctree and summary.ctree run without error", {
  data("ctree_input")

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

  expect_output(print(trees[[1]]), "Tree shape")
  expect_output(summary(trees[[1]]), "Pigeonhole")
})
