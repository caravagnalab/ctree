test_that("root finds the node with no parents", {
  # A -> B, A -> C
  M = matrix(c(0, 1, 1,
                0, 0, 0,
                0, 0, 0), nrow = 3, byrow = TRUE)
  colnames(M) = rownames(M) = c("A", "B", "C")

  expect_equal(root(M), "A")
})

test_that("children and parent_of are consistent", {
  # A -> B, A -> C
  M = matrix(c(0, 1, 1,
                0, 0, 0,
                0, 0, 0), nrow = 3, byrow = TRUE)
  colnames(M) = rownames(M) = c("A", "B", "C")

  expect_equal(children(M, "A"), c("B", "C"))
  expect_length(children(M, "B"), 0)

  expect_equal(parent_of(M, "B"), "A")
  expect_equal(parent_of(M, "C"), "A")
  expect_null(parent_of(M, "A"))
})

test_that("leaves finds nodes with no children", {
  M = matrix(c(0, 1, 1,
                0, 0, 0,
                0, 0, 0), nrow = 3, byrow = TRUE)
  colnames(M) = rownames(M) = c("A", "B", "C")

  expect_equal(leaves(M), c("B", "C"))
})

test_that("is_tree validates a proper tree structure", {
  M = matrix(c(0, 1, 1,
                0, 0, 0,
                0, 0, 0), nrow = 3, byrow = TRUE)
  colnames(M) = rownames(M) = c("A", "B", "C")

  expect_true(is_tree(M))

  # Two roots is not a tree
  M2 = matrix(0, ncol = 2, nrow = 2)
  colnames(M2) = rownames(M2) = c("A", "B")
  expect_false(is_tree(M2))

  # An empty matrix is only a tree if explicitly allowed
  expect_false(is_tree(M2, allow.empty = FALSE))
  expect_true(is_tree(M2, allow.empty = TRUE))
})

test_that("reach computes the transitive closure of a node", {
  df = data.frame(
    from = c("A", "A", "B"),
    to = c("B", "C", "D"),
    stringsAsFactors = FALSE
  )

  expect_equal(sort(reach(df, "A")), c("B", "C", "D"))
  expect_equal(reach(df, "B"), "D")
  expect_null(reach(df, "D"))
})
