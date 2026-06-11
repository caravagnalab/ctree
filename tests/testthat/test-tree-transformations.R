test_that("edgesToDataFrame and DataFrameToEdges round-trip", {
  edges = c("A~B", "A~C", "B~D")

  df = edgesToDataFrame(edges)
  expect_equal(df$from, c("A", "A", "B"))
  expect_equal(df$to, c("B", "C", "D"))

  expect_equal(DataFrameToEdges(df), edges)
})

test_that("edgesToDataFrame returns an empty data frame for no edges", {
  df = edgesToDataFrame(character(0))

  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 0)
})

test_that("edgesToMatrix and MatrixToEdges round-trip", {
  edges = c("A~B", "A~C")

  M = edgesToMatrix(edges)
  expect_equal(M["A", "B"], 1)
  expect_equal(M["A", "C"], 1)
  expect_equal(M["B", "C"], 0)

  expect_equal(sort(MatrixToEdges(M)), sort(edges))
})

test_that("DataFrameToMatrix and MatrixToDataFrame round-trip", {
  df = data.frame(
    from = c("A", "A"),
    to = c("B", "C"),
    stringsAsFactors = FALSE
  )

  M = DataFrameToMatrix(df)
  expect_true(all(c("A", "B", "C") %in% colnames(M)))
  expect_equal(M["A", "B"], 1)
  expect_equal(M["A", "C"], 1)

  df2 = MatrixToDataFrame(M)
  expect_equal(nrow(df2), nrow(df))
  expect_setequal(paste(df2$from, df2$to), paste(df$from, df$to))
})
