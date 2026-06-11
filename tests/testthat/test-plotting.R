test_that("plot.ctree returns a ggplot object", {
  x = make_test_tree()

  p = plot(x)

  expect_s3_class(p, "ggplot")
})

test_that("plot_CCF_clusters returns a ggplot object", {
  x = make_test_tree()

  p = plot_CCF_clusters(x)

  expect_s3_class(p, "ggplot")
})

test_that("plot_clone_size returns a ggplot object", {
  x = make_test_tree()

  p = plot_clone_size(x)

  expect_s3_class(p, "ggplot")
})

test_that("plot_information_transfer returns a ggplot object", {
  x = make_test_tree()

  p = plot_information_transfer(x)

  expect_s3_class(p, "ggplot")
})

test_that("plot_icon returns a ggplot object", {
  x = make_test_tree()

  p = plot_icon(x)

  expect_s3_class(p, "ggplot")
})
