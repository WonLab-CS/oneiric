test_that("find_tinkerbell generates parameter sets", {
  # Test basic functionality
  result <- find_tinkerbell(time = 1, n_maps = 5, plot = FALSE, export = FALSE)

  expect_type(result, "matrix")
  expect_equal(ncol(result), 6)  # Should have 6 parameters: a, b, c, d, x_0, y_0
  expect_equal(nrow(result), 5)  # Should have 5 parameter sets
})

test_that("chaos_map creates chaotic territories", {
  # Test chaos map creation
  result <- chaos_map(n_cells = 50, expanse = 0.1, chaos = "tinkerbell", layers = 0)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 50)
  expect_true("barcodes" %in% colnames(result))
  expect_true("x" %in% colnames(result))
  expect_true("y" %in% colnames(result))
  expect_true("Territory" %in% colnames(result))
})
