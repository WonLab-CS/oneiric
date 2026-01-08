
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
