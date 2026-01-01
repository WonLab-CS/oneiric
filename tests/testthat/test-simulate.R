test_that("simulate_spatial creates expected output", {
  # Test basic functionality
  result <- simulate_spatial(n_cells = 100, n_territories = 2, n_samples = 1,
                           pattern = "circle", layers = 0)

  expect_type(result, "list")
  expect_length(result, 1)
  expect_true("barcodes" %in% colnames(result[[1]]))
  expect_true("x" %in% colnames(result[[1]]))
  expect_true("y" %in% colnames(result[[1]]))
  expect_true("Territory" %in% colnames(result[[1]]))
  expect_true("sample" %in% colnames(result[[1]]))

  # Check that we have the expected number of cells
  expect_equal(nrow(result[[1]]), 100)
})

test_that("simulate_cells creates expected output", {
  # Create spatial data first
  spatial <- simulate_spatial(n_cells = 50, n_territories = 2, n_samples = 1,
                            pattern = "circle", layers = 0)

  # Test simulate_cells
  result <- simulate_cells(spatial, n_genes = 100, cell_composition = 1)

  expect_type(result, "list")
  expect_true("counts" %in% names(result))
  expect_true("spatial" %in% names(result))

  # Check counts structure
  expect_type(result$counts, "list")
  expect_equal(length(result$counts), 1)  # One sample
  expect_equal(nrow(result$counts[[1]]), 100)  # n_genes
  expect_equal(ncol(result$counts[[1]]), 50)   # n_cells

  # Check spatial structure
  expect_type(result$spatial, "list")
  expect_equal(length(result$spatial), 1)  # One sample
  expect_true("cell_labels" %in% colnames(result$spatial[[1]]))
})
