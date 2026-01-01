test_that("export_simulation creates files", {
  # Create test data
  spatial <- simulate_spatial(n_cells = 20, n_territories = 2, n_samples = 1,
                            pattern = "circle", layers = 0)
  sim_data <- simulate_cells(spatial, n_genes = 50, cell_composition = 1)

  # Create temporary directory for testing
  temp_dir <- tempdir()

  # Test export_simulation
  expect_error(export_simulation(sim_data$spatial, sim_data$counts,
                               out_dir = temp_dir, file_tag = "test"), NA)

  # Check that files were created
  expected_files <- c(
    "test_spatial_coordinates_sample_1.csv",
    "test_gene_counts_sample_1.csv"
  )

  for (file in expected_files) {
    expect_true(file.exists(file.path(temp_dir, file)))
  }
})

test_that("min_max normalization works", {
  # Test with varying values
  x <- c(1, 5, 10, 3, 8)
  result <- min_max(x)

  expect_length(result, 5)
  expect_equal(min(result), 0)
  expect_equal(max(result), 1)

  # Test with constant values
  x_const <- c(5, 5, 5, 5)
  result_const <- min_max(x_const)

  expect_length(result_const, 4)
  expect_true(all(result_const == x_const))
})
