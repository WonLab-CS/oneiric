test_that("add_interactions works correctly", {
  # Create spatial data
  spatial <- simulate_spatial(n_cells = 50, n_territories = 2, n_samples = 1,
                            pattern = "circle", layers = 0)

  # Add cell labels
  spatial <- simulate_cells(spatial, n_genes = 50, cell_composition = 2)

  # Test add_interactions
  result <- add_interactions(spatial$spatial, k = 3)

  expect_type(result, "list")
  expect_equal(length(result), 1)  # One sample
  expect_true("interactions" %in% colnames(result[[1]]))
})

test_that("add_contact_deg works correctly", {
  # Create spatial data with cells
  spatial <- simulate_spatial(n_cells = 50, n_territories = 2, n_samples = 1,
                            pattern = "circle", layers = 0)
  sim_data <- simulate_cells(spatial, n_genes = 100, cell_composition = 2)

  # Test add_contact_deg
  result <- add_contact_deg(sim_data, k = 3, n_genes = 5)

  expect_type(result, "list")
  expect_true("simulated" %in% names(result))
  expect_true("genes" %in% names(result))
  expect_true("interactions" %in% names(result))

  # Check that genes were added
  expect_true(length(result$genes) > 0)
  expect_true(length(result$interactions) > 0)
})
