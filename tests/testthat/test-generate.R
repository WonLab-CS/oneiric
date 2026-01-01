test_that("generate_sim_data runs without error", {
  # Create temporary directory for testing
  temp_dir <- tempdir()

  # Test that the function runs without error
  # Note: This function has side effects (creates files) so we test it doesn't error
  expect_error(generate_sim_data(output = temp_dir, seed = 42, run_mem = FALSE, simple = FALSE), NA)
})
