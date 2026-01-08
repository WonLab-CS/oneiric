#' Generate simulated data for analysis
#'
#' This function generates simulated spatial transcriptomics datasets for
#' benchmarking spatial mapping tools. It creates multiple simulation scenarios
#' with different cell type compositions and territory patterns.
#'
#' @param output Character string specifying the path to the output directory
#' @param seed Integer seed for consistent data generation
#' @param run_mem Logical indicating whether to run memory-intensive computational benchmarks
#' @param simple Logical indicating whether to generate simplified test datasets
#'
#' @details
#' This function generates various simulated datasets including:
#' - Random territory patterns with one or two cell types per territory
#' - Circular territories with varying cell numbers (when run_mem = TRUE)
#' - Layered territories with differential expression patterns
#'
#' The function is primarily designed to reproduce the simulation scenarios
#' used in the original manuscript.
#'
#' @export
generate_sim_data <- function(output,
    seed = 1453,
    run_mem = FALSE,
    simple = FALSE){
    #-------------------------------------------------------------------------#
    # Fix parameters internally
    #-------------------------------------------------------------------------#
    set.seed(seed)
    #-------------------------------------------------------------------------#
    # run spatial sims with random territory types
    #-------------------------------------------------------------------------#
    random <- simulate_spatial(n_cells = 5000,
        n_territories = 3,
        n_samples = 12,
        pattern = "random",
        expanse = c(0.15, 0.3),
        width_range = c(0, 0.3),
        length_range = c(0.15,0.3),
        layers = seq(0,4),
        force_cells = 5)
    #-------------------------------------------------------------------------#
    # Add only one 1 cell per territory
    #-------------------------------------------------------------------------#
    random_one_cell <- simulate_cells(random,
        de_prob = 0.5,
        de_layer = 0.05,
        cell_composition = 1,
        no_label = FALSE,
        randomize_cells = TRUE)
    
    random_one_cell$spatial <- add_interactions(random_one_cell$spatial, k = 6)
    export_simulation(spatial = random_one_cell$spatial,
        cells = random_one_cell$counts,
        out_dir = output,
        file_tag = "one_cell")
    #-------------------------------------------------------------------------#
    # Add 2 cell types per territory
    #-------------------------------------------------------------------------#
    random_two_cell <- simulate_cells(random,
        de_prob = 0.5,
        de_layer = 0.05,
        cell_composition = 2,
        no_label = FALSE,
        randomize_cells = TRUE)
    random_two_cell$spatial <- add_interactions(random_two_cell$spatial, k = 6)
    export_simulation(spatial = random_two_cell$spatial,
        cells = random_two_cell$counts,
        out_dir = output,
        file_tag = "two_cell")
   
    #-------------------------------------------------------------------------#
    # Create a simple circle territory set but with varying number of cells
    # Create thsi data sets with interactions for consitency 
    #-------------------------------------------------------------------------#
    if(run_mem) {
        n_cells <- c(100, 500, 1000, 2500, 5000, 10000, 20000, 50000, 100000, 200000)
        for (i in seq_along(n_cells)) {
            circle <- simulate_spatial(n_cells = n_cells[i], 
                n_territories = 5, 
                n_samples = 3, 
                pattern = "circle", 
                expanse = c(0.1, 0.3), 
                force_cells = 5)
            circle <- simulate_cells(circle,
                cell_composition = 2,
                no_label = TRUE)
            circle$spatial <- add_interactions(circle$spatial)
            export_simulation(spatial = circle$spatial,
                cells = circle$counts,
                out_dir = output,
                file_tag = paste0("computational_performance_", n_cells[i], "_cells"))
        }
    }
    #-------------------------------------------------------------------------#
    # Get old set up as well just in case but with new code
    #-------------------------------------------------------------------------#
    if (simple){
        circular <- simulate_spatial(n_cells = 5000, 
            n_territories = 5, 
            n_samples = 12, 
            pattern = "circle", 
            expanse = c(0.1, 0.25), 
            force_cells = 10)
        circular <- simulate_cells(circular,
            cell_composition = 2,
            no_label = TRUE)
        circular$spatial <- add_interactions(circular$spatial, k = 6)
        export_simulation(spatial = circular$spatial,
            cells = circular$counts,
            out_dir = output,
            file_tag = "circle")


        circular_layer <- simulate_spatial(n_cells = 5000,
            n_territories = 1,
            n_samples = 12,
            pattern = "circle",
            expanse = c(0.4, 0.5),
            layers = 5,
            force_cells = 12)
        circular_layer <- simulate_cells(circular_layer,
            de_prob = 0.5,
            de_layer = 0.05,
            cell_composition = 2,
            no_label = TRUE)
        circular_layer$spatial <- add_interactions(circular_layer$spatial, k = 6)
        export_simulation(spatial = circular_layer$spatial,
            cells = circular_layer$counts,
            out_dir = output,
            file_tag = "layered")


        circular_dropped <- simulate_spatial(n_cells = 5000,
            n_territories = 2,
            n_samples = 12,
            pattern = "circle",
            expanse = c(0.15, 0.3),
            layers = 4,
            force_cells = 5)
        circular_dropped <- simulate_cells(circular_dropped,
            de_prob = 0.5,
            de_layer = 0.05,
            cell_composition = 2,
            no_label = TRUE,
            randomize_cells = TRUE)
        circular_dropped$spatial <- add_interactions(circular_dropped$spatial, k = 6)
        export_simulation(spatial = circular_dropped$spatial,
            cells = circular_dropped$counts,
            out_dir = output,
            file_tag = "dropped")
    }
}