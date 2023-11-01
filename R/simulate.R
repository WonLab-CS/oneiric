#-----------------------------------------------------------------------------#
################################# ONEIRIC #####################################
#-----------------------------------------------------------------------------#

#' Create spatial territories
#' @param n_cells int - number of cells to include in siluated data
#' @param n_territories int - number of territories to create 
#' @param n_samples int - number of samples to generate
#' @param pattern character - string defining the type of territories to 
#' generate ("circle", "rod", "chaos")
#' @param layers int - number of layers in a territory
#' @param max_expanse numeric - expansion proportion for circular and chaotic
#' territories as a proportion of max size (1000 "pixels")
#' @param max_width numeric - max width of rod territories as a proportion of
#' max size (1000 "pixels")
#' @param max_length numeric - max length of rod territories as a proportion of
#' max size (1000 "pixels")
#' @param border logical - should border cells have specific gene expression?
#' @export 
simulate_spatial <- function(n_cells = 6000,
    n_territories = 5,
    n_samples = 10,
    pattern = "circle",
    max_expanse = 0.3,
    max_width = 0.3,
    max_length = 0.5,
    border = TRUE,
    layers = 0) {

    sample_coord <- vector("list", n_samples)
    for (i in seq_len(n_samples)) {
        tmp_coord <- switch(pattern,
            "circle" = circular_map(n_cells, n_territories, max_expanse, layers),
            "rod" = rod_map(n_cells, n_territories, max_width, max_length, layers),
            "chaos" = chaos_map(n_cells, max_expanse, "tinkerbell", layers))
        tmp_coord$sample <- i
        sample_coord[[i]] <- tmp_coord
    }
    return(sample_coord)
}

#' @importFrom splatter newSplatParams setParams splatSimulateGroups
#' @export 
simulate_cells <- function(spatial,
    cell_composition = 1,
    n_genes = 2000,
    seed = 1729) {
    if (!is(spatial,"list")){
        spatial <- list(spatial)
    }
    spatial <- lapply(spatial, function(spatial,cell_composition,n_genes,seed){
        n_ters <- table(spatial$Territory)
        n_cells <- unlist(lapply(n_ters, function(cell, cell_comp){
                return(rep(ceiling(cell / cell_comp), times = cell_comp))
        }, cell_composition))
        cell_types <- make.unique(rep(names(n_ters), each = cell_composition))
        params <- newSplatParams(batchCells = nrow(spatial), nGenes = n_genes)
        params <- setParam(params, "seed", seed)
        sim <- splatSimulateGroups(params, group.prob = n_cells / sum(n_cells), verbose = FALSE)
        sim <- counts(sim)
        return(sim)
    }, cell_composition = cell_composition,
    n_genes = n_genes,
    seed = seed)
    return(spatial)
}
