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
    expanse = c(0.01, 0.5),
    width_range = c(0, 0.3),
    length_range = c(0.1,0.5),
    border = TRUE,
    layers = 0) {

    sample_coord <- vector("list", n_samples)
    for (i in seq_len(n_samples)) {
        tmp_coord <- switch(pattern,
            "circle" = circular_map(n_cells, n_territories, expanse, layers),
            "rod" = rod_map(n_cells, n_territories, width_range, length_range, layers),
            "chaos" = chaos_map(n_cells, expanse, "tinkerbell", layers))
        tmp_coord$sample <- i
        tmp_coord$barcodes <- paste0(tmp_coord$barcodes, "_",i)
        sample_coord[[i]] <- tmp_coord
    }
    return(sample_coord)
}

#' @importFrom splatter newSplatParams setParams splatSimulateGroups
#' @export 
simulate_cells <- function(spatial,
    cell_composition = 1,
    n_genes = 2000,
    as_layer = FALSE,
    de_prob = 0.25,
    de_layer = 0.05,
    seed = 1729) {
    if (!is(spatial,"list")){
        spatial <- list(spatial)
    }
    
    spatial <- do.call("rbind", spatial)
    ters <- sort(unique(spatial$Territory))
    total_cells <- nrow(spatial)
    for (i in seq_along(ters)){
        bars <- spatial$barcodes[spatial$Territory == ters[i]]
        split_size <- ceiling(length(bars) / cell_composition)
        split_names <- make.unique(as.character(rep(ters[i], times = cell_composition)))
        for (j in seq_len(cell_composition)){
            samp <- sample(bars, size = min(c(split_size,length(bars))))
            spatial$Territory[spatial$barcodes %in% samp] <-
                split_names[j]
            bars <- bars[!bars %in% samp]
        }
    }

    n_cells <- table(spatial$Territory)
    params <- splatter::newSplatParams(batchCells = total_cells, nGenes = n_genes)
    params <- splatter::setParam(params, "seed", seed)
    if (as_layer) {
        de <- c(de_prob, de_prob, rep(de_layer, times = length(grep("_", names(n_cells)))))
    } else {
        de <- rep(de_prob, times = length(n_cells))
    }
    sim <- splatter::splatSimulateGroups(params,
        group.prob = as.numeric(n_cells / sum(n_cells)),
        de.prob = de,
        verbose = FALSE)
    spatial <- assign_barcodes(spatial, sim)
    counts <- retrieve_counts(spatial, sim)
    spatial <- split(spatial, spatial$sample)
    return(list("counts" = counts, "spatial" = spatial))
}

