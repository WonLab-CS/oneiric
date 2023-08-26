
#' Create chaotic territoies
#' @param n_cells int - number of cells to create
#' @param max_expanse numeric -  max proportion of total width to use for
#' territory expansion.
#' @param layers int - number of layers in each rod territory
#' @param chaos character - string specifying which chaos map to use

#' @details Create a a territory using chaos mapping means that there
#' be only one territory but it will have a chaotic shape. 
#' @return coordinate data frame with barcodes, x, y, and Territories
#' @importFrom dplyr %>%
#' @importFrom vesalius build_vesalius_assay generate_tiles territory_morphing layer_territory
chaos_map <- function(n_cells = 6000,
    max_expanse = 0.1,
    layers = 0,
    chaos = "tinkerbell") {
    #-------------------------------------------------------------------------#
    # Make rods
    #-------------------------------------------------------------------------#
    max_width <- max_width * 1000
    max_length <- max_length * 1000
    x <- runif(n_cells, min = 1, max = 1000)
    y <- runif(n_cells, min = 1, max = 1000)
    barcodes <- paste0("cell_", seq(1, n_cells))
    coord <- data.frame(barcodes, x, y, "Territory" = 0)
    coord <- switch(chaos,
        "tinkerbell" = tinkerbell_map(coord))
    
    #-------------------------------------------------------------------------#
    # We will use vesalius if there is a need for layers
    #-------------------------------------------------------------------------#
    ves <- build_vesalius_assay(coordinates = coord, verbose = FALSE) %>%
            generate_tiles(filter_grid = 1, filter_threshold = 1, verbose = FALSE)
    ves@territories <- coord
    ter_to_layer <- unique(coord$Territory[coord$Territory != 0])
    ves <- territory_morphing(ves,
        territory = ter_to_layer,
        trial = "Territory",
        morphology_factor = max_expanse * 1000,
        verbose = FALSE)
    if (layers > 0) {
        ves <- layer_territory(ves,
            territory = ter_to_layer,
            trial = "last",
            layer_depth = layers,
            verbose = FALSE)
        coord$Territory[coord$Territory == ter_to_layer] <-
                paste0(ter_to_layer, "_",
                    ves@territories[coord$Territory ==
                    ter_to_layer, "Layer"])
    } else {
        coord <- ves@territories[, c("barcodes", "x", "y", "Morphology")]
        colnames(coord) <- gsub("Morphology", "Territory", colnames(coord))
    }
    return(coord)
    
}


#' @importFrom RANN nn2
tinkerbell_map <- function(coord, chaos = FALSE) {
    tinker <- tinkerbell(time = 8000)
    tinker$x <- min_max(tinker$x) * max(coord$x)
    tinker$y <- min_max(tinker$y) * max(coord$y)
    knn <- RANN::nn2(coord[, c("x", "y")],
        tinker,
        k = 1)$nn.idx[, 1]
    coord$Territory[knn] <- 1
    return(coord)
}

tinkerbell <- function(time = 6000,
    a = 0.9,
    b = -0.6013,
    c = 2,
    d = 0.5,
    x_0 = -0.72,
    y_0 = -0.64){
    x <- c(x_0, rep(0, time - 1))
    y <- c(y_0, rep(0, time - 1))

    for (i in seq(1, time - 1)) {
        x[i + 1] <- x[i]^2 - y[i]^2 + a * x[i] + b * y[i]
        y[i + 1] <- 2 * x[i] * y[i] + c * x[i] + d * y[i]
    }
    return(data.frame(x, y))
}