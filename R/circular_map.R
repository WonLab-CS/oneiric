#' Create circular  territoies
#' @param n_cells int - number of cells to create
#' @param n_territories int -  max number of territories to generate
#' @param layers int - number of layers in each rod territory
#' @param expanse numeric -  max proportion of total width to use as 
#' circle radius
#' @return coordinate data frame with barcodes, x, y, and Territories 
#' @importFrom vesalius build_vesalius_assay generate_tiles territory_morphing layer_territory
circular_map <- function(n_cells = 6000,
    n_territories = 5,
    expanse = c(0.01, 0.5),
    layers = 0) {
    if (length(expanse) != 2) {
        stop("Expanse need 2 values: min expanse and max expanse!")
    }
    x <- runif(n_cells, min = 1, max = 1000)
    y <- runif(n_cells, min = 1, max = 1000)
    barcodes <- paste0("Cell", seq(1, n_cells))
    coord <- data.frame(barcodes, x, y, "Territory" = 0)
    
    centers <- sample(seq(1, n_cells), size = n_territories, replace = FALSE)
    distances <- lapply(centers, function(idx, x, y, min_expanse, max_expanse) {
        d <- sqrt(((x - x[idx])^2 + (y - y[idx])^2))
        max_expanse <- runif(1, min = min_expanse, max = max_expanse)
        return(d <= max_expanse)
    }, x = x, y = y, min_expanse = expanse[1] * 1000, max_expanse = expanse[2] * 1000)

    for (i in seq_along(distances)) {
        coord$Territory[distances[[i]]] <- i
    }
    #-------------------------------------------------------------------------#
    # We will use vesalius if there is a need for layers
    #-------------------------------------------------------------------------#
    
    if (layers > 0) {
        ves <- build_vesalius_assay(coordinates = coord, verbose = FALSE) %>%
            generate_tiles(filter_grid = 1, filter_threshold = 1, verbose = FALSE)
        ves@territories <- coord
        ter_to_layer <- unique(coord$Territory[coord$Territory != 0])
        for (i in seq_along(ter_to_layer)) {
            ves <- layer_territory(ves,
                territory = ter_to_layer[i],
                trial = "Territory",
                layer_depth = layers,
                verbose = FALSE)
            coord$Territory[coord$Territory == ter_to_layer[i]] <-
                paste0(ter_to_layer[i], "_",
                    ves@territories[coord$Territory ==
                    ter_to_layer[i],ncol(ves@territories)])
        }

    }
    return(coord)
}