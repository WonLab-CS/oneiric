#' Create circular territories
#'
#' Generates spatial territories with circular shapes by randomly placing centers
#' and assigning cells within specified radius ranges. Supports multi-layer territories
#' using the vesalius package for morphological processing.
#'
#' @param n_cells Integer specifying the number of cells to simulate
#' @param n_territories Integer specifying the maximum number of circular territories to create
#' @param expanse Numeric vector of length 2 specifying the minimum and maximum radius
#'   as proportions of the total spatial area (default: c(0.01, 0.5))
#' @param layers Integer specifying the number of morphological layers to create
#'   within each territory (default: 0 for no layers)
#'
#' @return Data frame with columns: barcodes, x, y, Territory
#'
#' @details
#' The function creates circular territories by:
#' 1. Randomly distributing cells in 2D space
#' 2. Selecting random center points for territories
#' 3. Assigning cells to territories based on distance from centers
#' 4. Optionally creating layered territories using vesalius morphological operations
#'
#' @examples
#' # Create simple circular territories
#' coords <- circular_map(n_cells = 1000, n_territories = 3, expanse = c(0.1, 0.3))
#'
#' # Create layered circular territories
#' coords <- circular_map(n_cells = 1000, n_territories = 2, expanse = c(0.2, 0.4), layers = 3)
#'
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
    comment(coord) <- "circle"
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
        comment(coord) <- paste0("circle_",layers)
    }
    
    return(coord)
}