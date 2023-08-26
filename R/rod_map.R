
#' Create Rod like territoies
#' @param n_cells int - number of cells to create
#' @param n_territories int -  max number of territories to generate
#' @param layers int - number of layers in each rod territory
#' @param max_width numeric -  max proportion of total width to use as 
#' rod width
#' @param max_length numeric - max proportion of total length to use as
#' rod length
#' @return coordinate data frame with barcodes, x, y, and Territories
#' @importFrom vesalius build_vesalius_assay generate_tiles territory_morphing layer_territory
rod_map <- function(n_cells = 6000,
    n_territories = 5,
    layers = 0,
    max_width = 0.05,
    max_length = 0.5) {
    #-------------------------------------------------------------------------#
    # Make rods
    #-------------------------------------------------------------------------#
    max_width <- max_width * 1000
    max_length <- max_length * 1000
    x <- runif(n_cells, min = 1, max = 1000)
    y <- runif(n_cells, min = 1, max = 1000)
    barcodes <- paste0("cell_", seq(1, n_cells))
    coord <- data.frame(barcodes, x, y, "Territory" = 0)
    start <- sample(seq(1, n_cells), size = n_territories, replace = FALSE)
    distances <- lapply(start, rod_it, x, y, max_width, max_length)
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

rod_it <- function(idx, x, y, max_width, max_length) {
    d_start <- sqrt(((x - x[idx])^2 + (y - y[idx])^2))
    max_length <- runif(1, min = min(d_start), max_length)
    d_end <- sample(which(d_start <= max_length), size = 1)
    x_0 <- x[idx]
    y_0 <- y[idx]
    x_1 <- x[d_end]
    y_1 <- y[d_end]
    d <- dist_line(x_0, y_0, x_1, y_1, x, y)
    return(d <= max_width)
}


dist_line <- function(x_0, y_0, x_1, y_1, x, y) {
    px <- x_1 - x_0
    py <- y_1 - y_0
    norm <- px * px + py * py
    u <- ((x - x_0) * px + (y - y_0) * py) / norm
    u[u > 1] <- 1
    u[u < 0] <- 0
    x_2 <- x_0 + u * px
    y_2 <- y_0 + u * py
    dx <- x_2 - x
    dy <- y_2 - y
    dist <- sqrt(dx * dx + dy * dy)
    return(dist)
}
