
#' Create rod-like territories
#'
#' Generates spatial territories with elongated rod-like shapes by creating
#' line segments of varying lengths and assigning cells within specified
#' perpendicular distance ranges. Supports multi-layer territories using
#' morphological operations.
#'
#' @param n_cells Integer specifying the number of cells to simulate
#' @param n_territories Integer specifying the maximum number of rod territories to create
#' @param width_range Numeric vector of length 2 specifying the minimum and maximum
#'   rod width as proportions of total spatial area (default: c(0.0, 0.05))
#' @param length_range Numeric vector of length 2 specifying the minimum and maximum
#'   rod length as proportions of total spatial area (default: c(0.1, 0.5))
#' @param layers Integer specifying the number of morphological layers to create
#'   within each territory (default: 0 for no layers)
#'
#' @return Data frame with columns: barcodes, x, y, Territory
#'
#' @details
#' The function creates rod territories by:
#' 1. Randomly distributing cells in 2D space
#' 2. Selecting random starting points for rods
#' 3. Creating line segments of random length extending from start points
#' 4. Assigning cells within perpendicular distance ranges to form rod shapes
#' 5. Optionally creating layered territories using vesalius morphological operations
#'
#' @examples
#' # Create simple rod territories
#' coords <- rod_map(n_cells = 1000, n_territories = 3,
#'                   width_range = c(0.01, 0.03), length_range = c(0.2, 0.4))
#'
#' # Create layered rod territories
#' coords <- rod_map(n_cells = 1000, n_territories = 2,
#'                   width_range = c(0.02, 0.04), length_range = c(0.3, 0.5),
#'                   layers = 2)
#'
#' @importFrom vesalius build_vesalius_assay generate_tiles territory_morphing layer_territory
rod_map <- function(n_cells = 6000,
    n_territories = 5,
    width_range = c(0.0, 0.05),
    length_range = c(0.1, 0.5),
    layers = 0) {
    if (length(width_range) != 2){
        stop("width_range requires a min width and max width!")
    }
    if (length(length_range) != 2){
        stop("length_range requires a min width and max width!")
    }
    #-------------------------------------------------------------------------#
    # Make rods
    #-------------------------------------------------------------------------#
    width_range <- width_range * 1000
    length_range <- length_range * 1000
    x <- runif(n_cells, min = 1, max = 1000)
    y <- runif(n_cells, min = 1, max = 1000)
    barcodes <- paste0("Cell", seq(1, n_cells))
    coord <- data.frame(barcodes, x, y, "Territory" = 0)
    start <- sample(seq(1, n_cells), size = n_territories, replace = FALSE)
    distances <- lapply(start, rod_it, x, y, width_range, length_range)
    for (i in seq_along(distances)) {
        coord$Territory[distances[[i]]] <- i
    }
    comment(coord) <- "rod"
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
        comment(coord) <- paste0("rod_",layers)
    }
    return(coord)
    
}

#' Create a rod-shaped territory around a starting point
#'
#' This internal function generates a rod-shaped territory by selecting cells
#' within a specified distance range along a randomly chosen line segment.
#'
#' @param idx Integer index of the starting cell position
#' @param x Numeric vector of x-coordinates for all cells
#' @param y Numeric vector of y-coordinates for all cells
#' @param width_range Numeric vector specifying min and max rod width
#' @param length_range Numeric vector specifying min and max rod length
#'
#' @return Logical vector indicating which cells belong to the rod territory
#'
#' @keywords internal
rod_it <- function(idx, x, y, width_range, length_range) {
    d_start <- sqrt(((x - x[idx])^2 + (y - y[idx])^2))
    max_length <- runif(1, min = length_range[1], length_range[2])
    d_end <- sample(which(d_start <= max_length), size = 1)
    x_0 <- x[idx]
    y_0 <- y[idx]
    x_1 <- x[d_end]
    y_1 <- y[d_end]
    d <- dist_line(x_0, y_0, x_1, y_1, x, y)
    width <- runif(2, min = width_range[1], max = width_range[2])
    while (max(width) - min(width) < max(width)/2) {
        width <- runif(2, min = width_range[1], max = width_range[2])
    }
    return(d >= min(width) & d <= max(width))
}


#' Calculate perpendicular distance from points to a line segment
#'
#' This internal function computes the perpendicular distance from multiple points
#' to a line segment defined by two endpoints.
#'
#' @param x_0 Numeric value for x-coordinate of line start point
#' @param y_0 Numeric value for y-coordinate of line start point
#' @param x_1 Numeric value for x-coordinate of line end point
#' @param y_1 Numeric value for y-coordinate of line end point
#' @param x Numeric vector of x-coordinates for points to measure
#' @param y Numeric vector of y-coordinates for points to measure
#'
#' @return Numeric vector of perpendicular distances from each point to the line segment
#'
#' @keywords internal
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
