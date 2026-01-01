
#' Create chaotic territories using strange attractors
#'
#' Generates spatial territories with chaotic, irregular shapes by mapping
#' trajectories from strange attractors (such as the tinkerbell map) onto
#' spatial coordinates. Creates complex, non-geometric territory boundaries.
#'
#' @param n_cells Integer specifying the number of cells to simulate
#' @param expanse Numeric value specifying the expansion factor for territory
#'   morphology as a proportion of spatial area (default: 0.1)
#' @param chaos Character string specifying which chaotic map to use
#'   (currently only "tinkerbell" is supported)
#' @param layers Integer specifying the number of morphological layers to create
#'   within the territory (default: 0 for no layers)
#'
#' @return Data frame with columns: barcodes, x, y, Territory
#'
#' @details
#' The function creates chaotic territories by:
#' 1. Randomly distributing cells in 2D space
#' 2. Generating a trajectory using a strange attractor (tinkerbell map)
#' 3. Mapping the chaotic trajectory onto spatial coordinates using nearest neighbors
#' 4. Optionally expanding and layering the territory using vesalius morphological operations
#'
#' This creates territories with fractal-like, irregular boundaries that are
#' biologically more realistic than simple geometric shapes.
#'
#' @examples
#' # Create a simple chaotic territory
#' coords <- chaos_map(n_cells = 1000, expanse = 0.2, chaos = "tinkerbell")
#'
#' # Create a layered chaotic territory
#' coords <- chaos_map(n_cells = 1000, expanse = 0.3, chaos = "tinkerbell", layers = 3)
#'
#' @importFrom dplyr %>%
#' @importFrom vesalius build_vesalius_assay generate_tiles territory_morphing layer_territory
chaos_map <- function(n_cells = 6000,
    expanse = 0.1,
    chaos = "tinkerbell",
    layers = 0) {
    #-------------------------------------------------------------------------#
    # Make rods
    #-------------------------------------------------------------------------#
    expanse <- expanse[1] * 1000
    x <- runif(n_cells, min = 1, max = 1000)
    y <- runif(n_cells, min = 1, max = 1000)
    barcodes <- paste0("Cell", seq(1, n_cells))
    coord <- data.frame(barcodes, x, y, "Territory" = 0)
    coord <- switch(chaos,
        "tinkerbell" = tinkerbell_map(coord))
    comment(coord) <- "chaos"
    #-------------------------------------------------------------------------#
    # We will use vesalius if there is a need for layers
    #-------------------------------------------------------------------------#
    if ( expanse > 0 && layers > 0) {
        ves <- build_vesalius_assay(coordinates = coord, verbose = FALSE) %>%
            generate_tiles(filter_grid = 1, filter_threshold = 1, verbose = FALSE)
        ves@territories <- coord
        ter_to_layer <- unique(coord$Territory[coord$Territory != 0])
        ves <- territory_morphing(ves,
            territory = ter_to_layer,
            trial = "Territory",
            morphology_factor = expanse,
            verbose = FALSE)
        ves <- layer_territory(ves,
            territory = ter_to_layer,
            trial = "last",
            layer_depth = layers,
            verbose = FALSE)
        coord$Territory[coord$Territory == ter_to_layer] <-
                paste0(ter_to_layer, "_",
                    ves@territories[coord$Territory ==
                    ter_to_layer, "Layer"])
        comment(coord) <- paste0("chaos_",layers)
        return(coord)
    } else if (expanse > 0 && layers == 0) {
        ves <- build_vesalius_assay(coordinates = coord, verbose = FALSE) %>%
            generate_tiles(filter_grid = 1, filter_threshold = 1, verbose = FALSE)
        ves@territories <- coord
        ter_to_layer <- unique(coord$Territory[coord$Territory != 0])
        ves <- territory_morphing(ves,
            territory = ter_to_layer,
            trial = "Territory",
            morphology_factor = expanse,
            verbose = FALSE)
        coord <- ves@territories[, c("barcodes", "x", "y", "Morphology")]
        colnames(coord) <- gsub("Morphology", "Territory", colnames(coord))
        
        return(coord) 
    } else {
        return(coord)
    }
}


#' @importFrom RANN nn2

#' Apply tinkerbell chaotic map to spatial coordinates
#'
#' This internal function generates a chaotic territory pattern using the tinkerbell
#' strange attractor and maps it onto spatial coordinates.
#'
#' @param coord Data frame containing spatial coordinates (barcodes, x, y, Territory)
#' @param chaos Logical flag (unused parameter, kept for compatibility)
#'
#' @return Data frame with Territory column updated to reflect chaotic pattern
#'
#' @keywords internal
tinkerbell_map <- function(coord, chaos = FALSE) {
    data(oneiric)
    params <- map_params[sample(seq(1,nrow(map_params)),1), ]
    tinker <- tinkerbell(time = 8000,
        a = params[1],
        b = params[2],
        c = params[3],
        d = params[4],
        x_0 = params[5],
        y_0 = params[6])
    tinker$x <- min_max(tinker$x) * max(coord$x)
    tinker$y <- min_max(tinker$y) * max(coord$y)
    knn <- RANN::nn2(coord[, c("x", "y")],
        tinker,
        k = 1)$nn.idx[, 1]
    coord$Territory[knn] <- 1
    return(coord)
}

#' Generate tinkerbell strange attractor trajectory
#'
#' This internal function simulates the tinkerbell strange attractor by iteratively
#' applying the tinkerbell map equations to generate a chaotic trajectory.
#'
#' @param time Integer specifying the number of iterations to perform
#' @param a Numeric parameter for the tinkerbell map equation
#' @param b Numeric parameter for the tinkerbell map equation
#' @param c Numeric parameter for the tinkerbell map equation
#' @param d Numeric parameter for the tinkerbell map equation
#' @param x_0 Numeric initial x-coordinate
#' @param y_0 Numeric initial y-coordinate
#'
#' @return Data frame with 'x' and 'y' columns containing the trajectory coordinates
#'
#' @keywords internal
tinkerbell <- function(time = 8000,
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

#' @export
find_tinkerbell <- function(time = 120,
    n_maps = 96,
    plot = FALSE,
    export = FALSE,
    file_name = "maps") {
    p <- 1
    maps <- vector("list", 96)
    start <- Sys.time()

    if(plot) {
        file <- paste0(file_name, ".pdf")
        pdf(file, width = 20, height = 20)
        par(mfrow = c(4, 4))
    }
    while(p <= n_maps){
        if(as.numeric(difftime(Sys.time(), start, units = "secs")) >= time) {
            break
        }
        a <- runif(1, -5, 5)
        b <- runif(1, -5, 5)
        c <- runif(1, -5, 5)
        d <- runif(1, -5, 5)
        x_0 <- runif(1, -5, 5)
        y_0 <- runif(1, -5, 5)
        tinker <- oneiric:::tinkerbell(time = 6000,
            a = a,b=b,c=c,d=d,x_0=x_0,y_0=y_0)
        if (any(is.na(tinker$x)) || any(is.na(tinker$y))) {
            next
        }

        if (max(tinker$x) > 10 || max(tinker$y) > 10 || min(tinker$x) < (-10) || min(tinker$y) < (-10)){
            next
        }
        maps[[p]] <- c(a,b,c,d,x_0,y_0)
        names(maps[[p]]) <- c("a","b","c","d","x_0","y_0")
        if (plot) {
            plot(tinker$x,tinker$y,pch = 19, col = "cadetblue", main = paste("Sample p = ",p))
        }
        p <- p + 1
    }
    if(plot) {dev.off()}
    maps <- do.call("rbind", maps)
    if(export){
        file <- paste0(file_name, ".csv")
        write.csv(maps, file = file, row.names = FALSE)
    }
    return(maps)
}

