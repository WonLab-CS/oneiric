
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
    chaos = "tinkerbell",
    layers = 0) {
    #-------------------------------------------------------------------------#
    # Make rods
    #-------------------------------------------------------------------------#
    max_expanse <- max_expanse * 1000
    x <- runif(n_cells, min = 1, max = 1000)
    y <- runif(n_cells, min = 1, max = 1000)
    barcodes <- paste0("cell_", seq(1, n_cells))
    coord <- data.frame(barcodes, x, y, "Territory" = 0)
    coord <- switch(chaos,
        "tinkerbell" = tinkerbell_map(coord),
        "bogdanov" = bogdanov_map(coord))
    
    #-------------------------------------------------------------------------#
    # We will use vesalius if there is a need for layers
    #-------------------------------------------------------------------------#
    if ( max_expanse > 0 && layers > 0) {
        ves <- build_vesalius_assay(coordinates = coord, verbose = FALSE) %>%
            generate_tiles(filter_grid = 1, filter_threshold = 1, verbose = FALSE)
        ves@territories <- coord
        ter_to_layer <- unique(coord$Territory[coord$Territory != 0])
        ves <- territory_morphing(ves,
            territory = ter_to_layer,
            trial = "Territory",
            morphology_factor = max_expanse * 1000,
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
        return(coord)
    } else if (max_expanse > 0 && layers == 0) {
        ves <- build_vesalius_assay(coordinates = coord, verbose = FALSE) %>%
            generate_tiles(filter_grid = 1, filter_threshold = 1, verbose = FALSE)
        ves@territories <- coord
        ter_to_layer <- unique(coord$Territory[coord$Territory != 0])
        ves <- territory_morphing(ves,
            territory = ter_to_layer,
            trial = "Territory",
            morphology_factor = max_expanse * 1000,
            verbose = FALSE)
        coord <- ves@territories[, c("barcodes", "x", "y", "Morphology")]
        colnames(coord) <- gsub("Morphology", "Territory", colnames(coord))
        return(coord) 
    } else {
        return(coord)
    }
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

bogdanov_map <- function(coord, chaos = FALSE) {
    bog <- bogdanov(time = 8000)
    browser()
    bog$x <- min_max(bog$x) * max(coord$x)
    bog$y <- min_max(bog$y) * max(coord$y)
    knn <- RANN::nn2(coord[, c("x", "y")],
        bog,
        k = 1)$nn.idx[, 1]
    coord$Territory[knn] <- 1
    return(coord)
}

bogdanov <- function(time,
    e = 0,
    k = 1.2,
    m = 0){
    # Initial start might not be optimal here 
    x <- rep(0.3, length(time))
    y <- rep(0.3, length(time))
    for (i in seq(1, time - 1)){
        y[i + 1] <- y[i] + e * y[i] + k * x[i] * (x[i] - 1) + m * x[i] * y[i]
        x[i + 1] <- x[i] + y[i + 1]
    }
    return(data.frame(x, y))
}


# Bogdanov map function in R
dimx <- 1024    # Image dimension
dimy <- 1024

I <- integer(dimx * dimy)
for (i in 1:(dimx * dimy)) I[i] <- 0

e <- 0          # birth and grow
k <- 1.2        # Discretization
m <- 0          # stability
step <- 0.03    # Initial point step

for (oy in seq(-1.0, 1.0, by = step)) {
  for (ox in seq(-0.7, 1.3, by = step)) {
    if (ox == 0 || oy == 0) next     # Avoid origin
    x <- ox
    y <- oy
    for (q in 1:100000) {
      y <- y + e * y + k * x * (x - 1) + m * x * y
      x <- x + y
      if (abs(x) > 2 || abs(y) > 2) break    # the map diverges
      imx <- (x - 0.3) * dimx / 2 + dimx / 2
      imy <- y * dimy / 2 + dimy / 2
      if (abs(x) > 1e-6 || abs(y) > 1e-6) {    # Avoid origin
        ix <- as.integer(imx + 0.5)
        iy <- as.integer(imy + 0.5)
        if (ix >= 0 && ix < dimx && iy >= 0 && iy < dimy)
          I[dimx * (dimy - 1 - iy) + ix + 1] <- I[dimx * (dimy - 1 - iy) + ix + 1] + 1
      }
    }
  }
}

maxd <- max(I)    # Compute max frequency

fo <- file("bogdanov_map.pgm", "wb")
writeBin(charToRaw(paste("P5\n", dimx, " ", dimy, "\n255\n", sep = "")), fo)
for (i in 1:(dimx * dimy)) {
  v <- ifelse(I[i] == 0, 0, 32 + (256 - 32) * I[i] / maxd + 0.5)
  v <- min(255, v)
  b <- as.raw(255 - v)
  writeBin(b, fo)
}
close(fo)
